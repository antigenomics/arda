"""Multi-node (SLURM) sharding for very large inputs.

Single-node runs already stream in bounded memory (see ``annotate.mapper``). For
cluster scale we split the input once into N shards, annotate each shard as an
independent SLURM array task, then concatenate the per-shard AIRR TSVs:

    arda cluster split-fasta  big.fastq work/shards --shards 50
    # SLURM array task i:
    arda annotate -i work/shards/shard_<i>.fasta -o work/out/out_<i>.tsv ...
    arda cluster merge  work/out  big.airr.tsv

``arda slurm`` writes (and optionally submits) a single ``submit.sh`` that chains
all three with SLURM job dependencies — see ``render_submit_script``.
"""

from __future__ import annotations

import gzip
import shutil
from itertools import islice
from pathlib import Path

from .annotate import io as seqio

__all__ = ["split", "split_pairs", "merge", "render_submit_script",
           "render_rnaseq_submit_script", "SHARD_GLOB", "RNASEQ_SHARD_GLOB"]

SHARD_GLOB = "shard_*.fasta"
RNASEQ_SHARD_GLOB = "shard_*.airr.tsv"

# 5 digits, so `sorted()` on the names is numeric order. Plain `shard_10` sorts BEFORE
# `shard_2`, and the RNA-seq merge depends on shard order being the original read order.
_SHARD_WIDTH = 5
_MAX_SHARDS = 10 ** _SHARD_WIDTH - 1


def split(input: str | Path, out_dir: str | Path, shards: int,
          *, prefix: str = "shard") -> list[Path]:
    """Round-robin split a FASTA/FASTQ into ``shards`` FASTA files (one pass).

    Round-robin (record ``k`` → shard ``k % shards``) balances load even when
    record sizes vary, and every record lands in exactly one shard.
    """
    if shards < 1:
        raise ValueError("shards must be >= 1")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / f"{prefix}_{i}.fasta" for i in range(shards)]
    handles = [open(p, "w") for p in paths]
    try:
        for k, (sid, seq) in enumerate(seqio.read_sequences(input)):
            handles[k % shards].write(f">{sid}\n{seq}\n")
    finally:
        for h in handles:
            h.close()
    return paths


def _open_bytes(path: Path):
    """Open a (possibly gzipped) file in BINARY mode.

    Binary matters: the shard writer copies FASTQ records through verbatim. Text mode would
    decode and re-encode, and newline translation could rewrite a quality string.
    """
    return gzip.open(path, "rb") if path.suffix == ".gz" else open(path, "rb")


def _count_records(path: Path) -> int:
    """FASTQ records in *path*, by counting newlines in 1 MiB blocks."""
    n = 0
    with _open_bytes(path) as fh:
        while True:
            block = fh.read(1 << 20)
            if not block:
                break
            n += block.count(b"\n")
    return n // 4


def split_pairs(r1: str | Path, out_dir: str | Path, *, shards: int,
                r2: str | Path | None = None,
                prefix: str = "shard") -> list[tuple[Path, Path | None]]:
    """Split FASTQ into ``shards`` CONTIGUOUS blocks of read *pairs*, byte for byte.

    Not :func:`split`. That one writes FASTA — dropping the quality string ``merge_pair``'s
    per-base tie-break needs under ``--reconstruct`` — and round-robins *records*, which puts a
    fragment's two mates in different shards. Mate separation is not a hypothetical defect
    here: it produced a published false discovery in this project's own data (a spurious
    "R2-only blind spot") that had to be retracted.

    **Contiguous, not round-robin**, and that is load-bearing. Concatenating the per-shard
    Stage-1 AIRR in shard order then reproduces the single-node row order *exactly*, so Stage 2
    and Stage 3 see byte-identical input and the sharded result is byte-identical to a
    single-node run. Round-robin would only give a permutation, and the clonotype fold is not
    permutation-invariant (`correct` collapses error children onto the parent it meets first,
    and coverage assignment is first-with-longest-overlap-wins).

    R1 and R2 are cut at the same record boundaries, so mate *k* always lands in the same shard
    as mate *k*; no read ids are parsed to achieve it. ``map`` re-checks the pairing on every
    shard anyway.

    Args:
        r1: FASTQ (optionally gzipped). Single-end if ``r2`` is None.
        out_dir: written as ``<prefix>_00000_R1.fastq`` (+ ``_R2`` when paired).
        shards: number of blocks; a shard receiving no records is not written.
        r2: the mate file.
        prefix: shard file stem.

    Returns:
        ``(r1_path, r2_path | None)`` per non-empty shard, in shard order.

    Raises:
        ValueError: on FASTA input, a bad shard count, an empty input, or mates of
            different lengths (a truncated R2, caught here rather than after N wasted tasks).
    """
    r1, out_dir = Path(r1), Path(out_dir)
    r2 = Path(r2) if r2 is not None else None
    if not 1 <= shards <= _MAX_SHARDS:
        raise ValueError(f"shards must be between 1 and {_MAX_SHARDS}, got {shards}")
    for f in (r1, r2):
        if f is not None and seqio.detect_format(f) != "fastq":
            raise ValueError(
                f"{f} is not FASTQ. `arda cluster split` shards reads with their quality "
                f"strings; use `arda cluster split-fasta` for the FASTA/amplicon path."
            )

    total = _count_records(r1)
    if total == 0:
        raise ValueError(f"no FASTQ records in {r1}")
    if r2 is not None:
        n2 = _count_records(r2)
        if n2 != total:
            raise ValueError(
                f"R1 and R2 differ in length ({total} vs {n2} records); one file is truncated."
            )

    block = -(-total // shards)  # ceil, so `shards` blocks cover everything
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[tuple[Path, Path | None]] = []
    srcs = [(r1, "R1"), *([(r2, "R2")] if r2 is not None else [])]
    handles = [_open_bytes(p) for p, _ in srcs]
    try:
        for i in range(shards):
            paths: list[Path] = []
            counts: list[int] = []
            for (_, tag), fh in zip(srcs, handles):
                dest = out_dir / f"{prefix}_{i:0{_SHARD_WIDTH}d}_{tag}.fastq"
                got = 0
                with open(dest, "wb") as out:
                    remaining = block
                    while remaining:
                        take = min(remaining, 65536)
                        lines = list(islice(fh, take * 4))
                        if not lines:
                            break
                        out.writelines(lines)
                        got += len(lines) // 4
                        remaining -= take
                if got == 0:
                    dest.unlink()          # never leave an empty shard: `detect_format` raises on one
                    continue
                paths.append(dest)
                counts.append(got)
            if not paths:
                break
            if len(set(counts)) != 1:
                raise ValueError(
                    f"R1 and R2 disagree at shard {i} ({counts}); one file is truncated."
                )
            written.append((paths[0], paths[1] if len(paths) > 1 else None))
    finally:
        for fh in handles:
            fh.close()
    if not written:
        raise ValueError(f"no shards written from {r1}")
    return written


def merge(shard_outputs: str | Path | list[Path], output: str | Path) -> Path:
    """Concatenate per-shard AIRR TSVs into one, keeping a single header.

    ``shard_outputs`` may be a directory (its ``*.tsv`` are merged in sorted
    order) or an explicit list of files.
    """
    output = Path(output)
    if isinstance(shard_outputs, (str, Path)) and Path(shard_outputs).is_dir():
        files = sorted(Path(shard_outputs).glob("*.tsv"))
    else:
        files = [Path(p) for p in shard_outputs]  # type: ignore[union-attr]
    if not files:
        raise FileNotFoundError("no shard outputs to merge")
    with open(output, "w") as out:
        for i, p in enumerate(files):
            with open(p) as fh:
                header = fh.readline()
                if i == 0:
                    out.write(header)
                shutil.copyfileobj(fh, out)
    return output


def render_submit_script(
    input: str | Path,
    output: str | Path,
    work_dir: str | Path,
    *,
    shards: int,
    organism: str = "human",
    seqtype: str = "nt",
    threads: int = 8,
    strand: str = "both",
    map_d: bool = True,
    partition: str | None = None,
    time: str = "04:00:00",
    mem: str = "8G",
    arda_mmseqs: str | None = None,
) -> str:
    """Render a ``submit.sh`` that chains split → array-annotate → merge on SLURM.

    Uses ``sbatch --array`` with ``--wrap`` and an ``afterok`` dependency so the
    merge runs only once every shard succeeds. ``arda_mmseqs`` (if given) is
    exported so array tasks find the binary.
    """
    input, output, work_dir = Path(input), Path(output), Path(work_dir)
    shards_dir = work_dir / "shards"
    out_dir = work_dir / "out"
    last = shards - 1
    sbatch_common = [f"--cpus-per-task={threads}", f"--time={time}", f"--mem={mem}"]
    if partition:
        sbatch_common.append(f"--partition={partition}")
    common = " ".join(sbatch_common)
    map_d_flag = "--map-d" if map_d else "--no-map-d"
    env_line = f'export ARDA_MMSEQS="{arda_mmseqs}"' if arda_mmseqs else "# ARDA_MMSEQS from environment"
    return f"""#!/usr/bin/env bash
# Generated by `arda slurm`. Annotate {input} across {shards} SLURM array tasks.
set -euo pipefail
{env_line}

mkdir -p "{shards_dir}" "{out_dir}"

# 1. Split once (cheap, single pass).
arda cluster split-fasta "{input}" "{shards_dir}" --shards {shards}

# 2. Annotate each shard as an array task.
ARRAY_JID=$(sbatch --parsable {common} --job-name=arda-map --array=0-{last} \\
  --output="{out_dir}/slurm-%A_%a.log" \\
  --wrap 'arda annotate -i "{shards_dir}/shard_${{SLURM_ARRAY_TASK_ID}}.fasta" \\
          -o "{out_dir}/out_${{SLURM_ARRAY_TASK_ID}}.tsv" \\
          --organism {organism} --seqtype {seqtype} --threads {threads} \\
          --strand {strand} {map_d_flag}')

# 3. Merge once all shards succeed.
sbatch {common} --job-name=arda-merge --dependency=afterok:$ARRAY_JID \\
  --output="{out_dir}/slurm-merge.log" \\
  --wrap 'arda cluster merge "{out_dir}" "{output}"'
"""


def render_rnaseq_submit_script(
    r1: str | Path,
    out_prefix: str,
    work_dir: str | Path,
    *,
    shards: int,
    r2: str | Path | None = None,
    out_dir: str | Path = ".",
    organism: str = "human",
    threads: int = 8,
    kmer: int = 12,
    min_score: float = 75.0,
    reconstruct: bool = False,
    assemble: bool = True,
    complete_only: bool = True,
    map_d: bool = True,
    partition: str | None = None,
    time: str = "04:00:00",
    mem: str = "8G",
    reduce_time: str = "08:00:00",
    reduce_mem: str = "16G",
    arda_mmseqs: str | None = None,
) -> str:
    """Render a ``submit.sh`` chaining split → array-``map`` → reduce for paired RNA-seq.

    A sibling of :func:`render_submit_script` rather than a generalisation of it. That one
    chains split → array-``annotate`` → merge over a single FASTA, and its last step is a pure
    concatenation. This chain differs in every step: the shard unit is a read *pair*, the files
    are FASTQ with quality, and the last step is a **reduce** — merge, then assemble and correct
    once over the whole merged AIRR. Folding both into one renderer would mean a parameter
    deciding which of two unrelated pipelines you get.

    Only Stage 1 is distributed. `correct` collapses error variants and counts distinct
    fragments globally, and `assemble` grows contigs across reads, so sharding either would
    double-count clones and silently drop exactly the long-CDR3 contigs Stage 3 exists to build.

    Two details in the array body that are not cosmetic:

    * ``printf "%05d"`` — shard names are zero-padded so ``sorted()`` is numeric. The merge
      concatenates in name order and that order *is* read order.
    * ``[ -s "$f" ] || exit 0`` — a shard with no reads must not fail its task, or the
      ``afterok`` dependency drops the whole reduce step. (``split_pairs`` does not write empty
      shards; this is the belt to its braces, for a resubmitted or hand-edited array range.)
    """
    r1, work_dir, out_dir = Path(r1), Path(work_dir), Path(out_dir)
    shards_dir, map_dir = work_dir / "shards", work_dir / "out"
    last = shards - 1
    common = " ".join([f"--cpus-per-task={threads}", f"--time={time}", f"--mem={mem}"]
                      + ([f"--partition={partition}"] if partition else []))
    reduce_common = " ".join([f"--cpus-per-task={threads}", f"--time={reduce_time}",
                              f"--mem={reduce_mem}"]
                             + ([f"--partition={partition}"] if partition else []))
    env_line = f'export ARDA_MMSEQS="{arda_mmseqs}"' if arda_mmseqs else "# ARDA_MMSEQS from environment"
    r2_split = f' --r2 "{r2}"' if r2 else ""
    r2_map = f' --r2 "{shards_dir}/shard_${{i}}_R2.fastq"' if r2 else ""
    flags = f"--organism {organism} --threads {threads} --kmer {kmer} --min-score {min_score}"
    if reconstruct:
        flags += " --reconstruct"
    if not map_d:
        flags += " --no-map-d"
    reduce_flags = f"--organism {organism} --threads {threads}"
    reduce_flags += " --assemble" if assemble else " --no-assemble"
    reduce_flags += " --complete-only" if complete_only else " --all-junctions"
    if not map_d:
        reduce_flags += " --no-map-d"
    return f"""#!/usr/bin/env bash
# Generated by `arda cluster submit`.
#
# Stage 1 (map) is sharded across an array; Stages 2-3 (assemble, correct) are GLOBAL and run
# ONCE, in the reduce step, over the merged Stage-1 AIRR. Sharding them would double-count
# clonotypes and never build the cross-read contigs Stage 3 exists for.
#
# Shards are CONTIGUOUS blocks of read pairs, so concatenating the per-shard AIRR in shard
# order reproduces the single-node row order exactly and the result is byte-identical to a
# single-node run.
set -euo pipefail
{env_line}

mkdir -p "{shards_dir}" "{map_dir}" "{out_dir}"

# 1. Split once into contiguous pair blocks (quality preserved; mates never separated).
SPLIT_JID=$(sbatch --parsable {common} --job-name=arda-split \\
  --output="{map_dir}/slurm-split.log" \\
  --wrap 'arda cluster split --r1 "{r1}"{r2_split} --out-dir "{shards_dir}" --shards {shards}')

# 2. Stage 1 per shard.
ARRAY_JID=$(sbatch --parsable {common} --job-name=arda-map --array=0-{last} \\
  --dependency=afterok:$SPLIT_JID \\
  --output="{map_dir}/slurm-%A_%a.log" \\
  --wrap 'i=$(printf "%05d" ${{SLURM_ARRAY_TASK_ID}}); \\
          f="{shards_dir}/shard_${{i}}_R1.fastq"; [ -s "$f" ] || exit 0; \\
          arda map --r1 "$f"{r2_map} \\
            -o "{map_dir}/shard_${{i}}.airr.tsv" \\
            --report "{map_dir}/shard_${{i}}.map.json" {flags}')

# 3. Reduce ONCE: merge -> assemble -> correct.
sbatch {reduce_common} --job-name=arda-reduce --dependency=afterok:$ARRAY_JID \\
  --output="{map_dir}/slurm-reduce.log" \\
  --wrap 'arda cluster reduce --shard-dir "{map_dir}" --out-dir "{out_dir}" \\
            --out-prefix {out_prefix} {reduce_flags}'
"""
