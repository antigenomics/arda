"""Stages 2-3 as one function, shared by the single-node and sharded delivery paths.

The RNA-seq pipeline is map -> assemble -> correct. Stage 1 (`map`) is per-read and shards
perfectly; Stages 2 and 3 are **global** and do not shard at all:

* `correct` collapses sequencing-error variants onto a parent clonotype and counts distinct
  fragments. Run per shard, one clone split across N shards is counted N times, and error
  children collapse against a fraction of their parent's depth.
* `assemble` grows contigs across reads. Reads that tile one long CDR3 in different shards
  never meet, so the contig is never built -- which is precisely the class Stage 3 exists for.

So the sharded path runs `map` per shard, concatenates the Stage-1 AIRR **in shard order**,
and then calls exactly the same :func:`finish` the single-node path calls. Not "the same
steps" -- the same function, so the two cannot drift apart in a parameter.

That plus contiguous shards (:func:`arda.cluster.split_pairs`) is what makes a sharded run
**byte-identical** to a single-node one, rather than merely similar.
"""

from __future__ import annotations

import json
from pathlib import Path

from .. import __version__
from .._log import logger
from ._res import Stage

__all__ = ["finish", "run", "reduce", "write_stats_for", "OUTPUTS"]

#: Output basenames, relative to ``out_dir`` and given ``prefix``.
OUTPUTS = {
    "airr": "{prefix}.airr.tsv",
    "assembled_airr": "{prefix}.assembled.airr.tsv",
    "clones": "{prefix}.clones.tsv",
    "report": "{prefix}.arda.json",
    "stats": "{prefix}.stats.tsv",
}


def _provenance() -> dict:
    """What was actually used — so a cross-mode divergence is diagnosable, not mysterious.

    Two modes disagreeing is otherwise a long hunt. The usual causes are a different aligner
    build (conda's mmseqs vs the static one) or a different reference (a conda env that fetched
    its own), and both are invisible in the outputs. Cheap to record, so record it.
    """
    info: dict = {"arda_version": __version__}
    try:
        from .. import mmseqs

        info["mmseqs_version"] = mmseqs.version()
    except Exception:  # noqa: BLE001 — provenance must never break a run
        info["mmseqs_version"] = None
    try:
        from ..annotate.mapper import Reference  # noqa: F401 — import guard only
        from ..paths import vdj_dir

        alleles = vdj_dir() / "human" / "alleles.fasta"
        if alleles.exists():
            st = alleles.stat()
            info["reference"] = {"path": str(alleles), "bytes": st.st_size,
                                 "mtime": int(st.st_mtime)}
        else:
            info["reference"] = None
    except Exception:  # noqa: BLE001
        info["reference"] = None
    return info


def finish(airr: str | Path, out_dir: str | Path, out_prefix: str, *,
           organism: str = "human", threads: int = 0, assemble: bool = True,
           complete_only: bool = True, map_d: bool = True,
           d_max_evalue: float | None = None,
           ec_mode: str = "fast", min_junction_q: int | None = None,
           clonotype_key: str = "full", call_level: str = "allele", isotype: bool = True,
           map_report: dict | None = None, write_qc: bool = True, echo=None) -> dict:
    """Run Stages 2-3 over a Stage-1 AIRR and write the clonotype table + merged report.

    Called by both :func:`run` (single node) and :func:`reduce` (after a sharded Stage 1).

    Args:
        airr: Stage-1 mapped-reads AIRR TSV.
        out_dir, out_prefix: where the outputs land (see :data:`OUTPUTS`).
        map_report: the Stage-1 report to embed; for a sharded run, the merged one.
        echo: optional ``print``-like callback for progress lines.

    Returns:
        The merged report dict, as written to ``<prefix>.arda.json``.
    """
    from .correct import correct_airr

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    airr = Path(airr)
    paths = {k: out_dir / v.format(prefix=out_prefix) for k, v in OUTPUTS.items()}
    say = echo or logger.info

    stage = Stage()
    arep = None
    extra = None
    if assemble:
        from .assemble import assemble_contigs

        extra = paths["assembled_airr"]
        arep = assemble_contigs(airr, extra, organism=organism, threads=threads, map_d=map_d,
                                d_max_evalue=d_max_evalue)
        say(f"assemble: {arep.contigs_complete}/{arep.contigs} complete contigs from "
            f"{arep.seeds} seeds; rescued {arep.reads_rescued} reads")

    crep = correct_airr(airr, paths["clones"], organism=organism, map_d=map_d,
                        d_max_evalue=d_max_evalue, ec_mode=ec_mode,
                        min_junction_q=min_junction_q, clonotype_key=clonotype_key,
                        call_level=call_level, isotype=isotype,
                        complete_only=complete_only, extra_airr=extra)
    say(f"correct: {crep.clonotypes_in} -> {crep.clonotypes_out} clonotypes "
        f"({crep.collapsed} collapsed) over {crep.reads} reads")

    report = {
        **_provenance(),
        "wall_seconds": round(stage.wall_seconds, 3),
        "map": map_report,
        "assemble": arep.as_dict() if arep else None,
        "correct": crep.as_dict(),
    }
    paths["report"].write_text(json.dumps(report, indent=2) + "\n")
    if write_qc:
        write_stats_for(out_dir, out_prefix, organism=organism, say=say)
    return report


def write_stats_for(out_dir: str | Path, out_prefix: str, *, organism: str = "human",
                    say=None) -> int:
    """Write ``<prefix>.stats.tsv`` from the run's own artifacts. Returns the row count.

    ⛔ Written unconditionally, not behind a flag. It reads only files that already exist and
    costs one pass over each; the alternative is that the numbers an operator needs to decide
    whether a sample is usable exist only if they knew to ask for them BEFORE the run.

    ⛔ Called AFTER the report JSON is final. :func:`run` rewrites it with the whole-run wall time
    once Stages 2-3 return, so collecting inside :func:`finish` would put the Stage-2/3 time in
    the ``run`` scope under the name ``wall_seconds`` -- a wrong number that looks like a right one.
    """
    from ..stats import collect, write_stats

    out_dir = Path(out_dir)
    paths = {k: out_dir / v.format(prefix=out_prefix) for k, v in OUTPUTS.items()}
    rows = collect(airr=paths["airr"] if paths["airr"].exists() else None,
                   clones=paths["clones"] if paths["clones"].exists() else None,
                   report=paths["report"], organism=organism)
    write_stats(rows, paths["stats"])
    (say or logger.info)(f"stats: {len(rows)} rows -> {paths['stats']}")
    return len(rows)


def run(r1: str | Path, out_dir: str | Path, out_prefix: str, *,
        r2: str | Path | None = None, organism: str = "human", threads: int = 0,
        reconstruct: bool = False, min_score: float = 75.0, kmer: int | None = 12,
        assemble: bool = True, complete_only: bool = True, map_d: bool = True,
        d_max_evalue: float | None = None,
        limit: int | None = None, two_pass: bool = False, adaptive: bool = False,
        fast_segments: bool = False, prefilter: bool = False,
        segment_only_v: bool = False, indel_rescue: bool = False,
        ec_mode: str = "fast", min_junction_q: int | None = None,
        clonotype_key: str = "full", call_level: str = "allele", isotype: bool = True,
        shm: str = "framework",
        complete_junction_nt: int = 0,
        echo=None) -> dict:
    """Single-node map -> assemble -> correct."""
    from .map import map_rnaseq

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    airr = out_dir / OUTPUTS["airr"].format(prefix=out_prefix)
    say = echo or logger.info

    whole = Stage()
    mrep = map_rnaseq(r1, airr, r2=r2, organism=organism, threads=threads,
                      reconstruct=reconstruct, min_score=min_score, map_d=map_d,
                      d_max_evalue=d_max_evalue,
                      kmer=kmer, limit=limit, two_pass=two_pass, adaptive=adaptive,
                      fast_segments=fast_segments, prefilter=prefilter,
                      segment_only_v=segment_only_v, indel_rescue=indel_rescue,
                      # ⛔ The quality gate reads a column Stage 1 only writes when asked. In `run`
                      # both stages happen in one call, so the user cannot wire that up by hand --
                      # asking for the gate has to imply producing its input, or `--ec-mode
                      # accurate` would silently do nothing here.
                      with_junction_quality=(ec_mode != "fast" or min_junction_q is not None),
                      shm=shm, complete_junction_nt=complete_junction_nt)

    report = finish(airr, out_dir, out_prefix, organism=organism, threads=threads,
                    assemble=assemble, complete_only=complete_only, map_d=map_d,
                    ec_mode=ec_mode, min_junction_q=min_junction_q,
                    clonotype_key=clonotype_key, call_level=call_level, isotype=isotype,
                    d_max_evalue=d_max_evalue,
                    map_report=mrep.as_dict(), write_qc=False, echo=echo)
    report["wall_seconds"] = round(whole.wall_seconds, 3)
    (out_dir / OUTPUTS["report"].format(prefix=out_prefix)).write_text(
        json.dumps(report, indent=2) + "\n")
    write_stats_for(out_dir, out_prefix, organism=organism, say=say)
    return report


def _merge_map_reports(paths: list[Path]) -> dict:
    """Combine per-shard Stage-1 reports into one.

    Counts sum. Wall time and RSS deliberately do **not** collapse to a single
    ``wall_seconds`` / ``peak_rss_mb``: adding up 40 array tasks' wall time and calling it
    "wall seconds" would be a lie, and the max is what actually sized the job. Both are
    reported under explicit names so neither can be mistaken for a single-node number.
    """
    shards = [json.loads(Path(p).read_text()) for p in paths]
    if not shards:
        return {}
    per_locus: dict[str, int] = {}
    for s in shards:
        for locus, n in (s.get("per_locus") or {}).items():
            per_locus[locus] = per_locus.get(locus, 0) + int(n)
    total = sum(int(s.get("total_reads", 0)) for s in shards)
    mapped = sum(int(s.get("mapped_reads", 0)) for s in shards)
    return {
        "shards": len(shards),
        "organism": shards[0].get("organism"),
        "total_reads": total,
        "mapped_reads": mapped,
        "mapped_fraction": (mapped / total) if total else 0.0,
        "per_locus": dict(sorted(per_locus.items())),
        "constant_only_fragments": sum(int(s.get("constant_only_fragments", 0)) for s in shards),
        "isotype_from_mate": sum(int(s.get("isotype_from_mate", 0)) for s in shards),
        "min_score": shards[0].get("min_score"),
        "threads": shards[0].get("threads"),
        "wall_seconds_max": max(float(s.get("wall_seconds", 0.0)) for s in shards),
        "wall_seconds_sum": round(sum(float(s.get("wall_seconds", 0.0)) for s in shards), 3),
        "peak_rss_mb_max": max(float(s.get("peak_rss_mb", 0.0)) for s in shards),
    }


def reduce(shard_dir: str | Path, out_dir: str | Path, out_prefix: str, *,
           organism: str = "human", threads: int = 0, assemble: bool = True,
           complete_only: bool = True, map_d: bool = True,
           d_max_evalue: float | None = None, echo=None) -> dict:
    """Merge sharded Stage-1 output, then run Stages 2-3 once over the whole thing.

    The shard AIRRs are merged from an **explicit sorted list**, not a bare ``*.tsv`` glob:
    shard names are zero-padded so ``sorted()`` is numeric (``shard_10`` must not precede
    ``shard_2`` -- concatenation order *is* read order here), and naming the glob means
    `reduce` can never swallow its own ``clones.tsv`` if someone points ``--out-dir`` at
    ``--shard-dir``.
    """
    from ..cluster import RNASEQ_SHARD_GLOB, merge

    shard_dir, out_dir = Path(shard_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    say = echo or logger.info

    shards = sorted(shard_dir.glob(RNASEQ_SHARD_GLOB))
    if not shards:
        raise FileNotFoundError(f"no {RNASEQ_SHARD_GLOB} under {shard_dir}")
    airr = out_dir / OUTPUTS["airr"].format(prefix=out_prefix)
    merge(shards, airr)
    say(f"merged {len(shards)} shard AIRRs -> {airr}")

    mrep = _merge_map_reports(sorted(shard_dir.glob("shard_*.map.json")))
    if mrep:
        say(f"map (summed over {mrep['shards']} shards): "
            f"{mrep['mapped_reads']}/{mrep['total_reads']} reads mapped; loci={mrep['per_locus']}")
    return finish(airr, out_dir, out_prefix, organism=organism, threads=threads,
                  assemble=assemble, complete_only=complete_only, map_d=map_d,
                  d_max_evalue=d_max_evalue, map_report=mrep or None, echo=echo)
