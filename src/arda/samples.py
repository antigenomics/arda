"""What a *sample* is, when it is not one FASTQ pair.

A sample arrives split across files more often than not: Illumina writes one FASTQ per lane
(``SampleName_S1_L001_R1_001.fastq.gz``, ``..._L002_...``), and in-house pipelines chunk a run
into numbered parts. All of them are one repertoire and must produce **one** clonotype table.

The grouping is **declared, never inferred**. There is no filename pattern matching here and there
should not be: given ``RNA-SAMPLE_ID:12:00XX919:3_1.fastq.gz`` no rule can say which colon-separated
field is the sample without being told, and a rule that guesses wrong silently splits one repertoire
into four. So the id comes from ``--id`` or from a sample sheet, and nothing else.

**A multi-file sample is an already-sharded sample.** :func:`arda.cluster.split_pairs` exists to cut
one pair into contiguous blocks precisely so the per-block Stage-1 AIRR can be concatenated in block
order and fed to Stages 2-3 once; a sample delivered as several files has simply been sharded by
whoever wrote it. So the files of one sample are its shards, in declared order, and nothing in
`correct` or `assemble` needs to know samples exist.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from ._log import logger

__all__ = ["Sample", "load", "read_sheet", "SHEET_COLUMNS"]

#: The sample-sheet columns arda reads. Deliberately nf-core's spelling, so an existing
#: nf-core samplesheet works here unmodified.
SHEET_COLUMNS = ("sample", "fastq_1", "fastq_2")

#: Characters a sample id may not contain: it becomes an output *filename*.
_ID_FORBIDDEN = set('/\\:*?"<>|')


@dataclass(frozen=True)
class Sample:
    """One repertoire and the read-group pairs it was delivered in, **in declared order**.

    Declared order is shard order, and shard order is row order in the merged Stage-1 AIRR, which
    is what makes a multi-file run byte-identical to the same reads concatenated. Do not sort.
    """

    id: str
    pairs: tuple[tuple[Path, Path | None], ...]


def _check_id(sid: str, where: str) -> str:
    sid = sid.strip()
    if not sid:
        raise ValueError(f"{where}: empty sample id")
    bad = sorted(_ID_FORBIDDEN & set(sid))
    if bad:
        raise ValueError(
            f"{where}: sample id {sid!r} contains {''.join(bad)!r}, and the id becomes an output "
            f"filename")
    return sid


def _group(rows: list[tuple[str, Path, Path | None]], where: str) -> list[Sample]:
    """Merge rows sharing an id, keeping first-appearance order for ids and row order within one.

    Never: a sample is single-end or paired, never both. Half its reads silently losing their mate
    is not a configuration, and it reaches `map` as a length mismatch several minutes later.
    """
    order: list[str] = []
    by_id: dict[str, list[tuple[Path, Path | None]]] = {}
    for sid, r1, r2 in rows:
        if sid not in by_id:
            order.append(sid)
            by_id[sid] = []
        by_id[sid].append((r1, r2))
    out = []
    for sid in order:
        pairs = by_id[sid]
        if len({p[1] is None for p in pairs}) != 1:
            raise ValueError(
                f"{where}: sample {sid!r} mixes paired and single-end read groups; give every "
                f"read group an R2, or none of them")
        out.append(Sample(sid, tuple(pairs)))
    return out


def _exists(p: Path, where: str) -> Path:
    if not p.exists():
        raise ValueError(f"{where}: no such file: {p}")
    return p


def read_sheet(path: str | Path) -> list[Sample]:
    """Read a sample sheet. Columns ``sample``, ``fastq_1`` and optionally ``fastq_2``.

    TSV unless the name ends ``.csv``. Repeated ``sample`` values merge into one sample in row
    order -- nf-core's re-sequencing rule, and the reason an nf-core samplesheet works here as
    written. Extra columns are ignored, but **named once in a warning**: a sheet whose header says
    ``fastq2`` is not a sheet with no R2, and silently treating it as single-end halves the data.

    Relative paths resolve against the **sheet's own directory**, so a sheet travels with its data.
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"--samples: no such file: {path}")
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    base = path.parent

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        header = [h.strip() for h in (reader.fieldnames or [])]
        missing = [c for c in SHEET_COLUMNS[:2] if c not in header]
        if missing:
            raise ValueError(
                f"{path}: sample sheet is missing {', '.join(missing)}. Expected the columns "
                f"{', '.join(SHEET_COLUMNS)} (fastq_2 optional), {'comma' if delim == ',' else 'tab'}"
                f"-separated; saw {', '.join(header) or '(no header)'}")
        extra = [h for h in header if h and h not in SHEET_COLUMNS]
        if extra:
            logger.warning("%s: ignoring unknown sample-sheet column(s): %s",
                           path.name, ", ".join(extra))

        rows: list[tuple[str, Path, Path | None]] = []
        for n, row in enumerate(reader, start=2):        # 1 is the header
            where = f"{path}:{n}"
            sid = _check_id(row.get("sample") or "", where)
            r1 = (row.get("fastq_1") or "").strip()
            if not r1:
                raise ValueError(f"{where}: empty fastq_1")
            r2 = (row.get("fastq_2") or "").strip()
            rows.append((sid,
                         _exists(_resolve(r1, base), where),
                         _exists(_resolve(r2, base), where) if r2 else None))
    if not rows:
        raise ValueError(f"{path}: sample sheet has a header but no rows")
    return _group(rows, str(path))


def _resolve(value: str, base: Path) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else (base / p)


def load(*, r1: list[Path] | None = None, r2: list[Path] | None = None,
         ids: list[str] | None = None, sheet: str | Path | None = None,
         out_prefix: str | None = None) -> list[Sample]:
    """Resolve the CLI's input options into samples, or explain exactly why they do not.

    One pair and no ``--id`` is the historical single-sample invocation and keeps using
    ``--out-prefix``. Anything else names its samples: ``--out-prefix`` cannot name two of them,
    so asking for both is rejected rather than silently overwriting one output with the other.
    """
    r1 = list(r1 or [])
    r2 = list(r2 or [])
    ids = list(ids or [])

    if sheet is not None:
        if r1 or r2 or ids:
            raise ValueError("--samples is exclusive with --r1 / --r2 / --id: "
                             "the sheet already says which files belong to which sample")
        if out_prefix:
            raise ValueError("--samples names every sample; --out-prefix cannot name them all")
        return read_sheet(sheet)

    if not r1:
        raise ValueError("give --r1 (repeatable) or --samples")
    if r2 and len(r2) != len(r1):
        raise ValueError(f"--r2 given {len(r2)} time(s) against {len(r1)} --r1: they are matched "
                         f"by position, so give one --r2 per --r1 or none at all")
    if ids and len(ids) != len(r1):
        raise ValueError(f"--id given {len(ids)} time(s) against {len(r1)} --r1: they are matched "
                         f"by position, so give one --id per --r1")
    if len(r1) > 1 and not ids:
        raise ValueError(
            f"{len(r1)} --r1 given with no --id. arda will not guess which files share a sample "
            f"from their names; give one --id per --r1 (repeat an id to merge those read groups "
            f"into one sample), or use --samples")
    if ids and out_prefix:
        raise ValueError("--id names every sample; --out-prefix is redundant and ambiguous")
    if not ids and not out_prefix:
        raise ValueError("give --out-prefix (one sample) or --id per --r1")

    if not ids:
        ids = [_check_id(out_prefix or "", "--out-prefix")]
    else:
        ids = [_check_id(s, f"--id #{i + 1}") for i, s in enumerate(ids)]

    rows = [(ids[i],
             _exists(Path(a).expanduser(), f"--r1 #{i + 1}"),
             _exists(Path(r2[i]).expanduser(), f"--r2 #{i + 1}") if r2 else None)
            for i, a in enumerate(r1)]
    return _group(rows, "--id")
