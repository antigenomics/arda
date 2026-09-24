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

__all__ = ["Sample", "load", "read_sheet", "SHEET_COLUMNS", "LABEL_COLUMNS"]

#: The sample-sheet columns arda reads. Deliberately nf-core's spelling, so an existing
#: nf-core samplesheet works here unmodified.
SHEET_COLUMNS = ("sample", "fastq_1", "fastq_2", "project", "batch")

#: The two sheet DIALECTS arda reads, keyed by the column that identifies each.
#:
#: ``nf-core`` is the generic spelling above. ``airrflow`` is `nf-core/airrflow
#: <https://github.com/nf-core/airrflow>`_'s AIRR-metadata samplesheet, which is the community
#: standard for AIRR-seq and carries the protocol and subject fields a repertoire actually needs.
#: Reading both HERE is what lets one samplesheet drive the Nextflow module, the Snakemake
#: workflow and ``arda cluster`` alike, instead of each integration inventing its own.
#:
#: Never: a sheet is one dialect or the other, never both. A header carrying ``sample`` *and*
#: ``sample_id`` does not say which one names the repertoire, and guessing splits or merges it
#: silently.
_DIALECTS: dict[str, dict[str, str]] = {
    "nf-core": {"id": "sample", "r1": "fastq_1", "r2": "fastq_2"},
    "airrflow": {"id": "sample_id", "r1": "filename_R1", "r2": "filename_R2"},
}

#: airrflow columns arda reads. ``species`` picks the reference organism and ``single_cell``
#: picks ``arda cells``, so those two CHANGE THE COMMAND; the rest are labels. ``subject_id``
#: becomes :attr:`Sample.project` so ``arda qc batch`` groups a cohort by donor without being
#: told twice.
_AIRRFLOW_META = ("species", "pcr_target_locus", "single_cell", "subject_id", "tissue", "sex",
                  "age", "biomaterial_provider", "filename_I1", "filename")

#: The two that are LABELS, not inputs. Nothing in the pipeline reads them and no output changes
#: because of them; they exist so ``arda qc batch`` can group a cohort the way it was collected,
#: and so a sheet that already carries them stops being warned about.
LABEL_COLUMNS = ("project", "batch")

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
    #: Free-text grouping labels from the sheet. Empty when it did not carry them.
    project: str = ""
    batch: str = ""
    #: From an airrflow-dialect sheet; empty/False when the sheet did not carry them. ``species``
    #: and ``single_cell`` are the two that change what arda is asked to run.
    species: str = ""
    locus: str = ""
    single_cell: bool = False


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


def _group(rows: list[tuple[str, Path, Path | None]], where: str,
           labels: dict[str, dict[str, str | bool]] | None = None) -> list[Sample]:
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
        lab = (labels or {}).get(sid, {})
        out.append(Sample(sid, tuple(pairs),
                          str(lab.get("project", "")), str(lab.get("batch", "")),
                          str(lab.get("species", "")), str(lab.get("locus", "")),
                          bool(lab.get("single_cell", False))))
    return out


def _exists(p: Path, where: str) -> Path:
    if not p.exists():
        raise ValueError(f"{where}: no such file: {p}")
    return p


def read_sheet(path: str | Path) -> list[Sample]:
    """Read a sample sheet in either supported dialect.

    * **nf-core** -- ``sample``, ``fastq_1``, optional ``fastq_2``, ``project``, ``batch``.
    * **nf-core/airrflow** -- ``sample_id``, ``filename_R1``, optional ``filename_R2``, plus the
      AIRR metadata its schema requires (``species``, ``pcr_target_locus``, ``single_cell``,
      ``subject_id``, ``tissue``, ...). ``species`` and ``single_cell`` are the two that change
      what arda is asked to run; ``subject_id`` becomes :attr:`Sample.project` and ``tissue``
      :attr:`Sample.batch`, so ``arda qc batch`` groups a cohort by donor and tissue for free.

    The dialect is chosen by which id column the header carries, and a header carrying **both** is
    refused: it does not say which column names the repertoire. Reading both dialects here is what
    lets ONE samplesheet drive the Nextflow module, the Snakemake workflow and ``arda cluster``.

    TSV unless the name ends ``.csv``. Repeated ids merge into one sample in row order --
    nf-core's re-sequencing rule, and the reason an nf-core or airrflow samplesheet works here as
    written. Extra columns are ignored, but **named once in a warning**: a sheet whose header says
    ``fastq2`` is not a sheet with no R2, and silently treating it as single-end halves the data.

    Labels are read from the sample's FIRST row -- a sample's read groups belong to one subject,
    species and batch by definition, and a sheet that says otherwise has a typo rather than a
    meaning.

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
        dialect, cols = _pick_dialect(header, path, delim)
        known = set(cols.values()) | set(LABEL_COLUMNS)
        if dialect == "airrflow":
            known |= set(_AIRRFLOW_META)
        extra = [h for h in header if h and h not in known]
        if extra:
            logger.warning("%s: ignoring unknown sample-sheet column(s): %s",
                           path.name, ", ".join(extra))

        rows: list[tuple[str, Path, Path | None]] = []
        labels: dict[str, dict[str, str | bool]] = {}
        for n, row in enumerate(reader, start=2):        # 1 is the header
            where = f"{path}:{n}"
            sid = _check_id(row.get(cols["id"]) or "", where)
            r1 = (row.get(cols["r1"]) or "").strip()
            if not r1:
                raise ValueError(f"{where}: empty {cols['r1']}")
            r2 = (row.get(cols["r2"]) or "").strip()
            rows.append((sid,
                         _exists(_resolve(r1, base), where),
                         _exists(_resolve(r2, base), where) if r2 else None))
            labels.setdefault(sid, _row_labels(row, dialect, where))
    if not rows:
        raise ValueError(f"{path}: sample sheet has a header but no rows")
    return _group(rows, str(path), labels)


def _pick_dialect(header: list[str], path: Path, delim: str) -> tuple[str, dict[str, str]]:
    """Which dialect this header is, or a message naming both.

    Never: a header carrying BOTH id columns is refused rather than resolved. It does not say
    which one names the repertoire, and picking one silently splits or merges a sample.
    """
    seen = [name for name, cols in _DIALECTS.items() if cols["id"] in header]
    if len(seen) > 1:
        raise ValueError(
            f"{path}: sample sheet carries {' and '.join(_DIALECTS[n]['id'] for n in seen)}, so it "
            f"does not say which column names the repertoire. Keep one dialect per sheet.")
    if not seen:
        raise ValueError(
            f"{path}: sample sheet has no recognised id column. Expected either "
            f"{_DIALECTS['nf-core']['id']} + {_DIALECTS['nf-core']['r1']} (nf-core) or "
            f"{_DIALECTS['airrflow']['id']} + {_DIALECTS['airrflow']['r1']} (nf-core/airrflow), "
            f"{'comma' if delim == ',' else 'tab'}-separated; "
            f"saw {', '.join(header) or '(no header)'}")
    cols = _DIALECTS[seen[0]]
    if cols["r1"] not in header:
        raise ValueError(
            f"{path}: sample sheet is missing {cols['r1']}. Its header carries {cols['id']}, so "
            f"it is a {seen[0]}-dialect sheet, and that dialect needs both; "
            f"saw {', '.join(header) or '(no header)'}")
    return seen[0], cols


#: airrflow writes booleans as TRUE/FALSE; a sheet written by hand may say yes/1/true.
_TRUE = {"true", "t", "yes", "y", "1"}


def _row_labels(row: dict, dialect: str, where: str) -> dict[str, str | bool]:
    """The non-input columns of one row, normalised.

    Read from the sample's FIRST row: a sample's read groups belong to one subject, species and
    protocol by definition, and a sheet that says otherwise has a typo rather than a meaning.
    """
    if dialect == "nf-core":
        return {"project": (row.get("project") or "").strip(),
                "batch": (row.get("batch") or "").strip()}
    sc = (row.get("single_cell") or "").strip()
    species = (row.get("species") or "").strip().lower()
    if species and species not in ("human", "mouse"):
        # airrflow's own schema is an enum of exactly these two; arda ships five organisms, so a
        # third is passed through rather than refused -- but say so, because a typo here silently
        # selects the wrong reference.
        logger.warning("%s: species %r is outside airrflow's human|mouse enum; passing it "
                       "through to --organism unchanged", where, species)
    return {"project": (row.get("subject_id") or "").strip(),   # qc batch groups by donor
            "batch": (row.get("tissue") or "").strip(),
            "species": species,
            "locus": (row.get("pcr_target_locus") or "").strip().upper(),
            "single_cell": sc.lower() in _TRUE}


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
