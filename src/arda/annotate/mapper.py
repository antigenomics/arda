"""Runtime annotation: map input sequences to the reference and transfer markup.

Pipeline: read input (FASTA/FASTQ) -> MMseqs2 search against the curated scaffold
DB -> best hit per query -> project reference region markup onto the query (C++
hot path) -> AIRR TSV.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import polars as pl

from .. import mmseqs
from . import io as seqio
from .reference import load_reference
from .transfer import transfer_hit, AIRR_COLUMNS
from .airr_out import write_airr

__all__ = ["annotate_file", "annotate_records"]

_SEARCH_TYPE = {"nt": mmseqs.SEARCH_TYPE_NUCLEOTIDE, "aa": mmseqs.SEARCH_TYPE_PROTEIN}


def _best_hits(tsv: Path) -> dict[str, dict]:
    """Parse the mmseqs TSV and return the top-scoring hit per query."""
    cols = mmseqs.DEFAULT_FORMAT_OUTPUT.split(",")
    if tsv.stat().st_size == 0:  # no hits at all
        return {}
    df = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=cols,
                     infer_schema_length=0)
    if df.height == 0:
        return {}
    df = df.with_columns(pl.col("bits").cast(pl.Float64, strict=False))
    df = df.sort("bits", descending=True).unique(subset="query", keep="first")
    return {row["query"]: row for row in df.iter_rows(named=True)}


def annotate_records(
    records: list[tuple[str, str]],
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    sensitivity: float = 5.7,
) -> list[dict]:
    """Annotate in-memory ``(id, sequence)`` records; return AIRR record dicts."""
    if seqtype not in _SEARCH_TYPE:
        raise ValueError(f"seqtype must be 'nt' or 'aa', got {seqtype!r}")
    ref = load_reference(organism, seqtype)
    threads = threads or (os.cpu_count() or 1)
    seqs = dict(records)

    with tempfile.TemporaryDirectory(prefix="arda_") as td:
        tmp = Path(td)
        query_fa = seqio.write_fasta(iter(records), tmp / "query.fasta")
        out_tsv = tmp / "hits.tsv"
        mmseqs.easy_search(
            query_fa, ref.target_fasta, out_tsv, tmp / "mmseqs_tmp",
            search_type=_SEARCH_TYPE[seqtype], sensitivity=sensitivity, threads=threads,
        )
        best = _best_hits(out_tsv)

    out: list[dict] = []
    for qid, qseq in records:
        hit = best.get(qid)
        entry = ref.get(hit["target"]) if hit else None
        if hit is None or entry is None:
            rec = {c: "" for c in AIRR_COLUMNS}
            rec["sequence_id"], rec["sequence"] = qid, qseq
            out.append(rec)
        else:
            out.append(transfer_hit(qid, qseq, hit, entry, seqtype))
    return out


def annotate_file(
    input: str | Path,
    output: str | Path,
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
) -> Path:
    """Annotate a FASTA/FASTQ file and write an AIRR TSV."""
    records = list(seqio.read_sequences(input))
    airr = annotate_records(records, organism=organism, seqtype=seqtype, threads=threads)
    return write_airr(airr, output)
