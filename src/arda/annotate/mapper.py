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
from ..paths import data_dir
from ..refbuild.translate import reverse_complement
from . import io as seqio
from .reference import load_reference
from .transfer import transfer_hit, AIRR_COLUMNS
from .airr_out import write_airr

__all__ = ["annotate_file", "annotate_records"]

_SEARCH_TYPE = {"nt": mmseqs.SEARCH_TYPE_NUCLEOTIDE, "aa": mmseqs.SEARCH_TYPE_PROTEIN}

# Tuned defaults (see memory/mmseqs-params.md). Short germline-similar queries:
# moderately high sensitivity, keep only a few best hits, no coverage filter so
# partial RNA-seq reads still map. The whole reference (all loci) is one DB, so a
# single search annotates mixed bulk RNA-seq across all loci at once.
_SENSITIVITY = {"nt": 7.0, "aa": 7.0}
_MAX_SEQS = 50


def _cached_target_db(target_fasta: Path, organism: str, seqtype: str) -> Path:
    """Build (once) and reuse an mmseqs DB for the reference scaffolds.

    Cached under ``data/mmseqs_db/<organism>_<seqtype>``; rebuilt if the source
    FASTA is newer than the cached DB. Avoids re-creating the target DB on every
    annotation call (the dominant overhead for small inputs).
    """
    cache = data_dir() / "mmseqs_db" / f"{organism}_{seqtype}"
    cache.mkdir(parents=True, exist_ok=True)
    db = cache / "db"
    dbtype = 2 if seqtype == "nt" else 1
    if not db.exists() or db.stat().st_mtime < Path(target_fasta).stat().st_mtime:
        mmseqs.createdb(target_fasta, db, dbtype=dbtype)
    return db


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
    sensitivity: float | None = None,
    strand: str = "both",
) -> list[dict]:
    """Annotate in-memory ``(id, sequence)`` records; return AIRR record dicts.

    Args:
        strand: ``"both"`` (default, nt only) searches both strands and re-orients
            reverse-complement hits; ``"forward"`` searches the plus strand only.
            Ignored for protein input.
    """
    if seqtype not in _SEARCH_TYPE:
        raise ValueError(f"seqtype must be 'nt' or 'aa', got {seqtype!r}")
    ref = load_reference(organism, seqtype)
    threads = threads or (os.cpu_count() or 1)
    sensitivity = _SENSITIVITY[seqtype] if sensitivity is None else sensitivity
    mm_strand = (2 if strand == "both" else 1) if seqtype == "nt" else None

    target_db = _cached_target_db(ref.target_fasta, organism, seqtype)

    with tempfile.TemporaryDirectory(prefix="arda_") as td:
        tmp = Path(td)
        query_fa = seqio.write_fasta(iter(records), tmp / "query.fasta")
        query_db = tmp / "queryDB"
        res_db = tmp / "resDB"
        out_tsv = tmp / "hits.tsv"
        mmseqs.createdb(query_fa, query_db, dbtype=2 if seqtype == "nt" else 1)
        mmseqs.search(
            query_db, target_db, res_db, tmp / "mmseqs_tmp",
            search_type=_SEARCH_TYPE[seqtype], sensitivity=sensitivity,
            max_seqs=_MAX_SEQS, threads=threads,
            extra=(["--strand", str(mm_strand)] if mm_strand is not None else None),
        )
        mmseqs.convertalis(query_db, target_db, res_db, out_tsv, threads=threads,
                           search_type=_SEARCH_TYPE[seqtype])
        best = _best_hits(out_tsv)

    out: list[dict] = []
    for qid, qseq in records:
        hit = best.get(qid)
        entry = ref.get(hit["target"]) if hit else None
        if hit is None or entry is None:
            rec = {c: "" for c in AIRR_COLUMNS}
            rec["sequence_id"], rec["sequence"] = qid, qseq
            out.append(rec)
            continue
        # mmseqs reports reverse-strand nt hits with qstart > qend and aligned
        # strings already on the coding strand. Re-orient: work on the revcomp,
        # remap the alignment start to forward coords on it.
        qs, qe = int(hit["qstart"]), int(hit["qend"])
        rev = qs > qe
        work = qseq
        if rev:
            work = reverse_complement(qseq)
            qlen = len(qseq)
            hit = dict(hit)
            hit["qstart"], hit["qend"] = qlen - qs + 1, qlen - qe + 1
        out.append(transfer_hit(qid, work, hit, entry, seqtype, rev_comp=rev))
    return out


def annotate_file(
    input: str | Path,
    output: str | Path,
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    strand: str = "both",
) -> Path:
    """Annotate a FASTA/FASTQ file and write an AIRR TSV."""
    records = list(seqio.read_sequences(input))
    airr = annotate_records(records, organism=organism, seqtype=seqtype,
                            threads=threads, strand=strand)
    return write_airr(airr, output)
