"""Runtime annotation: map input sequences to the reference and transfer markup.

Pipeline: read input (FASTA/FASTQ) -> MMseqs2 search against the curated scaffold
DB -> best hit per query -> project reference region markup onto the query (C++
hot path) -> AIRR TSV.
"""

from __future__ import annotations

import os
import queue
import tempfile
import threading
from pathlib import Path

import polars as pl

from .. import mmseqs
from ..paths import data_dir, vdj_dir
from ..refbuild.translate import reverse_complement
from . import io as seqio
from .reference import load_reference, Reference
from .transfer import transfer_hit, AIRR_COLUMNS
from .airr_out import airr_header, format_rows

__all__ = ["annotate_file", "annotate_records", "build_index"]

# Streaming defaults: process the input in bounded chunks so memory stays flat
# regardless of input size (a 30M-read FASTQ never fully loads).
_CHUNK_SIZE = 50_000

_SEARCH_TYPE = {"nt": mmseqs.SEARCH_TYPE_NUCLEOTIDE, "aa": mmseqs.SEARCH_TYPE_PROTEIN}

# Tuned defaults (see memory/mmseqs-params.md). Short germline-similar queries:
# moderately high sensitivity, keep only a few best hits, no coverage filter so
# partial RNA-seq reads still map. The whole reference (all loci) is one DB, so a
# single search annotates mixed bulk RNA-seq across all loci at once.
_SENSITIVITY = {"nt": 7.0, "aa": 7.0}
# Max MMseqs2 target hits kept per read, i.e. how many V-J scaffolds may compete before
# best-hit selection. This is NOT a cosmetic knob: the V(D)J loci are paralog-dense, and at 50
# the true scaffold was being truncated out of the candidate list *before* `_best_hits` ever saw
# it. Measured on PRJNA371303 TRA amplicon (20k reads, arda-benchmark RESULTS.md §16):
#
#   max_seqs   V-gene concord.   J-gene   junction_aa exact   wall(amplicon)  wall(RNA-seq)
#       50         83.4 %        76.0 %       84.77 %            7.6 s           48 s
#      300         96.1 %        88.4 %       87.37 %           13.4 s           58 s
#
# The mapped read SET is identical either way (bit scores can only rise with more candidates —
# verified: 68 reads scored higher at 300, 0 scored lower), so the filter and the `--min-score 75`
# calibration are untouched; only the V/J/junction CALLS improve. This is why the V+J re-ranking
# experiment (RESULTS §5) was a wash: no re-ranking can recover a scaffold that was never a
# candidate. Cost is +21 % wall and 2.3x peak RSS (820 -> 1908 MB at chunk 200k); lower it via
# `--max-seqs` if memory-bound.
_MAX_SEQS = 300

# MMseqs2's nucleotide prefilter allocates a k-mer index table of 4**k entries, so peak RSS tracks
# 4**k * 8 B almost exactly and is INDEPENDENT of database size, --max-seqs, --chunk-size and
# --threads. Measured on the V+J | J+C reference, 100 k reads, 8 threads:
#
#     k       11      12      13      14      15
#     RSS   202 MB  298 MB  697 MB  2313 MB  ~8.4 GB
#     wall  3.76 s  3.60 s  4.63 s  4.15 s      --
#
# and recall/precision are INVARIANT over k=11..14 (recall 1.0000, precision .9463-.9469, FP 269-270).
# So k is a pure memory/speed knob in this range: 12 is both the fastest measured and 2.3x smaller
# than 13. arda-benchmark OPTIMIZATION.md §6.3c. Nucleotide only: the aa prefilter is a different index.
_KMER = {"nt": 12, "aa": None}


# Columns of DEFAULT_FORMAT_OUTPUT that are genuinely text; everything else is numeric.
_STR_COLS = frozenset({"query", "target", "cigar", "qaln", "taln"})


def _dbtype(seqtype: str) -> int:
    return 2 if seqtype == "nt" else 1


def _committed_index(organism: str, seqtype: str, target_fasta: Path | None = None) -> Path | None:
    """A precompiled mmseqs DB shipped in ``database/`` — if it is version-matched **and current**.

    The staleness check is not optional. ``build-db`` rewrites ``alleles.fasta`` in place, and without
    this every later run silently searched the *previous* reference; the only symptom was a result that
    quietly refused to change. Found exactly that way: a rebuild adding 345 constant-region scaffolds
    produced zero constant-region hits. The ``data/mmseqs_db`` fallback below has always compared
    mtimes -- this path did not.
    """
    d = vdj_dir(organism) / "mmseqs" / seqtype
    db, ver = d / "db", d / "VERSION"
    if not (db.exists() and ver.exists()):
        return None
    if target_fasta is not None and db.stat().st_mtime < Path(target_fasta).stat().st_mtime:
        return None  # the reference was rebuilt after this index was compiled
    try:
        if ver.read_text().strip() == mmseqs.version():
            return db
    except Exception:  # noqa: BLE001 — mmseqs missing; fall through to rebuild
        return None
    return None


def _cached_target_db(target_fasta: Path, organism: str, seqtype: str) -> Path:
    """Resolve the mmseqs target DB for the reference scaffolds.

    Prefers the precompiled DB shipped in ``database/vdj/<org>/mmseqs/<seqtype>``
    (used out of the box when its mmseqs version matches). Otherwise builds once
    into a ``data/mmseqs_db`` cache and reuses it (rebuilt if the FASTA is newer).
    """
    committed = _committed_index(organism, seqtype, target_fasta)
    if committed is not None:
        return committed
    cache = data_dir() / "mmseqs_db" / f"{organism}_{seqtype}"
    cache.mkdir(parents=True, exist_ok=True)
    db = cache / "db"
    if not db.exists() or db.stat().st_mtime < Path(target_fasta).stat().st_mtime:
        mmseqs.createdb(target_fasta, db, dbtype=_dbtype(seqtype))
    return db


def build_index(organism: str = "all", *, force: bool = False) -> None:
    """(Re)build the precompiled mmseqs DBs shipped under ``database/``.

    Writes ``database/vdj/<org>/mmseqs/<seqtype>/db*`` + a ``VERSION`` marker so
    the runtime can use them out of the box (and detect a mmseqs-version mismatch).
    Skips up-to-date DBs unless ``force``.
    """
    from ..igblast import SUPPORTED_ORGANISMS
    orgs = SUPPORTED_ORGANISMS if organism == "all" else (organism,)
    ver = mmseqs.version()
    for org in orgs:
        for seqtype in ("nt", "aa"):
            fasta = vdj_dir(org) / ("alleles.aa.fasta" if seqtype == "aa" else "alleles.fasta")
            if not fasta.exists():
                continue
            out = vdj_dir(org) / "mmseqs" / seqtype
            db, vfile = out / "db", out / "VERSION"
            if db.exists() and not force and vfile.exists() and vfile.read_text().strip() == ver:
                continue
            out.mkdir(parents=True, exist_ok=True)
            for stale in out.glob("db*"):
                stale.unlink()
            mmseqs.createdb(fasta, db, dbtype=_dbtype(seqtype))
            vfile.write_text(ver + "\n")
            print(f"[arda] built mmseqs index {org}/{seqtype} ({ver})")


def _best_hits(tsv: Path) -> dict[str, dict]:
    """Parse the mmseqs TSV and return the top-scoring hit per query.

    Best = highest whole-scaffold bit score. (A V-segment-restricted re-ranking was
    tried to fix V-paralog mis-assignment on short 3'-V-anchored reads; it only traded
    V accuracy for J and cost ~17 % throughput — the disagreement is largely irreducible
    paralog ambiguity — so max-bits is kept. See arda-benchmark RESULTS.md.)
    """
    cols = mmseqs.DEFAULT_FORMAT_OUTPUT.split(",")
    if tsv.stat().st_size == 0:  # no hits at all
        return {}
    # Typed, not `infer_schema_length=0`. Reading 17 columns as Utf8 and then casting one of them was
    # a large share of peak RSS before `mmseqs.top_hit` shrank this file from 194 MB to 1 MB; keep the
    # schema explicit so it stays cheap if `top_hit` is ever bypassed.
    schema = {c: (pl.Utf8 if c in _STR_COLS else pl.Float64) for c in cols}
    df = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=cols,
                     schema_overrides=schema)
    if df.height == 0:
        return {}
    df = df.with_columns(pl.col("bits").cast(pl.Float64, strict=False))
    df = df.sort("bits", descending=True).unique(subset="query", keep="first")
    return {row["query"]: row for row in df.iter_rows(named=True)}


def _annotate_chunk(
    records: list[tuple[str, str]],
    ref: Reference,
    target_db: Path,
    seqtype: str,
    *,
    threads: int,
    sensitivity: float,
    mm_strand: int | None,
    map_d: bool = True,
    mapped_only: bool = False,
    max_seqs: int = _MAX_SEQS,
    kmer: int | None = -1,
) -> list[dict]:
    """Annotate one batch against a preloaded reference + cached target DB.

    ``mapped_only`` skips the empty record for non-hits (the RNA-seq filter path,
    where 95-99 % of reads have no hit and building throwaway records dominates).
    """
    if not records:
        return []
    if kmer == -1:  # sentinel: caller did not override, use the seqtype default
        kmer = _KMER[seqtype]
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
            max_seqs=max_seqs, threads=threads, kmer=kmer,
            extra=(["--strand", str(mm_strand)] if mm_strand is not None else None),
        )
        # Reduce to one alignment per query BEFORE materialising it as text. With --max-seqs 300 a
        # 100 k-read chunk emits ~804 k alignment rows -- 194 MB of cigar/qaln/taln -- of which we
        # keep 4 k. Parsing the other 800 k was arda's single largest memory consumer (877 MB peak,
        # vs 284 MB for mmseqs itself). Bit-identical: same target and score on all 4,101 queries.
        mmseqs.convertalis(query_db, target_db, mmseqs.top_hit(res_db, tmp / "bestDB"),
                           out_tsv, threads=threads, search_type=_SEARCH_TYPE[seqtype])
        best = _best_hits(out_tsv)

    out: list[dict] = []
    for qid, qseq in records:
        hit = best.get(qid)
        entry = ref.get(hit["target"]) if hit else None
        if hit is None or entry is None:
            if mapped_only:
                continue
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
        dg = ref.d_germlines.get(entry.locus) if (seqtype == "nt" and map_d) else None
        out.append(transfer_hit(qid, work, hit, entry, seqtype, rev_comp=rev,
                                d_germlines=dg))
    return out


def _prep(organism, seqtype, threads, sensitivity, strand):
    if seqtype not in _SEARCH_TYPE:
        raise ValueError(f"seqtype must be 'nt' or 'aa', got {seqtype!r}")
    # Reject an unknown strand rather than silently searching forward-only. On stranded paired
    # libraries R2 is antisense, so a typo here quietly discards ~40 % of the recoverable
    # repertoire (the R2-only fragments; arda-benchmark RESULTS.md §13d).
    if strand not in ("both", "forward"):
        raise ValueError(f"strand must be 'both' or 'forward', got {strand!r}")
    ref = load_reference(organism, seqtype)
    threads = threads or (os.cpu_count() or 1)
    sensitivity = _SENSITIVITY[seqtype] if sensitivity is None else sensitivity
    mm_strand = (2 if strand == "both" else 1) if seqtype == "nt" else None
    target_db = _cached_target_db(ref.target_fasta, organism, seqtype)
    return ref, target_db, threads, sensitivity, mm_strand


def annotate_records(
    records: list[tuple[str, str]],
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    sensitivity: float | None = None,
    strand: str = "both",
    map_d: bool = True,
) -> list[dict]:
    """Annotate in-memory ``(id, sequence)`` records; return AIRR record dicts.

    Args:
        strand: ``"both"`` (default, nt only) searches both strands and re-orients
            reverse-complement hits; ``"forward"`` searches the plus strand only.
            Ignored for protein input.
        map_d: ``True`` (default) maps D segments into the junction of VDJ-locus
            hits (``d_call``/``d2_call``/``np*``); ``False`` skips D mapping (nt
            input only — D mapping never runs for protein input).
    """
    ref, target_db, threads, sensitivity, mm_strand = _prep(
        organism, seqtype, threads, sensitivity, strand)
    return _annotate_chunk(records, ref, target_db, seqtype,
                           threads=threads, sensitivity=sensitivity,
                           mm_strand=mm_strand, map_d=map_d)


def annotate_file(
    input: str | Path,
    output: str | Path,
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    sensitivity: float | None = None,
    strand: str = "both",
    chunk_size: int = _CHUNK_SIZE,
    map_d: bool = True,
) -> Path:
    """Annotate a FASTA/FASTQ file and stream an AIRR TSV.

    The input is processed in bounded chunks with a background reader thread that
    prefetches the next chunk while the current one is annotated (mmseqs releases
    the GIL during its subprocess), so memory stays flat for arbitrarily large
    FASTQ and read parsing overlaps compute. The reference + target DB are loaded
    once and reused across all chunks.
    """
    output = Path(output)
    ref, target_db, threads, sensitivity, mm_strand = _prep(
        organism, seqtype, threads, sensitivity, strand)

    chunks: queue.Queue = queue.Queue(maxsize=2)

    def reader():
        try:
            for chunk in seqio.chunked(seqio.read_sequences(input), chunk_size):
                chunks.put(chunk)
        finally:
            chunks.put(None)  # sentinel

    t = threading.Thread(target=reader, daemon=True)
    t.start()
    with open(output, "w") as fh:
        fh.write(airr_header() + "\n")
        while True:
            chunk = chunks.get()
            if chunk is None:
                break
            recs = _annotate_chunk(chunk, ref, target_db, seqtype,
                                   threads=threads, sensitivity=sensitivity,
                                   mm_strand=mm_strand, map_d=map_d)
            fh.write(format_rows(recs))
    t.join()
    return output
