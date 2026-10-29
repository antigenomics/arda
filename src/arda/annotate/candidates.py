"""Batched preliminary segment alignment evidence for probabilistic consumers."""
from __future__ import annotations

import os
from itertools import islice
from pathlib import Path
from tempfile import TemporaryDirectory

from .. import _markup, mmseqs
from .._threads import import_polars, resolve_threads
from .mapper import _cached_segment_db
from .reference import load_reference


pl = import_polars()

_QUERY_SCHEMA = {"sequence_id": pl.String, "sequence": pl.String}
_INT_COLS = {"qstart", "qend", "tstart", "tend", "qlen", "tlen", "alnlen", "mismatch", "gapopen"}
_TEXT_COLS = {"query", "target", "cigar", "qaln", "taln"}
_MM_SCHEMA = {c: pl.String if c in _TEXT_COLS else pl.Int64 if c in _INT_COLS else pl.Float64
              for c in mmseqs.DEFAULT_FORMAT_OUTPUT.split(",")}
_SCHEMA = {("sequence_id" if c == "query" else c): t for c, t in _MM_SCHEMA.items()} | {
    "segment": pl.String, "locus": pl.String, "alleles": pl.List(pl.String),
    "source": pl.String, "score": pl.Int64}


def alignment_candidates(records, *, organism="human", threads=0, max_seqs=300,
                         evalue=1e-3, sensitivity=7.0, batch_size=10000):
    """Return independent nucleotide V/J/D alignment candidates for each input.

    These are preliminary evidence, not assignments or recombination boundaries.
    Returns flat Polars ``queries`` and ``alignments`` frames plus ``limits``.
    Calls are bounded to ``batch_size`` inputs; iterate batches for large cohorts.
    Coordinates are 1-based closed in the input query and segment germline.
    V/J uses one batched forward-strand MMseqs segment search; callers supply
    coding-strand sequences. ``max_seqs`` caps prefilter targets, not alleles.
    MMseqs returns one local alignment per target subject to ``evalue``.
    D retains each allele's best positive (+1 match/-1 mismatch) ungapped span
    per allele, without an E-value gate or greedy D/D2 selection. Local
    endpoints are evidence only; inference must reconsider trimming and N/P.
    ``threads`` budgets the sequential MMseqs and native D phases. D uses
    contiguous query spans; batches below 32 queries run serially.
    ``ARDA_MMSEQS_SEARCH_OPTS`` must be unset; search settings are explicit here.
    No claim of exhaustive alignments or recombination scenarios is made.
    """
    if os.environ.get("ARDA_MMSEQS_SEARCH_OPTS", "").strip():
        raise ValueError("ARDA_MMSEQS_SEARCH_OPTS is unsupported by alignment_candidates; use explicit arguments")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    records = list(islice(records, batch_size + 1))
    if len(records) > batch_size:
        raise ValueError("supply at most batch_size records per bounded call")
    if max_seqs < 1 or evalue <= 0:
        raise ValueError("max_seqs and evalue must be positive")
    ids = [name for name, _ in records]
    if len(set(ids)) != len(ids):
        raise ValueError("sequence IDs must be unique")
    for name, seq in records:
        if not name or any(c.isspace() for c in name) or not seq or set(seq.upper()) - set("ACGTN"):
            raise ValueError("expected nonempty whitespace-free IDs and nucleotide sequences")
    records = [(name, seq.upper()) for name, seq in records]
    threads = resolve_threads(threads)
    if pl.thread_pool_size() > threads:
        raise ValueError(f"Polars is already initialized with {pl.thread_pool_size()} threads, "
                         f"exceeding the CPU budget {threads}. Start a fresh process with "
                         "POLARS_MAX_THREADS no larger than the requested budget before importing Polars.")
    limits = dict(max_seqs=max_seqs, evalue=evalue, sensitivity=sensitivity, kmer=13,
                  strand="forward", rejected_reverse_hits=0, d_threads=1 if len(records) < 32 else min(threads, len(records)), d_serial_below=32, vj_spans="one per target", d_spans="best positive span per allele", batch_size=batch_size)
    queries = pl.DataFrame({"sequence_id": ids, "sequence": [seq for _, seq in records]}, schema=_QUERY_SCHEMA)
    hits = pl.DataFrame(schema=_SCHEMA)
    if not records:
        return dict(queries=queries, alignments=hits, limits=limits)
    ref = load_reference(organism, "nt")
    db = _cached_segment_db(ref, organism)
    if db is None:
        raise RuntimeError("segment reference unavailable for candidate search")
    with TemporaryDirectory(prefix="arda_candidates_") as tmp:
        tmp = Path(tmp)
        fasta = tmp / "query.fasta"
        fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        qdb = mmseqs.createdb(fasta, tmp / "query", dbtype=2)
        result = mmseqs.search(qdb, db, tmp / "result", tmp / "work",
                               search_type=3, threads=threads, max_seqs=max_seqs,
                               evalue=evalue, sensitivity=sensitivity, kmer=13,
                               extra=["--strand", "1"])
        tsv = tmp / "hits.tsv"
        mmseqs.convertalis(qdb, db, result, tsv, threads=threads, search_type=3)
        if tsv.stat().st_size:
            raw = pl.read_csv(tsv, separator="\t", has_header=False, schema=_MM_SCHEMA,
                              quote_char=None).rename({"query": "sequence_id"})
            raw = raw.join(queries.select("sequence_id", pl.col("sequence").str.len_chars().alias("expected_qlen")),
                           on="sequence_id", how="left", validate="m:1")
            invalid = (pl.col("expected_qlen").is_null() | (pl.col("qlen") != pl.col("expected_qlen")) |
                       (pl.col("tlen") < 1) |
                       (pl.col("qstart") < 1) | (pl.col("qend") < 1) |
                       (pl.col("qstart") > pl.col("qlen")) | (pl.col("qend") > pl.col("qlen")) |
                       (pl.col("tstart") < 1) | (pl.col("tend") < 1) |
                       (pl.col("tstart") > pl.col("tlen")) | (pl.col("tend") > pl.col("tlen")))
            if raw.filter(invalid.fill_null(True)).height:
                raise RuntimeError("invalid MMseqs query ID, length or alignment bounds")
            reverse = (pl.col("qstart") > pl.col("qend")) | (pl.col("tstart") > pl.col("tend"))
            limits["rejected_reverse_hits"] = raw.filter(reverse).height
            raw = raw.filter(~reverse).drop("expected_qlen")
            metadata = pl.DataFrame([
                dict(target=target, segment=target[0], locus=entry.locus,
                     alleles=(entry.v_call if target.startswith("V|") else entry.j_call).split(","))
                for target, entry in ref.entries.items() if target.startswith(("V|", "J|"))])
            raw = raw.filter(pl.col("target").str.starts_with("V|") | pl.col("target").str.starts_with("J|"))
            hits = raw.join(metadata, on="target", how="left", validate="m:1").with_columns(
                pl.lit("mmseqs").alias("source"), pl.lit(None, dtype=pl.Int64).alias("score"))
            if hits["locus"].null_count():
                raise RuntimeError("segment targets missing reference metadata")
    ds = [(locus, allele, seq) for locus, genes in ref.d_germlines.items() for allele, seq in genes]
    rows = _markup.d_candidates_batch([seq for _, seq in records], [seq for _, _, seq in ds], threads=threads)
    if rows:
        d = pl.DataFrame(rows, schema={c: pl.Int64 for c in
            ("qi", "di", "score", "qstart", "qend", "tstart", "tend", "mismatch")}, orient="row")
        dref = pl.DataFrame({"locus": [x[0] for x in ds], "target": [x[1] for x in ds],
                             "germline": [x[2] for x in ds]}).with_row_index("di")
        d = d.join(queries.with_row_index("qi"), on="qi").join(dref, on="di").with_columns(
            pl.col("sequence").str.slice(pl.col("qstart")-1, pl.col("qend")-pl.col("qstart")+1).alias("qaln"),
            pl.col("germline").str.slice(pl.col("tstart")-1, pl.col("tend")-pl.col("tstart")+1).alias("taln"),
            pl.concat_list("target").alias("alleles"), pl.lit("D").alias("segment"),
            pl.lit("d_local").alias("source"),
            (pl.col("qend")-pl.col("qstart")+1).alias("alnlen"),
            pl.col("sequence").str.len_chars().cast(pl.Int64).alias("qlen"),
            pl.col("germline").str.len_chars().cast(pl.Int64).alias("tlen"),
            pl.lit(0, dtype=pl.Int64).alias("gapopen"),
        ).with_columns(pl.concat_str(pl.col("alnlen"), pl.lit("M")).alias("cigar"))
        for column, dtype in _SCHEMA.items():
            if column not in d.columns:
                d = d.with_columns(pl.lit(None, dtype=dtype).alias(column))
        hits = pl.concat([hits.select(list(_SCHEMA)), d.select(list(_SCHEMA))], how="vertical")
    return dict(queries=queries, alignments=hits.select(list(_SCHEMA)).sort(
        "sequence_id", "segment", "target", "qstart", "tstart", "qend", "tend"), limits=limits)
