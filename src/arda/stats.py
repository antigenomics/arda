"""Run QC: one long-format TSV describing a library, its reads and its clonotypes.

``arda stats`` reads the artifacts a run already wrote -- the Stage-1 AIRR, the clonotype table,
the ``.arda.json`` report -- and emits every number an operator needs to decide whether a sample
is usable, **without re-reading the FASTQ**. It adds no alignment and no reference lookup beyond
the germline gene list.

**Long format, four columns**, ``scope / key / metric / value``. The scopes, and what keys them:

``run``
    keyed by stage (``map`` / ``correct`` / ``assemble``) -- the run report, flattened verbatim
``sample``
    unkeyed -- library-wide totals, junction lengths and quality, SHM rate, gene coverage
``chain``
    keyed by locus (``TRB``, ``IGH``, ...) -- reads AND clonotypes
``v_gene`` / ``j_gene``
    keyed by gene (``TRBV19``) -- reads and clonotypes per germline gene
``allele_candidate``
    keyed by ``allele:mutation`` -- a recurrent, high-quality V mutation
``junction_aa_len`` / ``read_len`` / ``clone_size`` / ``isotype``
    the distributions, all keyed ``locus:bucket`` -- junction length in residues, aligned read
    length in 10-nt buckets, clonotype size in powers of two, constant-region class

Never: Long, not wide, and deliberately: the metric set differs per scope (a gene has no junction
length, a chain has no allele frequency), so a wide table would be mostly empty cells, and the
one thing a QC table must support is ``grep`` / ``join`` / a per-metric plot across samples. One
value per cell, one row per fact -- no ``134/62`` hybrids, and integers stay integers.

**Never: The chimera, non-functional and stop-codon counts are FLAGS, not filters.** Nothing here
removes a row from any output; ``stats`` only reads. See ``correct --flag-chimeras`` for why the
chimera signature cannot separate a true PCR artefact from two real clones sharing a prefix and
a suffix.

⚠ **Alleles vs SHM is a heuristic, and it is reported as one.** A mutation seen in most of an
allele's reads at high Phred is far more likely a germline the reference does not carry than
somatic hypermutation or a miscall -- but arda does not genotype, and ``allele_candidate`` is a
shortlist to look at, never a call. The thresholds are exposed (``--allele-min-frac``,
``--allele-min-reads``) precisely so the number can be re-derived rather than trusted.
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

__all__ = ["collect", "write_stats", "write_stats_json", "STATS_COLUMNS", "ALLELE_MIN_FRAC",
           "ALLELE_MIN_READS"]

STATS_COLUMNS = ["scope", "key", "metric", "value"]

#: A V mutation is a *candidate allele* when it is carried by at least this fraction of the reads
#: calling that allele AND by at least :data:`ALLELE_MIN_READS` of them. The fraction is what
#: separates germline from SHM (hypermutation is per-clone, so it does not reach half an allele's
#: reads); the count is what keeps a 2-read allele from producing a candidate off one read.
ALLELE_MIN_FRAC = 0.5
ALLELE_MIN_READS = 10

#: Reads pulled from the AIRR. Named explicitly because a Stage-1 AIRR has 83+ columns and
#: ``to_list()`` on all of them costs ~2.4 KB/row against ~0.4 KB for these -- the same reason
#: ``correct`` narrows before it materialises.
_READ_COLS = (
    "sequence_id", "locus", "v_call", "j_call", "junction", "junction_aa", "productive",
    "stop_codon", "vj_in_frame", "v_identity", "v_mutations", "j_mutations",
    "junction_quality", "v_mutation_quality", "j_mutation_quality", "mmseqs2_qlen",
    "rev_comp", "junction_completed_nt",
)
_CLONE_COLS = ("locus", "v_call", "j_call", "junction", "junction_aa", "duplicate_count",
               "consensus_count", "chimera_parents", "c_call")


def _fmt(value) -> str:
    """One value, one cell. Integers stay exact; floats get 6 significant digits."""
    if isinstance(value, bool):
        return str(int(value))
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return str(int(value)) if value.is_integer() else f"{value:.6g}"
    return "" if value is None else str(value)


def _read_present(path: Path, wanted: tuple[str, ...]) -> pl.DataFrame:
    """Read only the wanted columns that the file actually has.

    ``infer_schema_length=0`` keeps everything Utf8: a junction column of pure digits, an empty
    ``v_identity`` and a ``duplicate_count`` are all handled by an explicit cast at the point of
    use, rather than by whatever polars guesses from the first 100 rows of one sample.
    """
    header = path.open().readline().rstrip("\n").split("\t")
    cols = [c for c in wanted if c in header]
    if not cols:
        return pl.DataFrame()
    return pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None,
                       columns=cols)


def _phred(s: str) -> list[int]:
    """Parse a ``v_mutation_quality`` cell: comma-separated Phred INTEGERS, not a Phred+33 string.

    Never: The two quality columns arda writes are different encodings -- ``junction_quality`` is the
    raw Phred+33 characters (so it lines up byte-for-byte with ``junction``) and
    ``v_mutation_quality`` is comma-joined integers (there is no string to line up with). Reading
    one as the other yields plausible numbers off by 33 and a length that is silently wrong.
    """
    try:
        return [int(x) for x in s.split(",") if x != ""]
    except ValueError:
        return []


def _gene(call: str) -> str:
    """First allele of a (possibly tied) call, stripped to the gene: ``TRBV19*01,…`` -> ``TRBV19``."""
    return (call or "").split(",")[0].split("*")[0]


#: Metric names a MERGED map report uses for facts a single-file one names differently.
#: ``_merge_map_reports`` renames ``wall_seconds`` -> ``wall_seconds_max``/``_sum`` and
#: ``peak_rss_mb`` -> ``peak_rss_mb_max`` for a multi-read-group sample, so the same fact lives at
#: two addresses depending on how the sample was DELIVERED. Never: a cross-sample join on
#: ``run/map/wall_seconds`` then loses exactly the lane-split samples, with no error anywhere. The
#: canonical name is restored here, in the one place a report becomes QC rows, so ``.arda.json``
#: itself is untouched. ``wall_seconds_sum`` is a different fact and keeps its own name.
_REPORT_ALIASES = {"wall_seconds_max": "wall_seconds", "peak_rss_mb_max": "peak_rss_mb"}


def _dotted(prefix: str, value) -> list[tuple[str, object]]:
    """``{"a": {"b": 1}}`` -> ``[("a.b", 1)]``, to any depth.

    Never: one level of flattening is not enough. ``segment_search.reasons`` is a dict, so a
    one-level flatten dropped it silently -- and ``reasons`` is the ONLY evidence for why
    ``fast_fraction`` is low on a sample, which is the first question asked of a slow run.
    """
    if isinstance(value, dict):
        pairs: list[tuple[str, object]] = []
        for k, v in sorted(value.items()):
            pairs += _dotted(f"{prefix}.{k}" if prefix else str(k), v)
        return pairs
    return [(prefix, value)]


def _flatten_report(rep: dict, out: list[tuple], section: str = "") -> None:
    """The run report, verbatim, as ``run / <section> / <field> / <value>`` rows.

    Flattened wholesale rather than through a whitelist: the report gains a field every couple of
    releases (`fast_fraction`, `prefilter_stats`, `reads_assigned`), and a whitelist is how a new
    one silently fails to appear in the QC table for four releases.
    """
    for k, v in sorted(rep.items()):
        if isinstance(v, dict):
            # Never: Only the three STAGE names become a `key`; every other nested dict (`per_locus`,
            # `segment_search`, `prefilter_stats`, `reference`) folds into the metric name. A bare
            # `map --report` JSON has no stage wrapper, so recursing on shape instead would put
            # `per_locus` in the key column there and in the metric column for a merged report --
            # the same number under two addresses depending on which file it came from.
            if not section and k in ("map", "correct", "assemble"):
                _flatten_report(v, out, k)
            else:
                for name, sv in _dotted(k, v):
                    if not isinstance(sv, (list, tuple)):
                        out.append(("run", section, name, _fmt(sv)))
        elif isinstance(v, (list, tuple)):
            # A list of scalars is comma-joined, the way AIRR encodes a tied `v_call`. `input` is
            # a plain string for one read group and a LIST for several, so dropping lists here
            # deleted the provenance row exactly when a sample had more than one input file.
            if v and not any(isinstance(x, (dict, list, tuple)) for x in v):
                out.append(("run", section, k, ",".join(str(x) for x in v)))
        else:
            name = k
            if k in _REPORT_ALIASES and _REPORT_ALIASES[k] not in rep:
                name = _REPORT_ALIASES[k]
            out.append(("run", section, name, _fmt(v)))


def _ratio(rep: dict, out: list[tuple], name: str, num: tuple[str, str],
           den: tuple[str, str]) -> None:
    """Emit ``sample / "" / name`` = ``num / den``, or nothing if either is missing or zero."""
    a = (rep.get(num[0]) or {}).get(num[1])
    b = (rep.get(den[0]) or {}).get(den[1])
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and b:
        out.append(("sample", "", name, _fmt(a / b)))


def _yield_rows(rep: dict, out: list[tuple]) -> None:
    """The ratios that say where a library's reads ended up.

    Never: **a count does not compare across samples; a fraction does.** `reads_assigned` of
    28 means nothing next to another sample's 28 until you know one was handed 330 reads and the
    other 3.3 million -- and today those two numbers live in different stage blocks of the report,
    so the reader divides by hand or does not ask. A batch table compares METRICS, so the ratio
    has to be a metric.

    Only ratios whose denominator is unambiguous ship. `reads_incomplete` and `reads_low_quality`
    are counted inside `correct`, whose input is the junction-bearing subset rather than every
    mapped read, and there is no counter that names that subset (`correct.reads` is explicitly
    NOT the conservation invariant, `correct.py:150-155`). Rather than divide by a denominator
    that is nearly right, they stay counts.
    """
    # `reads_assigned` is the read-conservation invariant -- the sum of `duplicate_count` over the
    # clonotype table, i.e. every read that ended up inside a clonotype, assembly rescues included.
    _ratio(rep, out, "reads_used_fraction", ("correct", "reads_assigned"), ("map", "total_reads"))
    _ratio(rep, out, "reads_used_of_mapped_fraction",
           ("correct", "reads_assigned"), ("map", "mapped_reads"))
    _ratio(rep, out, "reads_per_clonotype_mean",
           ("correct", "reads_assigned"), ("correct", "clonotypes_out"))
    _ratio(rep, out, "contigs_complete_fraction",
           ("assemble", "contigs_complete"), ("assemble", "contigs"))

    # The unmapped ledger, as fractions of the reads arda was handed. Same argument as above and
    # the same denominator for every bucket, stated once: `prefilter_rejected` of 788 says nothing
    # next to another sample's 788, but 0.597 against a batch median of 0.61 says the prefilter is
    # behaving and 0.99 says it is eating the library.
    total = (rep.get("map") or {}).get("total_reads")
    for bucket, n in sorted(((rep.get("map") or {}).get("unmapped") or {}).items()):
        if bucket != "accounted" and isinstance(n, int) and isinstance(total, int) and total:
            out.append(("sample", "", f"unmapped_{bucket}_fraction", _fmt(n / total)))


def _gene_universe(organism: str) -> dict[tuple[str, str], set[str]]:
    """``(locus, segment) -> {gene}`` from the reference's own anchor table.

    ``cdr3_anchors.tsv`` is the right source: it is per ALLELE with its locus, segment and IMGT
    functionality already resolved, so gene coverage is measured against the germline set arda
    actually maps to rather than against a hand-kept list that drifts from it.
    """
    from .paths import vdj_dir

    path = vdj_dir() / organism / "cdr3_anchors.tsv"
    if not path.exists():
        return {}
    df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
    universe: dict[tuple[str, str], set[str]] = {}
    for locus, seg, allele in zip(df["locus"], df["segment"], df["allele"]):
        universe.setdefault((locus, seg), set()).add(allele.split("*")[0])
    return universe


def _hist(work: pl.DataFrame, scope: str, col: str, metric: str, out: list[tuple],
          *, weight: str | None = None) -> None:
    """A sparse distribution as ``scope / <locus>:<bucket> / metric / count`` rows.

    Never: min/max/mean cannot show a SHAPE, and shape is the whole diagnostic. A bimodal junction
    length is two primer sets in one tube; a read-length cliff is an adapter left on; a clone-size
    distribution with no singletons is a library that was over-amplified before it was sequenced.
    Each reads as an ordinary mean.

    Only OCCUPIED buckets get a row -- the same rule the ``v_gene``/``j_gene`` scopes already use.
    An amplicon occupies ~30 of ~200 possible junction lengths, so a row per empty bucket would be
    most of the table. The key is ``locus:bucket`` for every distribution scope, so one parser
    reads all of them.
    """
    agg = pl.len() if weight is None else pl.col(weight).sum()
    counts = (work.filter(pl.col(col) > 0).group_by(["locus", col])
              .agg(agg.alias("n")).sort(["locus", col]))
    for locus, bucket, n in counts.iter_rows():
        out.append((scope, f"{locus or '?'}:{bucket}", metric, _fmt(n)))


# ── read-level ────────────────────────────────────────────────────────────────────────────────

def _read_rows(df: pl.DataFrame, out: list[tuple]) -> None:
    """Per-chain read statistics, plus their sample-wide totals."""
    from .rnaseq.correct import _CANONICAL_AA

    n = df.height
    have = set(df.columns)
    jn = df["junction"] if "junction" in have else pl.Series([""] * n)
    jaa = df["junction_aa"] if "junction_aa" in have else pl.Series([""] * n)
    locus = (df["locus"] if "locus" in have else pl.Series([""] * n)).fill_null("")

    # A junction is TRUNCATED when it does not span both conserved anchors -- a distinct defect
    # from a stop codon or a frameshift, which `_COMPLETE` folds together and which the caller
    # here needs separated (a truncated junction is a short read; a stop codon is biology).
    spans = jaa.fill_null("").str.contains(_CANONICAL_AA)
    work = pl.DataFrame({
        "locus": locus,
        "_has_jn": jn.fill_null("").str.len_chars() > 0,
        "_spans": spans,
        "_jnt": jn.fill_null("").str.len_chars(),
        "_jaa": jaa.fill_null("").str.len_chars(),
        "_stop": jaa.fill_null("").str.contains(r"\*"),
        "_inframe": (jn.fill_null("").str.len_chars() % 3) == 0,
        "_v_gene": pl.Series([_gene(c) for c in (df["v_call"] if "v_call" in have
                                                 else [""] * n)]),
        "_j_gene": pl.Series([_gene(c) for c in (df["j_call"] if "j_call" in have
                                                 else [""] * n)]),
    })
    if "rev_comp" in have:
        work = work.with_columns(_rc=(df["rev_comp"] == "T"))
    if "junction_completed_nt" in have:
        work = work.with_columns(
            _completed=(df["junction_completed_nt"].cast(pl.Int64, strict=False).fill_null(0) > 0))
    if "mmseqs2_qlen" in have:
        # 10-nt buckets. The read length arda actually aligned, which is NOT `read_length_mean`
        # in the run report -- that one is measured over every read before mapping.
        # Never: via Float64. mmseqs writes this column as `407.0`, and a direct Int64 cast of
        # "407.0" is null under `strict=False` -- so every bucket came out 0 and the whole
        # distribution silently vanished rather than failing.
        work = work.with_columns(
            _qlen_bin=((df["mmseqs2_qlen"].cast(pl.Float64, strict=False).fill_null(0)
                        .cast(pl.Int64)) // 10) * 10)
    if "productive" in have:
        work = work.with_columns(_productive=(df["productive"] == "T"),
                                 _nonfunctional=(df["productive"] == "F"))
    if "v_identity" in have:
        work = work.with_columns(
            _vid=df["v_identity"].cast(pl.Float64, strict=False))
    for col, name in (("v_mutations", "_vmut"), ("j_mutations", "_jmut")):
        if col in have:
            s = df[col].fill_null("")
            work = work.with_columns(
                pl.when(s.str.len_chars() == 0).then(0)
                .otherwise(s.str.count_matches(",") + 1).alias(name))
    if "junction_quality" in have:
        # `ponytail:` a Python pass over the Phred strings -- polars has no ord() over Utf8.
        # `sum(map(ord, s))` is C-speed and the column only exists under `map --junction-quality`.
        q = df["junction_quality"].fill_null("").to_list()
        work = work.with_columns(
            pl.Series("_qmean", [sum(map(ord, s)) / len(s) - 33 if s else None for s in q]),
            pl.Series("_qmin", [min(map(ord, s)) - 33 if s else None for s in q]))

    #: metric name -> how to aggregate it. One table, so a chain row and the sample row are the
    #: same computation over a different frame -- they cannot disagree.
    aggs = [
        ("reads", pl.len()),
        ("reads_with_junction", pl.col("_has_jn").sum()),
        ("reads_truncated_junction", (pl.col("_has_jn") & ~pl.col("_spans")).sum()),
        ("reads_stop_codon", (pl.col("_has_jn") & pl.col("_stop")).sum()),
        ("reads_out_of_frame", (pl.col("_has_jn") & ~pl.col("_inframe")).sum()),
        ("junction_nt_min", pl.col("_jnt").filter(pl.col("_spans")).min()),
        ("junction_nt_max", pl.col("_jnt").filter(pl.col("_spans")).max()),
        ("junction_nt_mean", pl.col("_jnt").filter(pl.col("_spans")).mean()),
        ("junction_aa_min", pl.col("_jaa").filter(pl.col("_spans")).min()),
        ("junction_aa_max", pl.col("_jaa").filter(pl.col("_spans")).max()),
        ("junction_aa_mean", pl.col("_jaa").filter(pl.col("_spans")).mean()),
        ("v_genes_observed", pl.col("_v_gene").filter(pl.col("_v_gene") != "").n_unique()),
        ("j_genes_observed", pl.col("_j_gene").filter(pl.col("_j_gene") != "").n_unique()),
    ]
    if "_productive" in work.columns:
        aggs += [("reads_productive", pl.col("_productive").sum()),
                 ("reads_nonfunctional", pl.col("_nonfunctional").sum())]
    if "_vid" in work.columns:
        aggs += [("v_identity_mean", pl.col("_vid").mean()),
                 ("shm_rate", 1.0 - pl.col("_vid").mean())]
    for name in ("_vmut", "_jmut"):
        if name in work.columns:
            aggs.append((f"{name[1]}_mutations_per_read", pl.col(name).mean()))
    if "_qmean" in work.columns:
        aggs += [("junction_quality_mean", pl.col("_qmean").mean()),
                 ("junction_quality_min_mean", pl.col("_qmin").mean())]
    if "_rc" in work.columns:
        # Strand balance. A stranded protocol that lost its orientation, or an unstranded one read
        # as stranded, shows up here and in no other number arda writes.
        aggs.append(("reads_rev_comp", pl.col("_rc").sum()))
    if "_completed" in work.columns:
        aggs.append(("reads_junction_completed", pl.col("_completed").sum()))

    exprs = [e.alias(name) for name, e in aggs]
    for row in work.group_by("locus").agg(exprs).sort("locus").iter_rows(named=True):
        for name, _ in aggs:
            out.append(("chain", row["locus"] or "?", name, _fmt(row[name])))
    for name, value in zip([a for a, _ in aggs],
                           work.select(exprs).row(0)):
        out.append(("sample", "", name, _fmt(value)))

    # Per-gene reads. Only genes that were SEEN get a row -- the reference universe is reported as
    # a coverage fraction below, and one row per unobserved gene would be 90 % of this table.
    for seg, col in (("v_gene", "_v_gene"), ("j_gene", "_j_gene")):
        counts = (work.filter(pl.col(col) != "").group_by(col).agg(pl.len().alias("reads"))
                  .sort(col))
        for gene, reads in counts.iter_rows():
            out.append((seg, gene, "reads", _fmt(reads)))

    # Junction length over SPANNING reads only, the same population the min/max/mean above use --
    # a truncated junction has a length, and counting it would put a second, fake mode on the left.
    _hist(work.filter(pl.col("_spans")), "junction_aa_len", "_jaa", "reads", out)
    if "_qlen_bin" in work.columns:
        _hist(work, "read_len", "_qlen_bin", "reads", out)


def _allele_rows(df: pl.DataFrame, out: list[tuple], *, min_frac: float, min_reads: int) -> None:
    """Split the V mutation entries into candidate alleles and SHM, and score both on quality.

    A mutation list, a novel allele and a miscall are the same string in the AIRR. What separates
    them is FREQUENCY WITHIN THE ALLELE (a germline difference is in every read of that allele;
    hypermutation is per-clone) and PHRED (a miscall is not). So both are reported, per variant,
    and the classification is a threshold the caller can move.
    """
    if "v_call" not in df.columns or "v_mutations" not in df.columns:
        return
    calls = df["v_call"].fill_null("").to_list()
    muts = df["v_mutations"].fill_null("").to_list()
    quals = (df["v_mutation_quality"].fill_null("").to_list()
             if "v_mutation_quality" in df.columns else [""] * len(calls))

    per_allele: dict[str, int] = {}
    hits: dict[tuple[str, str], list[int]] = {}   # (allele, mutation) -> its Phred scores
    counts: dict[tuple[str, str], int] = {}
    for call, mut, qual in zip(calls, muts, quals):
        allele = call.split(",")[0]
        if not allele:
            continue
        per_allele[allele] = per_allele.get(allele, 0) + 1
        if not mut:
            continue
        entries = mut.split(",")
        scores = _phred(qual) if qual else []
        # A quality string that does not match the mutation list 1:1 describes different bases;
        # drop it rather than pair entry i with score i and report a number that is not one.
        if len(scores) != len(entries):
            scores = []
        for i, entry in enumerate(entries):
            key = (allele, entry)
            counts[key] = counts.get(key, 0) + 1
            if scores:
                hits.setdefault(key, []).append(scores[i])

    n_allele, n_shm = 0, 0
    q_allele: list[int] = []
    q_shm: list[int] = []
    for (allele, entry), reads in sorted(counts.items()):
        denom = per_allele.get(allele, 0)
        frac = reads / denom if denom else 0.0
        scores = hits.get((allele, entry), [])
        candidate = frac >= min_frac and reads >= min_reads
        if candidate:
            n_allele += 1
            q_allele += scores
            key = f"{allele}:{entry}"
            out.append(("allele_candidate", key, "reads", _fmt(reads)))
            out.append(("allele_candidate", key, "allele_reads", _fmt(denom)))
            out.append(("allele_candidate", key, "frequency", _fmt(frac)))
            if scores:
                out.append(("allele_candidate", key, "mean_quality",
                            _fmt(sum(scores) / len(scores))))
        else:
            n_shm += 1
            q_shm += scores
    out.append(("sample", "", "allele_candidates", _fmt(n_allele)))
    out.append(("sample", "", "shm_variants", _fmt(n_shm)))
    out.append(("sample", "", "allele_candidate_min_frac", _fmt(min_frac)))
    out.append(("sample", "", "allele_candidate_min_reads", _fmt(min_reads)))
    if q_allele:
        out.append(("sample", "", "allele_candidate_mean_quality",
                    _fmt(sum(q_allele) / len(q_allele))))
    if q_shm:
        out.append(("sample", "", "shm_variant_mean_quality", _fmt(sum(q_shm) / len(q_shm))))


# ── clonotype-level ───────────────────────────────────────────────────────────────────────────

def _clone_rows(df: pl.DataFrame, out: list[tuple]) -> None:
    from .rnaseq.correct import _CANONICAL_AA

    n = df.height
    have = set(df.columns)
    jaa = (df["junction_aa"] if "junction_aa" in have else pl.Series([""] * n)).fill_null("")
    jn = (df["junction"] if "junction" in have else pl.Series([""] * n)).fill_null("")
    work = pl.DataFrame({
        "locus": (df["locus"] if "locus" in have else pl.Series([""] * n)).fill_null(""),
        "_dup": (df["duplicate_count"].cast(pl.Int64, strict=False)
                 if "duplicate_count" in have else pl.Series([0] * n, dtype=pl.Int64)),
        "_spans": jaa.str.contains(_CANONICAL_AA),
        "_stop": jaa.str.contains(r"\*"),
        "_inframe": (jn.str.len_chars() % 3) == 0,
        "_jnt": jn.str.len_chars(),
        "_jaa": jaa.str.len_chars(),
        "_chimera": ((df["chimera_parents"].fill_null("").str.len_chars() > 0)
                     if "chimera_parents" in have
                     else pl.Series([False] * n, dtype=pl.Boolean)),
        "_v_gene": pl.Series([_gene(c) for c in (df["v_call"] if "v_call" in have
                                                 else [""] * n)]),
        "_j_gene": pl.Series([_gene(c) for c in (df["j_call"] if "j_call" in have
                                                 else [""] * n)]),
        "_isotype": (df["c_call"] if "c_call" in have else pl.Series([""] * n)).fill_null(""),
    })
    # Power-of-two buckets, keyed by the bucket's LOWER BOUND (1, 2, 4, 8, ...). Clone sizes span
    # five orders of magnitude in a real library, so linear buckets are one occupied row and a
    # thousand empty ones.
    work = work.with_columns(
        _size_bin=pl.when(pl.col("_dup") > 0)
        .then((2 ** pl.col("_dup").log(2).floor()).cast(pl.Int64)).otherwise(0))
    aggs = [
        ("clonotypes", pl.len()),
        ("clonotype_reads", pl.col("_dup").sum()),
        ("clonotypes_truncated_junction", (~pl.col("_spans")).sum()),
        ("clonotypes_stop_codon", pl.col("_stop").sum()),
        ("clonotypes_out_of_frame", (~pl.col("_inframe")).sum()),
        ("clonotypes_chimeric", pl.col("_chimera").sum()),
        ("chimeric_reads", pl.col("_dup").filter(pl.col("_chimera")).sum()),
        ("clonotype_junction_nt_min", pl.col("_jnt").min()),
        ("clonotype_junction_nt_max", pl.col("_jnt").max()),
        ("clonotype_junction_nt_mean", pl.col("_jnt").mean()),
        ("clonotype_junction_aa_min", pl.col("_jaa").min()),
        ("clonotype_junction_aa_max", pl.col("_jaa").max()),
        ("clonotype_junction_aa_mean", pl.col("_jaa").mean()),
        ("clonotype_v_genes_observed",
         pl.col("_v_gene").filter(pl.col("_v_gene") != "").n_unique()),
        ("clonotype_j_genes_observed",
         pl.col("_j_gene").filter(pl.col("_j_gene") != "").n_unique()),
    ]
    exprs = [e.alias(name) for name, e in aggs]
    for row in work.group_by("locus").agg(exprs).sort("locus").iter_rows(named=True):
        for name, _ in aggs:
            out.append(("chain", row["locus"] or "?", name, _fmt(row[name])))
    for name, value in zip([a for a, _ in aggs], work.select(exprs).row(0)):
        out.append(("sample", "", name, _fmt(value)))

    for seg, col in (("v_gene", "_v_gene"), ("j_gene", "_j_gene")):
        counts = (work.filter(pl.col(col) != "")
                  .group_by(col).agg(pl.len().alias("clonotypes"),
                                     pl.col("_dup").sum().alias("reads_in_clonotypes"))
                  .sort(col))
        for gene, clones, reads in counts.iter_rows():
            out.append((seg, gene, "clonotypes", _fmt(clones)))
            out.append((seg, gene, "reads_in_clonotypes", _fmt(reads)))

    _hist(work, "clone_size", "_size_bin", "clonotypes", out)
    _hist(work, "clone_size", "_size_bin", "reads", out, weight="_dup")
    _hist(work.filter(pl.col("_spans")), "junction_aa_len", "_jaa", "clonotypes", out)

    # Isotype, from the clonotype's dominant `c_call`. Keyed like every other distribution scope
    # so one parser reads them all; a clonotype with no constant-region evidence gets no row.
    iso = (work.filter(pl.col("_isotype") != "").group_by(["locus", "_isotype"])
           .agg(pl.len().alias("clonotypes"), pl.col("_dup").sum().alias("reads"))
           .sort(["locus", "_isotype"]))
    for locus, isotype, clones, reads in iso.iter_rows():
        out.append(("isotype", f"{locus or '?'}:{isotype}", "clonotypes", _fmt(clones)))
        out.append(("isotype", f"{locus or '?'}:{isotype}", "reads", _fmt(reads)))


# ── single cell ───────────────────────────────────────────────────────────────────────────────

#: What `arda cells` writes, and what QC reads back out of it.
_CHAIN_COLS = ("cell_id", "locus", "v_call", "j_call", "junction_aa", "molecules", "reads",
               "status")


def _cell_rows(prefix: Path, out: list[tuple]) -> None:
    """``arda cells`` output as the SAME scopes a bulk run produces.

    Never: single cell was the one mode that wrote no QC table at all, and its ``.arda.json``
    shares no key with the bulk one -- so a cohort mixing the two had nothing to join on. The
    scopes are deliberately the bulk ones (``sample`` / ``chain`` / ``v_gene`` / ``j_gene`` /
    ``junction_aa_len``), because "which loci did this sample yield, and at what junction lengths"
    is the same question either way. What is genuinely single-cell -- cells, molecules, the knee,
    the contig N50 -- comes through the ``run`` scope from its own report, where it cannot be
    mistaken for a read count.

    ``pairing_rate`` and ``doublet_rate`` are derived here rather than left to the reader: they
    are the AIRR Community's chain-pairing QC, they are the two numbers a batch is compared on,
    and ``cell_summary`` has already assigned every cell the status they count.
    """
    report = prefix.parent / f"{prefix.name}.arda.json"
    if report.exists():
        _flatten_report(json.loads(report.read_text()), out)
        rep = json.loads(report.read_text())
        cells = rep.get("cells") or 0
        if cells:
            for status in ("paired", "doublet_candidate"):
                n = rep.get(f"cells_{status}")
                if n is not None:
                    name = "pairing_rate" if status == "paired" else "doublet_rate"
                    out.append(("sample", "", name, _fmt(n / cells)))
        placed, seen = rep.get("molecules_placed"), rep.get("molecules_in")
        if placed is not None and seen:
            out.append(("sample", "", "molecules_placed_fraction", _fmt(placed / seen)))

    path = prefix.parent / f"{prefix.name}.chains.tsv"
    if not path.exists():
        return
    df = _read_present(path, _CHAIN_COLS)
    if not df.height:
        return
    have = set(df.columns)
    n = df.height
    # `extra` chains are the ones the extra-chain gate REJECTED. Counting them as yield would
    # report ambient contamination as signal -- the whole point of the gate.
    work = pl.DataFrame({
        "locus": (df["locus"] if "locus" in have else pl.Series([""] * n)).fill_null(""),
        "cell_id": (df["cell_id"] if "cell_id" in have else pl.Series([""] * n)).fill_null(""),
        "_called": ((df["status"] != "extra") if "status" in have
                    else pl.Series([True] * n, dtype=pl.Boolean)),
        "_mol": (df["molecules"].cast(pl.Int64, strict=False).fill_null(0)
                 if "molecules" in have else pl.Series([0] * n, dtype=pl.Int64)),
        "_reads": (df["reads"].cast(pl.Int64, strict=False).fill_null(0)
                   if "reads" in have else pl.Series([0] * n, dtype=pl.Int64)),
        "_jaa": (df["junction_aa"] if "junction_aa" in have
                 else pl.Series([""] * n)).fill_null("").str.len_chars(),
        "_v_gene": pl.Series([_gene(c) for c in (df["v_call"] if "v_call" in have
                                                 else [""] * n)]),
        "_j_gene": pl.Series([_gene(c) for c in (df["j_call"] if "j_call" in have
                                                 else [""] * n)]),
    }).filter(pl.col("_called"))

    aggs = [("chains", pl.len()),
            ("cells", pl.col("cell_id").n_unique()),
            ("molecules", pl.col("_mol").sum()),
            ("reads", pl.col("_reads").sum()),
            ("junction_aa_min", pl.col("_jaa").filter(pl.col("_jaa") > 0).min()),
            ("junction_aa_max", pl.col("_jaa").filter(pl.col("_jaa") > 0).max()),
            ("junction_aa_mean", pl.col("_jaa").filter(pl.col("_jaa") > 0).mean())]
    exprs = [e.alias(name) for name, e in aggs]
    for row in work.group_by("locus").agg(exprs).sort("locus").iter_rows(named=True):
        for name, _ in aggs:
            out.append(("chain", row["locus"] or "?", name, _fmt(row[name])))
    for name, value in zip([a for a, _ in aggs], work.select(exprs).row(0)):
        out.append(("sample", "", name, _fmt(value)))

    for seg, col in (("v_gene", "_v_gene"), ("j_gene", "_j_gene")):
        counts = (work.filter(pl.col(col) != "").group_by(col)
                  .agg(pl.len().alias("chains")).sort(col))
        for gene, chains in counts.iter_rows():
            out.append((seg, gene, "chains", _fmt(chains)))

    _hist(work, "junction_aa_len", "_jaa", "chains", out)
    # Molecules per chain, powers of two -- the same axis `scplot`'s chain_support panel draws.
    work = work.with_columns(
        _mol_bin=pl.when(pl.col("_mol") > 0)
        .then((2 ** pl.col("_mol").log(2).floor()).cast(pl.Int64)).otherwise(0))
    _hist(work, "chain_support", "_mol_bin", "chains", out)


def _coverage_rows(out: list[tuple], organism: str) -> None:
    """``% of reference genes seen``, per locus and sample-wide, for reads and for clonotypes.

    Computed from the rows already in ``out`` rather than from a third pass over the frames: the
    per-gene counts are there, the per-locus totals are there, and re-deriving them would be a
    second definition of "observed" that can drift from the one the chain rows used.
    """
    universe = _gene_universe(organism)
    if not universe:
        return
    seen: dict[tuple[str, str], set[str]] = {}   # (segment, "reads"/"clonotypes") -> genes
    multi: dict[tuple[str, str], set[str]] = {}
    for scope, key, metric, value in out:
        if scope not in ("v_gene", "j_gene"):
            continue
        seg = "V" if scope == "v_gene" else "J"
        for label, wanted in (("reads", "reads"), ("clonotypes", "clonotypes"),
                              ("chains", "chains")):
            if metric == wanted:
                seen.setdefault((seg, label), set()).add(key)
                if int(value) > 1:
                    multi.setdefault((seg, label), set()).add(key)

    for seg in ("V", "J"):
        total = len({g for (loc, s), genes in universe.items() if s == seg for g in genes})
        prefix = "v" if seg == "V" else "j"
        out.append(("sample", "", f"{prefix}_genes_reference", _fmt(total)))
        for label in ("reads", "clonotypes", "chains"):
            got = seen.get((seg, label), set())
            if not got:
                continue
            more = multi.get((seg, label), set())
            out.append(("sample", "", f"{prefix}_gene_coverage_{label}",
                        _fmt(len(got) / total if total else 0.0)))
            out.append(("sample", "", f"{prefix}_gene_coverage_{label}_multi",
                        _fmt(len(more) / total if total else 0.0)))


# ── entry points ──────────────────────────────────────────────────────────────────────────────

def collect(*, airr: str | Path | None = None, clones: str | Path | None = None,
            report: str | Path | None = None, cells: str | Path | None = None,
            r1: str | Path | None = None,
            r2: str | Path | None = None, organism: str = "human",
            allele_min_frac: float = ALLELE_MIN_FRAC,
            allele_min_reads: int = ALLELE_MIN_READS) -> list[tuple]:
    """Every statistic arda can derive from a finished run, as ``(scope, key, metric, value)``.

    Every input is optional and each contributes its own scopes, so this works on a bare
    ``arda annotate`` output as well as on a full ``arda rnaseq`` run directory.

    Args:
        airr: Stage-1 (or ``annotate``) AIRR TSV -> the ``chain`` read rows, ``*_gene`` reads,
            and ``allele_candidate``.
        clones: clonotype table -> the ``chain`` clonotype rows and ``*_gene`` clonotypes.
        report: ``<prefix>.arda.json`` (or a bare ``--report`` JSON) -> the ``run`` scope, which
            is where total/mapped reads, threads, wall time and peak RSS come from.
        cells: an ``arda cells`` output PREFIX (not a file) -> its report plus ``.chains.tsv``,
            as the same ``sample`` / ``chain`` / ``*_gene`` / ``junction_aa_len`` scopes a bulk
            run produces, so one batch table holds both kinds of sample.
        r1, r2: the input FASTQs. Used ONLY for their size on disk and for whether the library is
            paired -- neither is recoverable from the AIRR, which holds the mapped subset.
    """
    out: list[tuple] = []
    if report is not None and Path(report).exists():
        rep = json.loads(Path(report).read_text())
        _flatten_report(rep, out)
        _yield_rows(rep, out)
    if r1 is not None:
        paths = [Path(p) for p in (r1, r2) if p is not None]
        out.append(("sample", "", "paired", _fmt(r2 is not None)))
        out.append(("sample", "", "input_files", _fmt(len(paths))))
        out.append(("sample", "", "input_bytes",
                    _fmt(sum(p.stat().st_size for p in paths if p.exists()))))
    if airr is not None:
        df = _read_present(Path(airr), _READ_COLS)
        if df.height:
            _read_rows(df, out)
            _allele_rows(df, out, min_frac=allele_min_frac, min_reads=allele_min_reads)
    if clones is not None:
        df = _read_present(Path(clones), _CLONE_COLS)
        if df.height:
            _clone_rows(df, out)
    if cells is not None:
        _cell_rows(Path(cells), out)
    _coverage_rows(out, organism)
    # Never: a metric with no input is OMITTED, not blank. An aggregation over an empty filtered
    # set returns None, which `_fmt` renders as "" -- and a reader casting the column silently
    # gets 0, so `chain/IGH/junction_nt_min` read as "the shortest IGH junction was 0 nt" when
    # the truth is "no IGH read spanned both anchors". Dropping the row is what `docs/qc.rst`
    # already promises and what a cross-sample join needs to see as null.
    return [r for r in out if r[3] != ""]


def write_stats(rows: list[tuple], output: str | Path) -> int:
    """Write ``rows`` as the QC TSV. Returns the row count."""
    pl.DataFrame(rows, schema=STATS_COLUMNS, orient="row").write_csv(
        output, separator="\t", quote_style="never")
    return len(rows)


def _typed(value: str):
    """``"453"`` -> 453, ``"0.1736"`` -> 0.1736, ``"TRBV19"`` -> ``"TRBV19"``."""
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def write_stats_json(rows: list[tuple], output: str | Path, *, arda_version: str = "") -> int:
    """The same rows as nested, TYPED JSON: ``{scope: {key: {metric: value}}}``.

    The TSV is what an operator greps and what a shell pipeline joins; this is what a program
    reads, and what ``arda qc report`` inlines into its HTML. Both are written from the SAME
    ``rows`` list rather than from two passes over the data, so they cannot disagree -- the cost
    is re-parsing a number that was just formatted, which is nothing against the frames it came
    from.

    Values are typed here and stringified there for the same reason the TSV drops blank rows: a
    consumer must be able to tell "no IGH read spanned both anchors" from "the minimum was 0".
    """
    stats: dict[str, dict[str, dict[str, object]]] = {}
    for scope, key, metric, value in rows:
        stats.setdefault(scope, {}).setdefault(key, {})[metric] = _typed(value)
    doc = {"arda_version": arda_version or _version(), "rows": len(rows), "stats": stats}
    Path(output).write_text(json.dumps(doc, indent=2, sort_keys=True))
    return len(rows)


def _version() -> str:
    from . import __version__
    return __version__
