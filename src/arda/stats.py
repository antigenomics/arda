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

⛔ Long, not wide, and deliberately: the metric set differs per scope (a gene has no junction
length, a chain has no allele frequency), so a wide table would be mostly empty cells, and the
one thing a QC table must support is ``grep`` / ``join`` / a per-metric plot across samples. One
value per cell, one row per fact -- no ``134/62`` hybrids, and integers stay integers.

⛔ **The chimera, non-functional and stop-codon counts are FLAGS, not filters.** Nothing here
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

__all__ = ["collect", "write_stats", "STATS_COLUMNS", "ALLELE_MIN_FRAC", "ALLELE_MIN_READS"]

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

    ⛔ The two quality columns arda writes are different encodings -- ``junction_quality`` is the
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


def _flatten_report(rep: dict, out: list[tuple], section: str = "") -> None:
    """The run report, verbatim, as ``run / <section> / <field> / <value>`` rows.

    Flattened wholesale rather than through a whitelist: the report gains a field every couple of
    releases (`fast_fraction`, `prefilter_stats`, `reads_assigned`), and a whitelist is how a new
    one silently fails to appear in the QC table for four releases.
    """
    for k, v in sorted(rep.items()):
        if isinstance(v, dict):
            # ⛔ Only the three STAGE names become a `key`; every other nested dict (`per_locus`,
            # `segment_search`, `prefilter_stats`, `reference`) folds into the metric name. A bare
            # `map --report` JSON has no stage wrapper, so recursing on shape instead would put
            # `per_locus` in the key column there and in the metric column for a merged report --
            # the same number under two addresses depending on which file it came from.
            if not section and k in ("map", "correct", "assemble"):
                _flatten_report(v, out, k)
            else:
                for sk, sv in sorted(v.items()):
                    if not isinstance(sv, (dict, list, tuple)):
                        out.append(("run", section, f"{k}.{sk}", _fmt(sv)))
        elif not isinstance(v, (list, tuple)):
            out.append(("run", section, k, _fmt(v)))


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
    })
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
        for label, wanted in (("reads", "reads"), ("clonotypes", "clonotypes")):
            if metric == wanted:
                seen.setdefault((seg, label), set()).add(key)
                if int(value) > 1:
                    multi.setdefault((seg, label), set()).add(key)

    for seg in ("V", "J"):
        total = len({g for (loc, s), genes in universe.items() if s == seg for g in genes})
        prefix = "v" if seg == "V" else "j"
        out.append(("sample", "", f"{prefix}_genes_reference", _fmt(total)))
        for label in ("reads", "clonotypes"):
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
            report: str | Path | None = None, r1: str | Path | None = None,
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
        r1, r2: the input FASTQs. Used ONLY for their size on disk and for whether the library is
            paired -- neither is recoverable from the AIRR, which holds the mapped subset.
    """
    out: list[tuple] = []
    if report is not None and Path(report).exists():
        _flatten_report(json.loads(Path(report).read_text()), out)
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
    _coverage_rows(out, organism)
    return out


def write_stats(rows: list[tuple], output: str | Path) -> int:
    """Write ``rows`` as the QC TSV. Returns the row count."""
    pl.DataFrame(rows, schema=STATS_COLUMNS, orient="row").write_csv(
        output, separator="\t", quote_style="never")
    return len(rows)
