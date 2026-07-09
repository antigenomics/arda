"""Stage 2 — CDR3 error correction (parent:child count ratio).

Collapses sequencing-error CDR3 variants onto their parent clonotype, a Python port
of vdjtools' ``Corrector`` (``com.antigenomics.vdjtools.preprocess.Corrector``) on top
of :mod:`seqtree` neighbour search (the fast substitution-bounded index that plays the
role of milib's ``SequenceTreeMap``). Rule (vdjtools defaults ``maxMismatches=2``,
``ratio=0.05``): a clonotype ``C`` is a **child** of a neighbour ``P`` that differs by
``m`` substitutions iff ``count[P] >= count[C] * (1/ratio)**m`` (i.e. the parent is at
least ``20**m`` times more abundant). Children route their reads to the parent; chains
(``C -> B -> A``) collapse to the ultimate ancestor. Counts strictly increase along
parent pointers, so there are no cycles.

seqtree is an optional dependency (``pip install 'arda-mapper[rnaseq]'``).
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import polars as pl

__all__ = ["correct_airr", "CorrectReport"]

_DNA = frozenset("ACGT")

# The generic heavy constant `isotype_class` returns when a read's C hit spans classes (IGHG1,IGHM
# -> IGHC): isotype unresolved. Deprioritised when aggregating a clonotype's dominant isotype.
_GENERIC_ISOTYPE = frozenset({"IGHC"})

# A clonotype requires a COMPLETE junction. Stage 1 reports a junction even when the read does
# not span it (see ``annotate.transfer``), so a raw per-read junction is a truncated fragment
# whenever the CDR3 runs off the end of the read. Aggregating those as clonotypes inflates the
# repertoire with prefixes of real junctions -- measured on PRJNA371303 RNA-seq (100 bp reads):
# 42 % of IGH and 29 % of TRB "clonotypes" were truncated, and ~10 % carried a stop codon.
# The effect is chain-dependent (IGH/TRB junctions are long; IGK/IGL fit in a 100 bp read),
# so it silently distorts chain composition, not just the total.
#
# `productive` alone is NOT sufficient: it flags stops and frameshifts but says nothing about
# truncation (640/5393 productive rows had no [FW] anchor). Require the canonical Cys104...
# [FW]118 anchors as well. Assembly-based tools only ever emit complete junctions,
# so this is also what makes arda's clonotype table comparable to theirs.
_CANONICAL_AA = r"^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$"
_COMPLETE = (
    pl.col("junction_aa").is_not_null()
    & pl.col("junction_aa").str.contains(_CANONICAL_AA)  # spans both anchors
    & ~pl.col("junction_aa").str.contains(r"\*")         # no stop codon
    & ~pl.col("junction_aa").str.contains("_")           # no frameshift-inserted N
    & (pl.col("junction").str.len_chars() % 3 == 0)      # in frame
)


@dataclass
class CorrectReport:
    clonotypes_in: int = 0
    clonotypes_out: int = 0
    reads: int = 0
    collapsed: int = 0  # clonotypes absorbed into a parent
    reads_with_junction: int = 0   # Stage-1 reads carrying any junction
    reads_incomplete: int = 0      # ...of which dropped as truncated/out-of-frame/stop

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["collapse_fraction"] = self.collapsed / self.clonotypes_in if self.clonotypes_in else 0.0
        return d


def _parents(junctions: list[str], counts: list[int], v: list[str], j: list[str],
             *, max_mismatches: int, ratio: float, require_vj: bool) -> list[int | None]:
    """For each clonotype, the strongest neighbour that is its parent (or None)."""
    import seqtree

    n = len(junctions)
    parent: list[int | None] = [None] * n
    safe = [i for i, s in enumerate(junctions) if s and set(s) <= _DNA]
    if not safe:
        return parent
    index = seqtree.Index.build([junctions[i] for i in safe], alphabet="nt")
    params = seqtree.SearchParams(
        max_subs=max_mismatches, max_total_edits=max_mismatches, engine="seqtm")
    inv_ratio = 1.0 / ratio
    for local_q, ci in enumerate(safe):
        best, best_count = None, -1
        for hit in index.search(junctions[ci], params):
            nj = safe[hit.ref_id]
            m = hit.n_subs
            if nj == ci or m == 0:
                continue
            if require_vj and (v[nj] != v[ci] or j[nj] != j[ci]):
                continue
            if counts[nj] >= counts[ci] * (inv_ratio ** m):
                # Strongest qualifying neighbour is the parent (deterministic tie-break).
                if counts[nj] > best_count or (counts[nj] == best_count
                                               and best is not None and junctions[nj] < junctions[best]):
                    best, best_count = nj, counts[nj]
        parent[ci] = best
    return parent


def _root(i: int, parent: list[int | None]) -> int:
    while parent[i] is not None:
        i = parent[i]  # type: ignore[assignment]
    return i


def correct_airr(
    airr_tsv: str | Path,
    output: str | Path,
    *,
    max_mismatches: int = 2,
    ratio: float = 0.05,
    require_vj: bool = False,
    complete_only: bool = True,
    read_map: str | Path | None = None,
    report_path: str | Path | None = None,
) -> CorrectReport:
    """Aggregate mapped reads into clonotypes and collapse CDR3 sequencing errors.

    Args:
        airr_tsv: Stage-1 mapped-reads AIRR TSV (needs ``junction``, ``sequence_id``).
        output: corrected clonotype table TSV (``junction``, ``junction_aa``, ``v_call``,
            ``j_call``, ``c_call``, ``locus``, ``duplicate_count``), sorted by abundance.
            A clonotype is keyed by ``(locus, v_call, j_call, junction)``; ``duplicate_count``
            (AIRR) is the number of distinct FRAGMENTS -- paired mates of one molecule counted
            once. ``c_call`` is the clonotype's dominant isotype CLASS (from ``c_class``: IGHG,
            IGHA, ...), preferring a resolved class over the ambiguous ``IGHC``; empty when no
            read carried a constant call.
        ratio: parent:child count ratio; must be in ``(0, 1)`` (vdjtools default 0.05).
        require_vj: only collapse neighbours sharing ``v_call`` and ``j_call``.
        complete_only: keep only reads whose junction spans both conserved anchors, is in
            frame, and has no stop codon (see :data:`_COMPLETE`). A read that stops short of
            the [FW]118 anchor yields a *prefix* of a junction, not a clonotype. Setting this
            ``False`` reproduces the raw per-read behaviour and is almost never what you want.
        read_map: optional TSV ``sequence_id -> junction`` (the corrected clonotype a
            read ends up in) — the read-id → junction map after correction.

    Returns:
        A :class:`CorrectReport`.
    """
    if not 0.0 < ratio < 1.0:
        # ratio == 0 -> 1/ratio divides by zero; ratio >= 1 -> a "parent" needs no more counts
        # than its child, which breaks the strictly-increasing-count invariant `_root` relies on
        # (mutual parents => infinite loop). vdjtools' default is 0.05.
        raise ValueError(f"ratio must be in (0, 1), got {ratio}")

    output = Path(output)
    df = pl.read_csv(airr_tsv, separator="\t", infer_schema_length=0)
    df = df.filter(pl.col("junction").is_not_null() & (pl.col("junction") != ""))
    n_with_junction = df.height
    if complete_only:
        df = df.filter(_COMPLETE)
    n_incomplete = n_with_junction - df.height

    # A clonotype is (locus, v_call, j_call, junction) -- NOT the junction alone. Two reads with the
    # same nucleotide junction but a different locus/V/J are different clonotypes; grouping on the
    # junction alone merged them and kept an arbitrary member's calls. Abundance is counted in
    # FRAGMENTS: paired mates `<id>/1` and `<id>/2` of one molecule both carry the junction, so
    # counting rows double-counts a fragment whenever both mates span it (insert-size dependent, so
    # the inflation is non-uniform across clonotypes -- exactly the count vector the ratio test
    # consumes). `read_ids` keeps every read so the read-map stays read-level.
    df = df.with_columns(pl.col("sequence_id").str.replace(r"/[12]$", "").alias("_frag"))
    keys = ["locus", "v_call", "j_call", "junction"]
    g = df.group_by(keys).agg(
        pl.col("_frag").n_unique().alias("count"),           # fragments, not reads
        pl.col("sequence_id").alias("read_ids"),
        pl.col("junction_aa").first().alias("junction_aa"),
    ).sort("count", descending=True)

    junctions = g["junction"].to_list()
    counts = [int(c) for c in g["count"].to_list()]
    v = [x or "" for x in g["v_call"].to_list()]
    j = [x or "" for x in g["j_call"].to_list()]
    read_ids = g["read_ids"].to_list()
    junction_aa = g["junction_aa"].to_list()
    locus = [x or "" for x in g["locus"].to_list()]
    # per-read isotype -> collapsed to the clonotype's dominant CLASS in the output (keyed by the full
    # sequence_id, mate suffix kept). Use ``c_class`` (the resolved class, e.g. IGHG, not the noisy
    # subclass IGHG1 -- IGHG1-4 are ~95% identical, so a per-read dominant subclass is a coin-flip);
    # fall back to ``c_call`` only if ``c_class`` is absent from the AIRR file.
    _iso_col = "c_class" if "c_class" in df.columns else ("c_call" if "c_call" in df.columns else None)
    cc_map = dict(zip(df["sequence_id"].to_list(), df[_iso_col].to_list())) if _iso_col else {}

    report = CorrectReport(clonotypes_in=len(junctions),
                           reads=sum(len(r) for r in read_ids),   # reads; `count` is fragments
                           reads_with_junction=n_with_junction,
                           reads_incomplete=n_incomplete)
    parent = _parents(junctions, counts, v, j, max_mismatches=max_mismatches,
                      ratio=ratio, require_vj=require_vj)

    # Accumulate each clonotype's count + reads into its ultimate ancestor.
    agg_count = counts[:]
    agg_reads: list[list[str]] = [list(r) for r in read_ids]
    order = sorted(range(len(junctions)), key=lambda i: counts[i])  # children first
    for i in order:
        p = parent[i]
        if p is None:
            continue
        r = _root(p, parent)
        agg_count[r] += agg_count[i]
        agg_reads[r].extend(agg_reads[i])
        agg_count[i] = 0
        agg_reads[i] = []

    roots = [i for i in range(len(junctions)) if parent[i] is None]
    roots.sort(key=lambda i: agg_count[i], reverse=True)
    report.clonotypes_out = len(roots)
    report.collapsed = report.clonotypes_in - report.clonotypes_out

    def _dominant_ccall(read_list: list[str]) -> str:
        # Prefer a RESOLVED isotype class: `isotype_class` emits the generic `IGHC` only when a read's
        # constant hit is ambiguous across classes (e.g. IGHG1,IGHM). Report IGHC only if NO read
        # resolves -- otherwise a handful of ambiguous reads could outvote the true class.
        calls = [c for sid in read_list if (c := cc_map.get(sid))]
        if not calls:
            return ""
        resolved = [c for c in calls if c not in _GENERIC_ISOTYPE]
        return Counter(resolved or calls).most_common(1)[0][0]

    out = pl.DataFrame({
        "junction": [junctions[i] for i in roots],
        "junction_aa": [junction_aa[i] for i in roots],
        "v_call": [v[i] for i in roots],
        "j_call": [j[i] for i in roots],
        "c_call": [_dominant_ccall(agg_reads[i]) for i in roots],
        "locus": [locus[i] for i in roots],
        "duplicate_count": [agg_count[i] for i in roots],
    })
    out.write_csv(output, separator="\t")

    if read_map is not None:
        rows = [(rid, junctions[i]) for i in roots for rid in agg_reads[i]]
        pl.DataFrame(rows, schema=["sequence_id", "junction"], orient="row").write_csv(
            read_map, separator="\t")
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report
