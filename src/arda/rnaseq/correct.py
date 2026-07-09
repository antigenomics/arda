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
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ..refbuild.translate import reverse_complement

__all__ = ["correct_airr", "CorrectReport"]

_DNA = frozenset("ACGT")

# The generic heavy constant `isotype_class` returns when a read's C hit spans classes (IGHG1,IGHM
# -> IGHC): isotype unresolved. Deprioritised when aggregating a clonotype's dominant isotype.
_GENERIC_ISOTYPE = frozenset({"IGHC"})


def _strip_mate(sid: str) -> str:
    """``<id>/1`` / ``<id>/2`` -> ``<id>`` (the fragment id shared by paired mates)."""
    return sid[:-2] if sid[-2:] in ("/1", "/2") else sid

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


def _assign_coverage(
    raw: pl.DataFrame,
    root_jn: list[str],
    root_loc: list[str],
    exact: dict[tuple[str, str, str, str], int],
    *,
    k: int = 12,
    min_ov: int = 20,
    max_mm: float = 0.12,
    cap: int = 64,
) -> list[list[str]]:
    """Assign every CDR3-overlapping read to the clonotype whose junction it belongs to.

    A clonotype's expression is *all* reads that encompass its junction, not only the reads that
    span it end-to-end: a long CDR3 is covered by many partial reads (V-side, J-side) that never
    reach both anchors. This is the read-counting the assembly tools (MiXCR/TRUST4) do, and it is
    what makes cross-tool frequencies comparable. Returns, parallel to ``root_jn``, the list of
    ``sequence_id`` s assigned to each clonotype. Each read is counted once.

    Pass 1 -- exact: a read whose ``(locus, v_call, j_call, junction)`` key is in ``exact`` (the
    clonotype key -> root-position map, collapsed children included) is assigned there
    (authoritative; also the path for a read with no ``sequence`` column). The key must match on
    V/J too, not the junction alone -- two clonotypes can share a junction nt under different V/J
    alleles (common in IGK), and keying on the junction alone routes both to one and leaves the
    other with zero reads.

    Pass 2 -- align: a remaining read's sequence is aligned to the root junctions (a shared k-mer
    fixes the offset; it joins the root with the longest ``>= min_ov`` overlap within the ``max_mm``
    per-base mismatch budget, tolerating SHM). Reads overlapping no junction by ``min_ov`` (pure
    germline V/C) stay unassigned -- ambiguous across every clonotype of that V.
    """
    assigned: list[list[str]] = [[] for _ in root_jn]
    n = raw.height
    cols = {c: raw[c].to_list() for c in raw.columns}

    def col(name):
        return cols.get(name, [None] * n)

    seqc, rcc, jnc, sidc = col("sequence"), col("rev_comp"), col("junction"), col("sequence_id")
    locc, vc, jc, cc = col("locus"), col("v_call"), col("j_call"), col("c_call")
    done: set[str] = set()

    # Pass 1: exact clonotype-key match (spanning + rescued reads; authoritative).
    for i in range(n):
        sid = sidc[i]
        if sid in done:                              # a rescued read appears in mapped + assembled rows
            continue
        rp = exact.get(((locc[i] or ""), (vc[i] or ""), (jc[i] or ""), (jnc[i] or "")))
        if rp is not None:
            assigned[rp].append(sid); done.add(sid)

    # Pass 2: align the rest (partial V-side / J-side reads that never reached a complete junction).
    index: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for ri, jn in enumerate(root_jn):
        for p in range(len(jn) - k + 1):
            lst = index[jn[p:p + k]]
            if len(lst) < cap:                       # bound the germline-shared k-mers
                lst.append((ri, p))
    for i in range(n):
        sid = sidc[i]
        if sid in done:
            continue
        seq = seqc[i]
        if not seq or len(seq) < k:
            continue
        loc = (locc[i] or "")[:3] or _gene3(vc[i]) or _gene3(jc[i]) or _gene3(cc[i])
        s = reverse_complement(seq) if str(rcc[i]).upper() in ("T", "TRUE", "1") else seq
        best_ri, best_ov, seen = None, min_ov - 1, set()
        L = len(s)
        for rp in range(0, L - k + 1):
            for (ri, jp) in index.get(s[rp:rp + k], ()):
                if loc and root_loc[ri] != loc:
                    continue
                d = jp - rp
                if (ri, d) in seen:
                    continue
                seen.add((ri, d))
                jr = root_jn[ri]
                lo, hi = (-d if d < 0 else 0), min(L, len(jr) - d)
                ov = hi - lo
                if ov <= best_ov:
                    continue
                budget = max_mm * ov
                mm = 0
                for kk in range(lo, hi):
                    if s[kk] != jr[kk + d]:
                        mm += 1
                        if mm > budget:
                            break
                if mm <= budget:
                    best_ov, best_ri = ov, ri
        if best_ri is not None:
            assigned[best_ri].append(sid); done.add(sid)
    return assigned


def _gene3(x: str | None) -> str:
    g = (x or "").split(",")[0].split("(")[0].strip()
    return g[:3] if g[:3] in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD") else ""


def correct_airr(
    airr_tsv: str | Path,
    output: str | Path,
    *,
    max_mismatches: int = 2,
    ratio: float = 0.05,
    require_vj: bool = False,
    complete_only: bool = True,
    coverage: bool = True,
    read_map: str | Path | None = None,
    extra_airr: str | Path | None = None,
    report_path: str | Path | None = None,
) -> CorrectReport:
    """Aggregate mapped reads into clonotypes and collapse CDR3 sequencing errors.

    Args:
        airr_tsv: Stage-1 mapped-reads AIRR TSV (needs ``junction``, ``sequence_id``).
        output: corrected clonotype table TSV (``junction``, ``junction_aa``, ``v_call``,
            ``j_call``, ``c_call``, ``locus``, ``duplicate_count``, ``consensus_count``), sorted
            by abundance. A clonotype is keyed by ``(locus, v_call, j_call, junction)``. Per the
            AIRR schema, ``duplicate_count`` is the number of READS supporting the clonotype (both
            paired mates of a molecule count) -- directly comparable to MiXCR ``readCount`` /
            TRUST4 ``#count`` -- and ``consensus_count`` is the number of distinct fragment
            consensuses (the two mates of one molecule are one consensus). ``c_call`` is the
            clonotype's dominant isotype CLASS (from ``c_class``: IGHG, IGHA, ...), preferring a
            resolved class over the ambiguous ``IGHC``; empty when no read carried a constant call.
        ratio: parent:child count ratio; must be in ``(0, 1)`` (vdjtools default 0.05).
        require_vj: only collapse neighbours sharing ``v_call`` and ``j_call``.
        complete_only: keep only reads whose junction spans both conserved anchors, is in
            frame, and has no stop codon (see :data:`_COMPLETE`). A read that stops short of
            the [FW]118 anchor yields a *prefix* of a junction, not a clonotype. Setting this
            ``False`` reproduces the raw per-read behaviour and is almost never what you want.
            (This governs which reads DEFINE clonotypes, not how they are counted -- see ``coverage``.)
        coverage: count a clonotype's abundance as EVERY read that encompasses its junction
            (aligns to it), not only the reads that span it end-to-end (default ``True``). A long
            CDR3 is covered by many partial V-side / J-side reads that never reach both anchors;
            counting only spanning reads under-reports it non-uniformly (the deficit scales with
            CDR3 length), which is why a spanning count correlates poorly with the assembly tools'
            read counts. Coverage counting (:func:`_assign_coverage`) matches MiXCR/TRUST4 and is
            the true expression estimate. ``False`` reverts to spanning-read counts.
        read_map: optional TSV ``sequence_id -> junction`` (the corrected clonotype a
            read ends up in) — the read-id → junction map after correction.
        extra_airr: optional Stage-3 assembled-reads AIRR (from
            :func:`~arda.rnaseq.assemble.assemble_contigs`), concatenated with ``airr_tsv``
            before aggregation. Its rows carry a contig's complete junction for reads whose own
            Stage-1 junction was incomplete, so a long-CDR3 clone no single read spans is counted
            once (the read's incomplete Stage-1 row is dropped by ``complete_only``).

    Returns:
        A :class:`CorrectReport`.
    """
    if not 0.0 < ratio < 1.0:
        # ratio == 0 -> 1/ratio divides by zero; ratio >= 1 -> a "parent" needs no more counts
        # than its child, which breaks the strictly-increasing-count invariant `_root` relies on
        # (mutual parents => infinite loop). vdjtools' default is 0.05.
        raise ValueError(f"ratio must be in (0, 1), got {ratio}")

    output = Path(output)
    raw = pl.read_csv(airr_tsv, separator="\t", infer_schema_length=0)
    if extra_airr is not None:
        extra = pl.read_csv(extra_airr, separator="\t", infer_schema_length=0)
        if extra.height:
            raw = pl.concat([raw, extra], how="diagonal")
    # Isotype lives on the CONSTANT-region reads: they carry ``c_class`` but no junction, so the
    # complete-only filter below drops them, and the JUNCTION reads that build the clonotype carry no
    # ``c_class`` of their own. Link the two by FRAGMENT id (paired mates share ``<id>``): map each
    # fragment -> the isotype class(es) any of its reads carried, before any filtering.
    frag_iso: dict[str, list[str]] = {}
    if "c_class" in raw.columns:
        for sid, cl in zip(raw["sequence_id"].to_list(), raw["c_class"].to_list()):
            if cl:
                frag_iso.setdefault(_strip_mate(sid), []).append(cl)
    df = raw.filter(pl.col("junction").is_not_null() & (pl.col("junction") != ""))
    n_with_junction = df.height
    if complete_only:
        df = df.filter(_COMPLETE)
    n_incomplete = n_with_junction - df.height

    # A clonotype is (locus, v_call, j_call, junction) -- NOT the junction alone. Two reads with the
    # same nucleotide junction but a different locus/V/J are different clonotypes. Per the AIRR schema
    # the table reports BOTH counts (see the output DataFrame): `duplicate_count` = READS (every row,
    # so a molecule whose two mates both span the junction counts twice) -- the same read-counting
    # convention as MiXCR `readCount` / TRUST4 `#count`, which is what makes the abundances directly
    # comparable across tools -- and `consensus_count` = distinct FRAGMENTS (`_frag`, the two mates of
    # one molecule collapsed to one consensus). The error-correction ratio test runs on the fragment
    # (consensus) count, which is insert-size-invariant (reads are inflated non-uniformly by insert
    # size). `read_ids` keeps every read so the read-map stays read-level.
    df = df.with_columns(pl.col("sequence_id").str.replace(r"/[12]$", "").alias("_frag"))
    keys = ["locus", "v_call", "j_call", "junction"]
    g = df.group_by(keys).agg(
        pl.col("_frag").n_unique().alias("count"),           # fragments (consensuses), not reads
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
    report.clonotypes_out = len(roots)
    report.collapsed = report.clonotypes_in - report.clonotypes_out

    # Per-clonotype read set. Coverage (default): every read that ENCOMPASSES the junction, assigned
    # by alignment to the final (post-collapse) root junctions -- the true expression. Spanning:
    # the reads that reached a complete junction, accumulated along parent pointers.
    if coverage:
        pos = {ri: p for p, ri in enumerate(roots)}                    # global root index -> roots position
        exact = {}                                                     # clonotype key (child too) -> root position
        for i in range(len(junctions)):
            exact.setdefault((locus[i], v[i], j[i], junctions[i]), pos[_root(i, parent)])
        read_sets = _assign_coverage(raw, [junctions[i] for i in roots], [locus[i] for i in roots], exact)
    else:
        read_sets = [agg_reads[i] for i in roots]
    dup = [len(rs) for rs in read_sets]                                 # AIRR duplicate_count: reads
    cons = [len({_strip_mate(x) for x in rs}) for rs in read_sets]      # AIRR consensus_count: fragments

    def _dominant_ccall(read_list: list[str]) -> str:
        # A clonotype's isotype = the dominant RESOLVED class over its fragments' constant mates.
        # `isotype_class` emits the generic `IGHC` only on cross-class ambiguity (IGHG1,IGHM), so
        # report IGHC only if NO read resolves -- a handful of ambiguous reads must not outvote it.
        calls: list[str] = []
        for sid in read_list:
            calls.extend(frag_iso.get(_strip_mate(sid), ()))
        if not calls:
            return ""
        resolved = [c for c in calls if c not in _GENERIC_ISOTYPE]
        return Counter(resolved or calls).most_common(1)[0][0]

    order = sorted(range(len(roots)), key=lambda r: (dup[r], cons[r]), reverse=True)
    out = pl.DataFrame({
        "junction": [junctions[roots[r]] for r in order],
        "junction_aa": [junction_aa[roots[r]] for r in order],
        "v_call": [v[roots[r]] for r in order],
        "j_call": [j[roots[r]] for r in order],
        "c_call": [_dominant_ccall(read_sets[r]) for r in order],
        "locus": [locus[roots[r]] for r in order],
        "duplicate_count": [dup[r] for r in order],              # reads encompassing the junction
        "consensus_count": [cons[r] for r in order],             # distinct fragment consensuses
    })
    out.write_csv(output, separator="\t")

    if read_map is not None:
        rows = [(rid, junctions[roots[r]]) for r in range(len(roots)) for rid in read_sets[r]]
        pl.DataFrame(rows, schema=["sequence_id", "junction"], orient="row").write_csv(
            read_map, separator="\t")
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report
