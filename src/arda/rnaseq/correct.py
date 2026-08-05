"""Stage 2 — CDR3 error correction (sequencing-error model).

Collapses sequencing-error CDR3 variants onto their parent clonotype, using :mod:`seqtree`
neighbour search (a fast edit-bounded index) to find substitution/indel neighbours. A clonotype
``C`` is an error **child** of a more-abundant neighbour ``P`` (differing by ``n_subs``
substitutions and ``n_indel`` inserted/deleted bases) iff the expected number of such misread
parent reads -- ``count[P] * p_sub**n_subs * p_ind**n_indel`` -- is at least ``count[C]``. The rates
are PER BASE and the per-mismatch probability is length-scaled (``p_sub = error_rate * L``): a single
mismatch over a longer junction sheds proportionally more error mass, so the default
``error_rate = 0.001`` reproduces vdjtools' ~1/20 at a 45 nt (15 aa) junction and scales elsewhere. A
multi-base (in-frame SHM) indel costs ``p_ind**len`` and is kept as a real clonotype. The count is
the SPANNING read depth -- reads that fully observe the
junction -- so the test is over the reads that actually saw the discriminating base (``"2/2, not
2/200"``); ``error_method`` in {binom, betabinom} instead piles up partial reads per position for
extra depth at very low coverage. Children route to the parent; chains collapse to the ultimate
ancestor; ``count[parent] * p_err >= count[child]`` with ``p_err < 1`` gives strictly increasing
counts along parent pointers, so there are no cycles.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ._res import Stage
import seqtree

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
    # See `_res.Stage`: peak is the WHOLE-PROCESS high-water mark as of this stage's end
    # (monotone -- getrusage offers no per-stage reset), gain is this stage's contribution.
    wall_seconds: float = 0.0
    peak_rss_mb: float = 0.0
    rss_gain_mb: float = 0.0

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["collapse_fraction"] = self.collapsed / self.clonotypes_in if self.clonotypes_in else 0.0
        return d


def _parents(junctions: list[str], counts: list[int], v: list[str], j: list[str],
             *, max_subs: int, max_indel: int, error_rate: float, indel_rate: float,
             require_vj: bool) -> list[int | None]:
    """For each clonotype, the strongest neighbour that is its sequencing-error parent (or None).

    ``counts`` is the SPANNING read count -- reads that fully observe the junction -- so the test is
    made only over reads that actually saw the discriminating position. A clonotype seen 2 times out
    of the 2 reads that reached its position is 100 % there (real), not a 1 % error of an abundant
    neighbour whose reads never covered that position ("2/2, not 2/200").

    A neighbour ``C`` is an error child of ``P`` (differing by ``n_subs`` substitutions and
    ``n_indel`` inserted/deleted bases) iff the expected number of such error reads,
    ``count[P] * p_sub**n_subs * p_ind**n_indel``, is at least ``count[C]``. ``error_rate`` and
    ``indel_rate`` are PER-BASE rates, and the per-substitution collapse probability scales with the
    junction length -- ``p_sub = error_rate * L`` -- because a single mismatch over a longer junction
    sheds proportionally more error mass (this is why vdjtools' 1/20 is calibrated for a 45 nt / 15 aa
    junction; ``error_rate = 0.001`` reproduces it there and scales for other lengths). Because a
    multi-base indel costs ``p_ind**len``, a 3-9 bp (in-frame SHM) indel is vanishingly unlikely as an
    error and kept as a real clonotype, while a 1 bp (instrument) indel collapses. Children route to
    the parent; chains collapse to the ancestor. ``count[parent] * p_err >= count[child]`` with
    ``p_err < 1`` makes counts strictly increase along parent pointers -> no cycles.
    """
    n = len(junctions)
    parent: list[int | None] = [None] * n
    safe = [i for i, s in enumerate(junctions) if s and set(s) <= _DNA]
    if not safe:
        return parent
    index = seqtree.Index.build([junctions[i] for i in safe], alphabet="nt")
    params = seqtree.SearchParams(max_subs=max_subs, max_ins=max_indel, max_dels=max_indel,
                                  max_total_edits=max_subs + max_indel, engine="seqtm")
    for ci in safe:
        L = len(junctions[ci])
        p_sub = min(0.5, error_rate * L)          # per-substitution collapse prob, scaled by length
        p_ind = min(0.5, indel_rate * L)
        best, best_count = None, -1
        for hit in index.search(junctions[ci], params):
            nj = safe[hit.ref_id]
            if nj == ci or (hit.n_subs + hit.n_ins + hit.n_dels) == 0:
                continue
            if require_vj and (v[nj] != v[ci] or j[nj] != j[ci]):
                continue
            p_err = p_sub ** hit.n_subs * p_ind ** (hit.n_ins + hit.n_dels)
            if counts[nj] * p_err >= counts[ci]:
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


def _binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), summed as 1 - CDF(k-1) (k terms; k is a small error count)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    from math import comb
    cdf = sum(comb(n, i) * p ** i * (1.0 - p) ** (n - i) for i in range(k))
    return max(0.0, 1.0 - cdf)


def _betabinom_sf(k: int, n: int, p: float, rho: float = 0.1) -> float:
    """P(X >= k) for a Beta-Binomial with mean ``p`` and overdispersion ``rho`` -- fatter tail than
    the binomial, so at very low depth a couple of correlated miscalls are not over-called as real."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    from math import lgamma
    s = (1.0 - rho) / rho if 0 < rho < 1 else 1e9
    a, b = p * s, (1.0 - p) * s
    lb = lgamma(a + b) - lgamma(a) - lgamma(b)
    cdf = 0.0
    for i in range(k):
        lcomb = lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1)
        cdf += (2.718281828459045 ** (lcomb + lgamma(i + a) + lgamma(n - i + b)
                                      - lgamma(n + a + b) + lb))
    return max(0.0, 1.0 - cdf)


def _error_pileup(
    raw: pl.DataFrame,
    junctions: list[str],
    v: list[str],
    j: list[str],
    locus: list[str],
    parent_simple: list[int | None],
    span_counts: list[int],
    *,
    error_rate: float,
    indel_rate: float,
    method: str,
    alpha: float = 1e-3,
    k: int = 12,
    cap: int = 64,
) -> list[int | None]:
    """Re-decide parentage from per-position read DEPTH -- the deep, low-coverage path.

    The simple test only counts reads that span the whole junction; at very low coverage there may be
    only one or two. Here every read (partial included) is aligned to the raw clonotype junctions and
    piled up per position, so a substitution difference is judged on the reads that actually covered
    THAT base. For candidate parent P of child C differing at substitution positions D, C is an error
    child iff at every ``d in D`` the child-allele depth is consistent with sequencing error of the
    parent-allele depth -- ``sf(child_depth; child_depth + parent_depth, error_rate) > alpha`` (a
    binomial tail, ``betabinom`` for an overdispersed one). Indel candidates (gapped, no fixed
    positions) keep the simple decision.
    """
    sf = _betabinom_sf if method == "betabinom" else _binom_sf
    index: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for ri, jn in enumerate(junctions):
        for p in range(len(jn) - k + 1):
            lst = index[jn[p:p + k]]
            if len(lst) < cap:
                lst.append((ri, p))
    depth = [[0] * len(jn) for jn in junctions]                       # per clonotype, per-position read depth
    cols = {c: raw[c].to_list() for c in raw.columns}
    n = raw.height

    def col(name):
        return cols.get(name, [None] * n)

    seqc, rcc, locc, vc, jc, cc = (col("sequence"), col("rev_comp"), col("locus"),
                                   col("v_call"), col("j_call"), col("c_call"))
    for i in range(n):
        seq = seqc[i]
        if not seq or len(seq) < k:
            continue
        loc = (locc[i] or "")[:3] or _gene3(vc[i]) or _gene3(jc[i]) or _gene3(cc[i])
        s = reverse_complement(seq) if str(rcc[i]).upper() in ("T", "TRUE", "1") else seq
        best_ri, best_ov, best_lo, best_hi, seen = None, 19, 0, 0, set()
        L = len(s)
        for rp in range(0, L - k + 1):
            for (ri, jp) in index.get(s[rp:rp + k], ()):
                if loc and locus[ri][:3] != loc:
                    continue
                d = jp - rp
                if (ri, d) in seen:
                    continue
                seen.add((ri, d))
                jr = junctions[ri]
                lo, hi = (-d if d < 0 else 0), min(L, len(jr) - d)
                ov = hi - lo
                if ov <= best_ov:
                    continue
                if sum(1 for kk in range(lo, hi) if s[kk] != jr[kk + d]) <= 0.12 * ov:
                    best_ri, best_ov, best_lo, best_hi = ri, ov, lo + d, hi + d
        if best_ri is not None:
            dep = depth[best_ri]
            for p in range(best_lo, best_hi):
                dep[p] += 1

    parent = list(parent_simple)
    safe = [i for i, s in enumerate(junctions) if s and set(s) <= _DNA]
    idx = seqtree.Index.build([junctions[i] for i in safe], alphabet="nt")
    params = seqtree.SearchParams(max_subs=2, max_ins=0, max_dels=0, max_total_edits=2, engine="seqtm")
    for ci in safe:
        jc_str = junctions[ci]
        best, best_count = None, -1
        for hit in idx.search(jc_str, params):
            nj = safe[hit.ref_id]
            if nj == ci or hit.n_subs == 0 or hit.n_ins or hit.n_dels:
                continue
            if v[nj] != v[ci] or j[nj] != j[ci] or len(junctions[nj]) != len(jc_str):
                continue
            disc = [p for p in range(len(jc_str)) if jc_str[p] != junctions[nj][p]]
            if not disc:
                continue
            ok = True
            for d in disc:
                cd, pd = depth[ci][d], depth[nj][d]
                if cd == 0:
                    ok = False
                    break
                if sf(cd, cd + pd, error_rate) <= alpha:    # child allele too deep to be error
                    ok = False
                    break
            if ok and span_counts[nj] > best_count:
                best, best_count = nj, span_counts[nj]
        if best is not None:
            parent[ci] = best
    return parent


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
    reach both anchors. This is the coverage read-counting an assembly-based extractor does, and it
    is the true expression estimate. Returns, parallel to ``root_jn``, the list of ``sequence_id`` s
    assigned to each clonotype. Each read is counted once.

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
            assigned[rp].append(sid)
            done.add(sid)

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
            assigned[best_ri].append(sid)
            done.add(sid)
    return assigned


def _gene3(x: str | None) -> str:
    g = (x or "").split(",")[0].split("(")[0].strip()
    return g[:3] if g[:3] in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD") else ""


def _clonotype_d(out: pl.DataFrame, organism: str) -> list[pl.Series]:
    """D (and tandem D-D) per clonotype, mapped into its error-corrected junction.

    Reads carry their own ``d_call``, but a read's D is called on a sequencing-error copy of
    the junction and, for a long CDR3, on a read that does not span the D at all. The clonotype
    is the first place the junction is both complete and corrected, and D is a deterministic
    function of ``(junction, v_call, j_call)`` -- so call it once here rather than voting over
    reads. Costs one gapless alignment per clonotype, not per read.

    VJ loci and organisms without D germlines come back empty, as does an unresolvable V/J.
    """
    from ..annotate.dmap import map_d_junction

    cols = {c: [] for c in ("d_call", "d2_call", "d_support", "d2_support")}
    for jn, vc, jc in zip(out["junction"], out["v_call"], out["j_call"]):
        try:
            call = map_d_junction(jn or "", vc or "", jc or "", organism)
        except (KeyError, ValueError):        # unknown allele / organism without anchors
            call = None
        cols["d_call"].append(call.d_call if call else "")
        cols["d2_call"].append(call.d2_call if call else "")
        cols["d_support"].append(call.d_support if call else "")
        cols["d2_support"].append(call.d2_support if call else "")
    return [pl.Series(k, v, dtype=pl.Utf8) for k, v in cols.items()]


def correct_airr(
    airr_tsv: str | Path,
    output: str | Path,
    *,
    organism: str = "human",
    map_d: bool = True,
    max_subs: int = 2,
    max_indel: int = 0,
    error_rate: float = 0.001,
    indel_rate: float = 0.001,
    require_vj: bool = True,
    error_method: str = "simple",
    complete_only: bool = True,
    coverage: bool = True,
    read_map: str | Path | None = None,
    extra_airr: str | Path | None = None,
    report_path: str | Path | None = None,
) -> CorrectReport:
    """Aggregate mapped reads into clonotypes and collapse CDR3 sequencing errors.

    Args:
        airr_tsv: Stage-1 mapped-reads AIRR TSV (needs ``junction``, ``sequence_id``).
        organism: reference organism, used only to map D into each clonotype's junction.
        map_d: append ``d_call``/``d2_call``/``d_support``/``d2_support``, called once per
            clonotype on its corrected junction (see :func:`_clonotype_d`). Default ``True``.
        output: corrected clonotype table TSV (``junction``, ``junction_aa``, ``v_call``,
            ``j_call``, ``c_call``, ``locus``, ``duplicate_count``, ``consensus_count``, and
            with ``map_d`` the four D columns), sorted
            by abundance. A clonotype is keyed by ``(locus, v_call, j_call, junction)``. Per the
            AIRR schema, ``duplicate_count`` is the number of READS supporting the clonotype (both
            paired mates of a molecule count) and ``consensus_count`` is the number of distinct
            fragment consensuses (the two mates of one molecule are one consensus). ``c_call`` is the
            clonotype's dominant isotype CLASS (from ``c_class``: IGHG, IGHA, ...), preferring a
            resolved class over the ambiguous ``IGHC``; empty when no read carried a constant call.
        max_subs: max substitutions between an error child and its parent (seqtree neighbour search).
        max_indel: max inserted/deleted bases searched for indel error children (default 0). A 1-2 bp
            instrument indel is a frameshift and is already dropped by ``complete_only``, so on
            complete junctions the indel search only costs time (~160x slower) and collapses nothing;
            a multi-base in-frame SHM indel costs ``(indel_rate*L)**len`` and is kept as a real
            clonotype either way. Set it > 0 only with ``--all-junctions`` (frameshift indels kept).
        error_rate: per-BASE substitution error rate (~Phred 30 = 0.001). The per-substitution
            collapse probability is length-scaled, ``error_rate * junction_len``, so the default
            reproduces vdjtools' ~1/20 at a 45 nt (15 aa) junction and scales for other lengths.
        indel_rate: per-BASE indel error rate (instrument-dependent; default 0.001, length-scaled).
        error_method: ``"simple"`` (default) tests on spanning read counts; ``"binom"`` /
            ``"betabinom"`` pile up partial reads per discriminating position for extra depth at
            very low coverage (:func:`_error_pileup`).
        require_vj: only collapse neighbours sharing ``v_call`` and ``j_call`` (default ``True`` -- a
            true sequencing error does not change the germline-anchored V/J call).
        complete_only: keep only reads whose junction spans both conserved anchors, is in
            frame, and has no stop codon (see :data:`_COMPLETE`). A read that stops short of
            the [FW]118 anchor yields a *prefix* of a junction, not a clonotype. Setting this
            ``False`` reproduces the raw per-read behaviour and is almost never what you want.
            (This governs which reads DEFINE clonotypes, not how they are counted -- see ``coverage``.)
        coverage: count a clonotype's abundance as EVERY read that encompasses its junction
            (aligns to it), not only the reads that span it end-to-end (default ``True``). A long
            CDR3 is covered by many partial V-side / J-side reads that never reach both anchors;
            counting only spanning reads under-reports it non-uniformly (the deficit scales with
            CDR3 length). Coverage counting (:func:`_assign_coverage`) is the true expression
            estimate. ``False`` reverts to spanning-read counts.
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
    for name, val in (("error_rate", error_rate), ("indel_rate", indel_rate)):
        # p_err < 1 keeps counts strictly increasing along parent pointers (no cycles); p_err == 0
        # would make a single mismatch collapse anything, p_err >= 1 never collapses.
        if not 0.0 < val < 1.0:
            raise ValueError(f"{name} must be in (0, 1), got {val}")
    if error_method not in ("simple", "binom", "betabinom"):
        raise ValueError(f"error_method must be simple|binom|betabinom, got {error_method!r}")

    stage = Stage()
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
    # same nucleotide junction but a different locus/V/J are different clonotypes. `read_ids` keeps
    # every SPANNING read (a complete-junction read fully observes the junction); its length is the
    # spanning read count the error-correction test runs on, and it also keeps the read-map read-level.
    # `count` is distinct fragments (`_frag`, the two mates of one molecule collapsed to one consensus),
    # reported as `consensus_count` when abundance is spanning (`coverage=False`).
    df = df.with_columns(pl.col("sequence_id").str.replace(r"/[12]$", "").alias("_frag"))
    keys = ["locus", "v_call", "j_call", "junction"]
    # Every ordering below is load-bearing; without them `correct` is not reproducible AT ALL.
    # Measured on 200k reads: three runs, same input, same flags -> three different clones.tsv.
    # polars' group_by is a multithreaded hash aggregation, so all three of these were arbitrary:
    #   * group order          -> `_parents` collapses error children onto whichever parent it met
    #                             first, so an identical junction flipped between paralogous V calls
    #                             (IGHV3-11*06 vs IGHV3-21*08) from run to run;
    #   * read order in a group -> `_assign_coverage` is first-with-longest-overlap-wins, so
    #                             `duplicate_count` moved (11/8 vs 9/6 for one IGK clonotype);
    #   * equal `count` rows    -> emitted in arbitrary order.
    # Sorting the input makes `.first()` well defined; sorting `read_ids` fixes coverage; and the
    # final sort carries the full group key, which is a TOTAL order -- so the row order no longer
    # depends on the input row order either. That last property is what lets a sharded or Nextflow
    # run be byte-identical to a single-node one.
    g = df.sort("sequence_id").group_by(keys, maintain_order=True).agg(
        pl.col("_frag").n_unique().alias("count"),           # fragments (consensuses), not reads
        pl.col("sequence_id").sort().alias("read_ids"),
        pl.col("junction_aa").first().alias("junction_aa"),
    ).sort(["count", *keys], descending=[True, False, False, False, False])

    junctions = g["junction"].to_list()
    counts = [int(c) for c in g["count"].to_list()]
    v = [x or "" for x in g["v_call"].to_list()]
    j = [x or "" for x in g["j_call"].to_list()]
    read_ids = g["read_ids"].to_list()
    junction_aa = g["junction_aa"].to_list()
    locus = [x or "" for x in g["locus"].to_list()]
    span_counts = [len(r) for r in read_ids]                 # spanning reads: the error-test count

    report = CorrectReport(clonotypes_in=len(junctions),
                           reads=sum(span_counts),
                           reads_with_junction=n_with_junction,
                           reads_incomplete=n_incomplete)
    parent = _parents(junctions, span_counts, v, j, max_subs=max_subs, max_indel=max_indel,
                      error_rate=error_rate, indel_rate=indel_rate, require_vj=require_vj)
    if error_method != "simple":
        parent = _error_pileup(raw, junctions, v, j, locus, parent, span_counts,
                               error_rate=error_rate, indel_rate=indel_rate, method=error_method)

    # Accumulate each clonotype's count + reads into its ultimate ancestor.
    agg_count = counts[:]
    agg_reads: list[list[str]] = [list(r) for r in read_ids]
    order = sorted(range(len(junctions)), key=lambda i: span_counts[i])  # children first
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

    # Sort by abundance, then break every tie deterministically. Ranking on counts alone left
    # tied clonotypes in read order, and read order comes from a threaded mmseqs search -- so
    # the same input produced the same rows in a different sequence from run to run.
    order = sorted(range(len(roots)),
                   key=lambda r: (-dup[r], -cons[r], junctions[roots[r]],
                                  v[roots[r]], j[roots[r]]))
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
    if map_d:
        out = out.with_columns(_clonotype_d(out, organism))
    out.write_csv(output, separator="\t")

    if read_map is not None:
        rows = [(rid, junctions[roots[r]]) for r in range(len(roots)) for rid in read_sets[r]]
        pl.DataFrame(rows, schema=["sequence_id", "junction"], orient="row").write_csv(
            read_map, separator="\t")
    stage.finish(report)
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report
