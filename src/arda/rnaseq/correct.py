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

from ..annotate.airr_out import read_airr as _read_airr

from ._res import Stage
import seqtree

from ..refbuild.translate import reverse_complement
from .denoise import (REGIMES, DenoiseReport, clonotype_quality, quality_rescue,
                      read_quality)

__all__ = ["correct_airr", "CorrectReport", "CLONOTYPE_KEYS"]



#: How a clonotype is identified.
#:
#: ``full`` -- ``(locus, v_call, j_call, junction)``, the historical key.
#: ``junction`` -- ``(locus, junction)``: the V/J calls are canonicalised to the junction's
#: majority before grouping, so **call splits collapse**. That class is invisible to every error
#: model -- a junction byte-identical to an abundant clone's under a different V or J call has no
#: discriminating base to score -- and on Jurkat it is the largest error class BY READS (130 of
#: 14,531), including an allele-level TRG split. Measured cost on a POLYCLONAL TRA amplicon: 132 of
#: 19,956 clonotypes merge (0.66 %), and in every ambiguous case inspected the minority call
#: carried ONE read against 4-10 for the majority on a short 30-39 nt junction -- a call error on a
#: low-abundance read, not a second clone. Benefit on Jurkat: TRB 35 -> 33 clonotypes at purity
#: .99096 -> .99696, reads unchanged.
CLONOTYPE_KEYS = ("full", "junction")

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


#: ``--ec-mode`` presets. Each names an ``error_method`` and a ``min_junction_q``; an explicitly
#: passed knob always wins (:func:`correct_airr` resolves ``None`` against the mode).
#:
#: ⛔ ``accurate`` is NOT ``binom``/``betabinom``. Those pile up partial reads per discriminating
#: position for extra depth at very low coverage, which sounds like the accurate answer and is not
#: one here: measured on a 302,172-read MIGEC library with one 293 k-read clone, `simple` takes
#: 0.73 s and 143 clonotypes, `binom` 197 s / 79 and `betabinom` 254 s / 78 — ~270x slower AND
#: more aggressive, i.e. they collapse MORE real low-frequency variants, which is the failure this
#: mode exists to avoid. On a monoclonal Jurkat library all three are byte-identical (90
#: clonotypes, 0.35/0.41/0.39 s). So the depth models earn no place in a shipped mode; what
#: `accurate` buys is the QUALITY gate, which is evidence the simple model never had.
#:
#: Q20 is the LOW end of the measured plateau, not the optimum on any one library. The gate's
#: effect is flat over Q20-32 and degrades by Q35 (on the MIGEC spike-ins the abundant published
#: variant loses 1,094 -> 612 reads there), and on that library Q25-30 is marginally better on
#: every axis than Q20. Shipping the conservative end is deliberate: `accurate` exists to protect
#: rare real variants, and one library is not enough to tune a default with.
#: ``amplicon`` and ``rnaseq`` additionally switch on the QUALITY-DIRECTED RESCUE
#: (:mod:`arda.rnaseq.denoise`), which reaches the class the abundance model structurally cannot:
#: a clonotype 4+ substitutions from anything has no ladder of observed intermediates behind it
#: (0 of 13 at k=4 on Jurkat, against 0.0019 predicted) and no discriminating base for
#: ``--min-junction-q`` to judge, but its reads are measurably bad. They are ROUTED to an abundant
#: parent, never deleted; one with no parent is kept.
EC_MODES: dict[str, dict] = {
    "fast": {"error_method": "simple", "min_junction_q": 0, "regime": "fast"},
    "accurate": {"error_method": "simple", "min_junction_q": 20, "regime": "accurate"},
    "amplicon": {"error_method": "simple", "min_junction_q": 20, "regime": "amplicon"},
    "rnaseq": {"error_method": "simple", "min_junction_q": 20, "regime": "rnaseq"},
}


@dataclass
class CorrectReport:
    clonotypes_in: int = 0
    clonotypes_out: int = 0
    #: ⛔ SPANNING reads entering Stage 2, counted BEFORE any correction runs. It is therefore
    #: invariant to everything ``correct`` does, and it is **NOT** the read-conservation quantity.
    #: Comparing it across ``--ec-mode`` shows 0 on every sample and reads exactly like the
    #: invariant holding -- which is how a 1.39 % leak on a full-depth TRA amplicon was missed.
    #: The invariant is :attr:`reads_assigned`.
    reads: int = 0
    #: **The read-conservation invariant**: ``sum(duplicate_count)`` over the emitted clonotype
    #: table, i.e. every read the correction actually placed. This is what must not fall when a
    #: denoising mode is switched on -- error correction MOVES reads onto a parent, it never
    #: discards them. Reported so a run is self-checking instead of relying on a docstring.
    reads_assigned: int = 0
    collapsed: int = 0  # clonotypes absorbed into a parent
    reads_with_junction: int = 0   # Stage-1 reads carrying any junction
    reads_incomplete: int = 0      # ...of which dropped as truncated/out-of-frame/stop
    # Reads MOVED onto their parent clonotype by --min-junction-q, not discarded (0 with the gate
    # off). They are still counted -- in the parent -- so this never subtracts from `reads`.
    reads_low_quality: int = 0
    clonotypes_low_quality: int = 0  # clonotypes that gave up EVERY read to the quality gate
    # The quality-directed rescue (`--ec-mode amplicon|rnaseq`). Reads MOVE to the parent; a
    # candidate with no parent is kept, and its reads are `orphan_reads` -- reported, never lost.
    rescued_clonotypes: int = 0
    rescued_reads: int = 0
    orphan_clonotypes: int = 0
    orphan_reads: int = 0
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


def _quality_gate(df: pl.DataFrame, *, min_q: int, max_subs: int, require_vj: bool,
                  error_rate: float = 1e-3) -> tuple[pl.DataFrame, int, int]:
    """Move reads whose junction differs from its putative parent ONLY at low-quality bases.

    ⛔ **The read is REASSIGNED to the parent, never discarded.** A read that reached a complete
    junction came off a real rearrangement of that locus; the evidence says its differing base is a
    miscall, which is a statement about one base, not about whether the molecule existed. So the
    read's clonotype key is rewritten to the parent's and it is counted there. Dropping it would
    understate the parent's expression by exactly the reads the correction says belong to it.

    (Before this, the gate filtered the reads out and coverage assignment happened to realign most
    of them onto the parent anyway -- 325 of 330 on a TRA amplicon, 5 lost to the ``min_ov`` floor
    in :func:`_assign_coverage`. Reassignment makes that structural: the reads route through the
    exact-key pass, and the count is conserved by construction, not by an alignment succeeding.)

    The error model in :func:`_parents` sees abundance and nothing else, so a sequencing miscall
    and a real low-frequency variant are the same object to it and the only lever is
    ``error_rate`` — which trades them off globally and cannot separate them. Phred does separate
    them, because it is a different measurement: a miscall is a detector artifact and reads low-Q,
    while a real base — a true variant, or a template error made before the UMI — reads high-Q.
    Measured at the mismatching base over 310,559 real MIGEC windows: the two published spike-in
    variants sit at median Q 34-35 (16.7 / 17.6 % below Q30), the surrounding 1-substitution error
    cloud at **median Q 24** (54.3 % below Q30) and the 2-substitution cloud at median Q 6
    (91.1 %). ⚠ An earlier pass reported the 1-sub cloud at median 16; that extraction required an
    exact interior seed and so could only see mismatches near the junction's ends, which is its
    low-Q half. The separation is a 14-point median gap, not 19.

    ⛔ Only the MISMATCHING bases are evidence. A junction agreeing with its parent everywhere
    else says nothing about whether the one differing base is real, so gating on the junction's
    minimum or mean quality asks the wrong question and mostly measures read length.

    The parent is the most abundant substitution-only neighbour within ``max_subs`` (sharing V/J
    when ``require_vj``). A clonotype with no more-abundant neighbour is not gated at all: there is
    no hypothesis "this is a misread of X" to test. A read whose quality string is missing or too
    short is KEPT — absent evidence is not evidence of error.

    Returns ``(df with reassigned rows, reads reassigned, clonotypes emptied, sequence_id ->
    parent clonotype key for every moved read, vacated clonotype key -> its parent's key)``.
    """
    if "junction_quality" not in df.columns:
        # Raise, never degrade. A silently unapplied gate is indistinguishable from a gate that
        # applied and found nothing to drop, and this project has already shipped that failure
        # once (IgBLAST's missing `_gl.aux` produced a truth with no junctions and exit 0).
        raise ValueError(
            "--min-junction-q needs a `junction_quality` column, which this AIRR does not have. "
            "Re-run Stage 1 with `arda rnaseq map --junction-quality` (Stage 1 is the only place "
            "the FASTQ quality is still in hand).")
    jn = df["junction"].to_list()
    jq = [x or "" for x in df["junction_quality"].to_list()]
    v = [x or "" for x in df["v_call"].to_list()]
    j = [x or "" for x in df["j_call"].to_list()]
    loc = [x or "" for x in df["locus"].to_list()]

    counts: Counter = Counter(zip(loc, v, j, jn))
    keys = sorted(counts)                                  # deterministic clonotype order
    pos = {k: i for i, k in enumerate(keys)}
    seqs = [k[3] for k in keys]
    cnt = [counts[k] for k in keys]

    # Discriminating positions against the most abundant qualifying neighbour, per clonotype --
    # and which neighbour that was, since a gated read is reassigned onto it.
    disc: list[tuple[int, ...]] = [() for _ in keys]
    par: list[int | None] = [None] * len(keys)
    safe = [i for i, s in enumerate(seqs) if s and set(s) <= _DNA]
    if safe:
        index = seqtree.Index.build([seqs[i] for i in safe], alphabet="nt")
        params = seqtree.SearchParams(max_subs=max_subs, max_ins=0, max_dels=0,
                                      max_total_edits=max_subs, engine="seqtm")
        for ci in safe:
            best, best_count = None, -1
            for hit in index.search(seqs[ci], params):
                nj = safe[hit.ref_id]
                if nj == ci or hit.n_subs == 0 or len(seqs[nj]) != len(seqs[ci]):
                    continue
                if require_vj and (keys[nj][1] != keys[ci][1] or keys[nj][2] != keys[ci][2]):
                    continue
                if cnt[nj] <= cnt[ci]:                     # a parent is strictly more abundant
                    continue
                if cnt[nj] > best_count or (cnt[nj] == best_count and seqs[nj] < seqs[best]):
                    best, best_count = nj, cnt[nj]
            if best is not None:
                a, b = seqs[ci], seqs[best]
                disc[ci] = tuple(p for p in range(len(a)) if a[p] != b[p])
                par[ci] = best

    # One representative amino-acid junction per clonotype, so a reassigned read carries the
    # parent's `junction_aa` and not its own (the group's is taken with `.first()` downstream).
    have_aa = "junction_aa" in df.columns
    aa = df["junction_aa"].to_list() if have_aa else [None] * df.height
    key_aa: dict[int, object] = {}
    for r in range(df.height):
        key_aa.setdefault(pos[(loc[r], v[r], j[r], jn[r])], aa[r])

    cut = min_q + 33                                       # Phred+33; `>= min_q` keeps the boundary
    out_loc, out_v, out_j = df["locus"].to_list(), df["v_call"].to_list(), df["j_call"].to_list()
    out_jn, out_aa, out_q = df["junction"].to_list(), list(aa), df["junction_quality"].to_list()
    moved: Counter = Counter()
    moved_rows: list[tuple[int, int]] = []                 # (row, parent clonotype index)
    for r in range(df.height):
        ci = pos[(loc[r], v[r], j[r], jn[r])]
        d = disc[ci]
        if not d:
            continue
        q = jq[r]
        if len(q) != len(jn[r]):                           # no quality for this read: no evidence
            continue
        if not any(ord(q[p]) < cut for p in d):
            continue
        # ⛔ AND the parent must be able to have PRODUCED this read. "More abundant" is far too
        # weak: one extra read made anything within `max_subs` a parent, and at 3 substitutions
        # that is a hypothesis nothing supports -- which is why the TRA amplicons, whose short
        # junctions (median 42 nt) have many 3-sub neighbours, lost 0.44 % and 1.39 % of their
        # reads while TRB gained.
        #
        # ⚠ The test uses THIS READ'S OWN PHRED, not a global `error_rate`. Using the global rate
        # was tried and is wrong in the other direction: it makes the gate a strict subset of the
        # abundance model, which is the thing the gate exists to reach past. Measured, MIGEC at
        # 1e-5 went 1,633 -> 1,633 error clonotypes (from 127) because a 2-sub parent scores
        # 293,327 * (4.8e-4)^2 = 0.068 < 1. Phred is the independent evidence here, so it is what
        # the plausibility is computed from.
        # ⚠ `dpos`, not `pos` -- `pos` is the clonotype-key index dict two lines up, and shadowing
        # it made the SECOND read of the loop fail with "'int' object is not subscriptable".
        p_read = 1.0
        for dpos in d:
            p_read *= 10.0 ** (-(ord(q[dpos]) - 33) / 10.0)
        if cnt[par[ci]] * p_read < cnt[ci]:
            continue
        out_loc[r], out_v[r], out_j[r], out_jn[r] = keys[par[ci]]
        out_aa[r] = key_aa.get(par[ci])
        # Its quality describes the junction it NO LONGER carries. A same-length quality string
        # belonging to a different junction is the one corruption nothing downstream can detect,
        # so it is blanked rather than carried over.
        out_q[r] = ""
        moved[ci] += 1
        moved_rows.append((r, par[ci]))

    # ⛔ Rewriting the row is not enough on its own. Coverage assignment re-derives every read's
    # clonotype from the UNFILTERED Stage-1 frame by its original ``(locus, v, j, junction)`` key,
    # so a moved read whose old clonotype still exists (because its other reads passed the gate)
    # would be routed straight back to it and the move silently undone. Name the moved reads by
    # ``sequence_id`` so that pass can override them read-by-read. The key a read is rewritten to
    # necessarily survives -- that read is now in it -- so no chain walk is needed.
    sids = df["sequence_id"].to_list()
    moved_sid = {sids[r]: keys[p_i] for r, p_i in moved_rows}

    cols = [pl.Series("locus", out_loc), pl.Series("v_call", out_v), pl.Series("j_call", out_j),
            pl.Series("junction", out_jn), pl.Series("junction_quality", out_q)]
    if have_aa:
        cols.append(pl.Series("junction_aa", out_aa))
    # Emptied = gave up every read AND received none back (a clonotype that took reads from a
    # child of its own is still there).
    received = {p_i for _, p_i in moved_rows}
    emptied_keys = {keys[ci]: keys[par[ci]] for ci, n in moved.items()
                    if n == cnt[ci] and ci not in received}
    return df.with_columns(cols), sum(moved.values()), len(emptied_keys), moved_sid, emptied_keys


def _root(i: int, parent: list[int | None]) -> int:
    while parent[i] is not None:
        i = parent[i]  # type: ignore[assignment]
    return i


def _binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), summed as 1 - CDF(k-1).

    ⛔ Computed in LOG space, and it has to be. The obvious form -- ``comb(n, i) * p**i *
    (1-p)**(n-i)`` -- raises ``OverflowError: int too large to convert to float`` the moment a
    clonotype is big: ``comb(293327, i)`` is an exact Python int with tens of thousands of digits,
    and multiplying it by a float converts it first. Measured on a real 302,172-read MIGEC library
    with one 293 k-read clone, ``--error-method binom`` crashed after 190 s. The overflow is not
    the whole story either: that formulation is O(k) with k = the child's read count, so it was
    also spending three minutes to reach the crash.

    So: `lgamma` for the log-binomial coefficient, `log1p(-p)` for the tail factor (accurate when
    p is tiny, which is exactly this caller's regime), and an early exit once the running CDF is
    within float epsilon of 1 -- the terms fall off geometrically past the mean, so a clone with
    293 k reads no longer costs 293 k iterations.
    """
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    from math import exp, lgamma, log, log1p
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    lp, lq = log(p), log1p(-p)
    lgn = lgamma(n + 1)
    cdf = 0.0
    for i in range(k):
        term = exp(lgn - lgamma(i + 1) - lgamma(n - i + 1) + i * lp + (n - i) * lq)
        cdf += term
        # Past the mean the terms decay geometrically; once the complement is below float
        # resolution the answer is 0.0 and every further term is noise.
        if cdf >= 1.0 - 1e-15 and i > n * p:
            return 0.0
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
    # Only the columns this function reads. A Stage-1 AIRR has 83, and `to_list()` on the
    # rest builds Python str objects for `sequence_alignment`, `germline_alignment` and
    # every region sequence -- measured 2.42 KB/row against 0.41 KB/row for the columns
    # actually used, i.e. ~2 KB wasted per mapped read (~7 GB at SRR5233639's full depth).
    # `col()` below already tolerates an absent column, so restricting the dict is safe.
    _USED = ("c_call", "j_call", "locus", "rev_comp", "sequence", "v_call")
    cols = {c: raw[c].to_list() for c in raw.columns if c in _USED}
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
        best_ri, best_ov, best_mm, best_lo, best_hi, seen = None, 19, 0, 0, 0, set()
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
                # Strictly worse overlap cannot win, so skip the mismatch scan. An EQUAL overlap
                # must still be scored: see `_assign_coverage` for why arrival order is not a
                # tie-break.
                if ov < best_ov or ov < 20:
                    continue
                mm = sum(1 for kk in range(lo, hi) if s[kk] != jr[kk + d])
                if mm <= 0.12 * ov and (best_ri is None or ov > best_ov or mm < best_mm):
                    best_ri, best_ov, best_mm, best_lo, best_hi = ri, ov, mm, lo + d, hi + d
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
            # A parent must be strictly MORE abundant than its child. `_parents` gets this for
            # free -- `count[parent] * p_err >= count[child]` with `p_err < 1` forces counts to
            # increase along parent pointers, which is the module's stated no-cycle proof. The
            # depth test above carries no such ordering: it asks only whether the child allele's
            # depth is consistent with sequencing error of the parent's, and that can be true in
            # BOTH directions for a pair of similar-abundance neighbours. Two clonotypes then
            # become each other's parent and `_root` walks the 2-cycle forever -- an unbounded
            # hang, not a wrong number, on `--error-method binom|betabinom`.
            if ok and span_counts[nj] > span_counts[ci] and span_counts[nj] > best_count:
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
    override: dict[str, int] | None = None,
    aliases: list[tuple[str, str, int]] | None = None,
    root_counts: list[int] | None = None,
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

    ``override`` maps a ``sequence_id`` straight to a root position, ahead of every key lookup: it
    carries the reads ``_quality_gate`` moved onto a parent, whose key in ``raw`` still names the
    error clonotype they were moved off.

    ``aliases`` are extra ``(junction, locus, root position)`` targets for the ALIGNMENT pass. A
    junction the quality gate vacated is still real sequence that partial reads legitimately cover
    -- it just is not its own clonotype any more -- so it stays in the index, pointing at the
    parent. Without it a read whose only >= ``min_ov`` overlap was with the vacated junction goes
    unassigned, which loses reads the correction says belong to the parent.

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
    # Only the columns this function reads. A Stage-1 AIRR has 83, and `to_list()` on the
    # rest builds Python str objects for `sequence_alignment`, `germline_alignment` and
    # every region sequence -- measured 2.42 KB/row against 0.41 KB/row for the columns
    # actually used, i.e. ~2 KB wasted per mapped read (~7 GB at SRR5233639's full depth).
    # `col()` below already tolerates an absent column, so restricting the dict is safe.
    _USED = ("c_call", "j_call", "junction", "locus", "rev_comp", "sequence",
             "sequence_id", "v_call")
    cols = {c: raw[c].to_list() for c in raw.columns if c in _USED}

    def col(name):
        return cols.get(name, [None] * n)

    seqc, rcc, jnc, sidc = col("sequence"), col("rev_comp"), col("junction"), col("sequence_id")
    locc, vc, jc, cc = col("locus"), col("v_call"), col("j_call"), col("c_call")
    done: set[str] = set()

    # Pass 1: exact clonotype-key match (spanning + rescued reads; authoritative).
    override = override or {}
    for i in range(n):
        sid = sidc[i]
        if sid in done:                              # a rescued read appears in mapped + assembled rows
            continue
        rp = override.get(sid)
        if rp is None:
            rp = exact.get(((locc[i] or ""), (vc[i] or ""), (jc[i] or ""), (jnc[i] or "")))
        if rp is not None:
            assigned[rp].append(sid)
            done.add(sid)

    # Pass 2: align the rest (partial V-side / J-side reads that never reached a complete junction).
    # Targets are the roots, then the vacated junctions; `tgt_root` maps a target back to the root
    # position its reads are credited to (identity for a root, the parent for an alias).
    tgt_jn = list(root_jn) + [a[0] for a in (aliases or ())]
    tgt_loc = list(root_loc) + [a[1] for a in (aliases or ())]
    tgt_root = list(range(len(root_jn))) + [a[2] for a in (aliases or ())]
    # ⛔ `cap` bounds how many roots one germline-shared k-mer may name, and WHICH ones it keeps is
    # therefore load-bearing. Keeping the first `cap` in target order made the survivors depend on
    # the root list's order, so any change to the root set silently re-shuffled them: merging call
    # splits under `--clonotype-key junction` cost 3 reads of 43,475 that had been placed by this
    # pass, none of which had a junction of its own to fall back on. Insert in DESCENDING ABUNDANCE
    # (then sequence, for a total order): the roots a partial read is most likely to belong to are
    # the ones that survive the cap, and the choice no longer depends on upstream grouping.
    # ⛔⛔ AND A ROOT OUTRANKS AN ALIAS, unconditionally. An alias is a FALLBACK -- a junction that
    # is no longer a clonotype, kept only so partial reads that covered it still reach the parent --
    # but it is ordered by its PARENT's abundance, which is high by construction. So aliases sorted
    # to the FRONT and evicted genuine low-abundance roots from the cap, and every partial read whose
    # only home was such a root went unassigned. Measured at full depth on a TRA amplicon
    # (SRR5233636, --ec-mode accurate, 23,360 aliases against 36,587 roots), all three arms emitting
    # an IDENTICAL clonotype table so this is purely read assignment:
    #     cap 64, aliases ordered by abundance   1,812,740   -25,473 vs --ec-mode fast
    #     cap 64, aliases OFF                    1,838,181       -32
    #     cap 1024, aliases ordered by abundance 1,869,556   +31,343
    # i.e. the alias mechanism -- added to rescue 5 reads of 9,208 on Ramos -- was costing 25,441
    # reads here, while the same aliases GAIN 31,343 once the index is big enough to hold both. They
    # are worth keeping; they must simply not outrank the roots. Roots first, then aliases into
    # whatever slots remain.
    order = sorted(range(len(tgt_jn)),
                   key=lambda ri: (ri >= len(root_jn),
                                   -(root_counts[tgt_root[ri]] if root_counts else 0),
                                   tgt_jn[ri], ri))
    index: dict[str, list[tuple[int, int]]] = defaultdict(list)
    placed = [False] * len(tgt_jn)
    for ri in order:
        jn = tgt_jn[ri]
        for p in range(len(jn) - k + 1):
            lst = index[jn[p:p + k]]
            if len(lst) < cap:                       # bound the germline-shared k-mers
                lst.append((ri, p))
                placed[ri] = True
    # ⛔ A target whose every k-mer was already full gets ZERO postings and is UNREACHABLE -- no read
    # can be assigned to it by this pass at all, however well it matches. That is where the cap's
    # read loss comes from, and raising `cap` is the expensive way to fix it: measured on a TRA
    # amplicon (SRR5233636, --ec-mode fast, identical clonotype table at every cap, so this is purely
    # read assignment) each doubling buys ~6,300 reads and costs ~1.5x wall --
    #     cap 64 1,838,213 / 207.4 s | 128 1,844,359 / 312.1 s | 256 1,850,917 / 485.9 s
    #     cap 1024 1,867,904 (+29,691 vs 64)
    # Instead give every unreachable target ONE slot, in its LEAST-LOADED k-mer -- the k-mer it
    # shares with the fewest others, i.e. the most specific one it has. The overflow is bounded by
    # the number of such targets and `cap` is unchanged.
    # ⛔ ROOTS ONLY. Forcing an ALIAS back in re-creates the very bug roots-first exists to fix --
    # the alias then competes with the rare root whose junction it carries, and at equal overlap and
    # equal mismatches the read is credited to the alias's abundant parent instead. Caught by
    # test_an_abundant_alias_must_not_steal_the_kmer_slots_of_a_rare_root. The guarantee wanted here
    # is that every real CLONOTYPE is reachable, not that every vacated junction is.
    for ri in range(len(root_jn)):
        jn = tgt_jn[ri]
        if placed[ri] or len(jn) < k:
            continue
        best_p, best_load = 0, None
        for p in range(len(jn) - k + 1):
            load = len(index[jn[p:p + k]])
            if best_load is None or load < best_load:
                best_p, best_load = p, load
        index[jn[best_p:best_p + k]].append((ri, best_p))
    for i in range(n):
        sid = sidc[i]
        if sid in done:
            continue
        seq = seqc[i]
        if not seq or len(seq) < k:
            continue
        loc = (locc[i] or "")[:3] or _gene3(vc[i]) or _gene3(jc[i]) or _gene3(cc[i])
        s = reverse_complement(seq) if str(rcc[i]).upper() in ("T", "TRUE", "1") else seq
        best_ri, best_ov, best_mm, seen = None, min_ov, 0, set()
        L = len(s)
        for rp in range(0, L - k + 1):
            for (ri, jp) in index.get(s[rp:rp + k], ()):
                if loc and tgt_loc[ri] != loc:
                    continue
                d = jp - rp
                if (ri, d) in seen:
                    continue
                seen.add((ri, d))
                jr = tgt_jn[ri]
                lo, hi = (-d if d < 0 else 0), min(L, len(jr) - d)
                ov = hi - lo
                # ⛔ `ov < best_ov`, NOT `ov <= best_ov`. Equal overlap is COMMON and was being
                # resolved by encounter order: the first root to reach a read kept it, and the
                # mismatch count computed right below was thrown away.
                #
                # That is not hypothetical. When a phantom clonotype (a junction window slid 10 nt
                # into V framework, arda <= 2.9.0) competed with the true Jurkat clone, ~47 % of the
                # 5,758 reads it took were exact 48-vs-48 overlap ties -- and on those the losing
                # TRUE root matched with 0 mismatches against the phantom's 1. Reads span the
                # junction end to end, so ties are the normal case, not the corner case.
                #
                # Overlap first, then FEWER MISMATCHES. Both were already computed.
                if ov < best_ov or ov < min_ov:
                    continue
                budget = max_mm * ov
                mm = 0
                for kk in range(lo, hi):
                    if s[kk] != jr[kk + d]:
                        mm += 1
                        if mm > budget:
                            break
                if mm <= budget and (best_ri is None or ov > best_ov or mm < best_mm):
                    best_ov, best_mm, best_ri = ov, mm, ri
        if best_ri is not None:
            assigned[tgt_root[best_ri]].append(sid)
            done.add(sid)
    return assigned


def _gene3(x: str | None) -> str:
    g = (x or "").split(",")[0].split("(")[0].strip()
    return g[:3] if g[:3] in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD") else ""


#: Clonotype D columns. The four calls/supports, then everything a D-D MARKUP consumer needs to
#: cut the junction up: where the V stops templating it, where the J starts, the D span(s), and
#: the non-templated stretches between them. ⛔ All coordinates are 1-based closed in JUNCTION
#: space -- the clonotype table has no read -- with -1 meaning "not located". Without them the
#: table named a second D and gave no way to find it; `d2_call` alone is not markup.
_D_STR_COLS = ("d_call", "d2_call", "d_support", "d2_support", "np1", "np2", "np3")
_D_INT_COLS = ("v_sequence_end", "d_sequence_start", "d_sequence_end",
               "d2_sequence_start", "d2_sequence_end", "j_sequence_start")


def _clonotype_d(out: pl.DataFrame, organism: str,
                 d_max_evalue: float | None = None) -> list[pl.Series]:
    """D (and tandem D-D) per clonotype, mapped into its error-corrected junction.

    Reads carry their own ``d_call``, but a read's D is called on a sequencing-error copy of
    the junction and, for a long CDR3, on a read that does not span the D at all. The clonotype
    is the first place the junction is both complete and corrected, and D is a deterministic
    function of ``(junction, v_call, j_call)`` -- so call it once here rather than voting over
    reads. Costs one gapless alignment per clonotype, not per read.

    ``d_max_evalue`` overrides the shipped E-value gate; see
    :func:`arda.annotate.transfer._map_d`.

    VJ loci and organisms without D germlines come back empty, as does an unresolvable V/J.
    """
    from ..annotate.dmap import map_d_junction

    cols: dict[str, list] = {c: [] for c in _D_STR_COLS + _D_INT_COLS}
    for jn, vc, jc in zip(out["junction"], out["v_call"], out["j_call"]):
        try:
            call = map_d_junction(jn or "", vc or "", jc or "", organism,
                                  d_max_evalue=d_max_evalue)
        except (KeyError, ValueError):        # unknown allele / organism without anchors
            call = None
        for c in _D_STR_COLS:
            cols[c].append(getattr(call, c) if call else "")
        for c in _D_INT_COLS:
            cols[c].append(getattr(call, c) if call else -1)
    return ([pl.Series(c, cols[c], dtype=pl.Utf8) for c in _D_STR_COLS]
            + [pl.Series(c, cols[c], dtype=pl.Int32) for c in _D_INT_COLS])


def correct_airr(
    airr_tsv: str | Path,
    output: str | Path,
    *,
    organism: str = "human",
    map_d: bool = True,
    d_max_evalue: float | None = None,
    max_subs: int = 3,
    max_indel: int = 0,
    error_rate: float = 0.001,
    indel_rate: float = 0.001,
    require_vj: bool = True,
    error_method: str | None = None,
    ec_mode: str = "fast",
    clonotype_key: str = "full",
    min_junction_q: int | None = None,
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
        map_d: append the D columns (``_D_STR_COLS`` + ``_D_INT_COLS``), called once per
            clonotype on its corrected junction (see :func:`_clonotype_d`). Default ``True``.
        d_max_evalue: E-value gate on the D call(s); ``None`` keeps the shipped 0.2. Lower is
            stricter -- 0.01 is the band where D agrees .9985 with IgBLAST on a TRB amplicon.
        output: corrected clonotype table TSV (``junction``, ``junction_aa``, ``v_call``,
            ``j_call``, ``c_call``, ``locus``, ``duplicate_count``, ``consensus_count``, and
            with ``map_d`` the D columns), sorted
            by abundance. A clonotype is keyed by ``(locus, v_call, j_call, junction)``. Per the
            AIRR schema, ``duplicate_count`` is the number of READS supporting the clonotype (both
            paired mates of a molecule count) and ``consensus_count`` is the number of distinct
            fragment consensuses (the two mates of one molecule are one consensus). ``c_call`` is the
            clonotype's dominant isotype CLASS (from ``c_class``: IGHG, IGHA, ...), preferring a
            resolved class over the ambiguous ``IGHC``; empty when no read carried a constant call.
        max_subs: max substitutions between an error child and its parent (seqtree neighbour search).
            This is a SEARCH RADIUS, not a threshold -- the accept/reject decision is the
            length-scaled probability model above, so widening it only lets the model SEE parents
            it would already have accepted. The default was 2 through 2.9.0, which truncated the
            search below what the model would take on a deep clone: on the two monoclonal cell
            lines in the arda-benchmark set 2 -> 3 collapses Jurkat 74 -> 57 clonotypes and Raji
            91 -> 58, while a polyclonal mouse spleen (7,942) and an oligoclonal B-LCL (13) are
            UNCHANGED at 2, 3 and 4 -- the model refuses those collapses on abundance regardless
            of radius. It saturates at 3 (4 gives the same four numbers), so 3 is the default.
        max_indel: max inserted/deleted bases searched for indel error children (default 0). A 1-2 bp
            instrument indel is a frameshift and is already dropped by ``complete_only``, so on
            complete junctions the indel search only costs time (~160x slower) and collapses nothing;
            a multi-base in-frame SHM indel costs ``(indel_rate*L)**len`` and is kept as a real
            clonotype either way. Set it > 0 only with ``--all-junctions`` (frameshift indels kept).
        error_rate: per-BASE substitution error rate (~Phred 30 = 0.001). The per-substitution
            collapse probability is length-scaled, ``error_rate * junction_len``, so the default
            reproduces vdjtools' ~1/20 at a 45 nt (15 aa) junction and scales for other lengths.
        indel_rate: per-BASE indel error rate (instrument-dependent; default 0.001, length-scaled).
        ec_mode: knob preset, ``"fast"`` (default, = today's shipped behaviour) or
            ``"accurate"``. See :data:`EC_MODES`; an explicitly passed ``error_method`` /
            ``min_junction_q`` overrides it.
        min_junction_q: reassign onto its parent a read whose junction differs from that parent
            below this Phred score (:func:`_quality_gate`). ``0`` disables it. Needs the
            ``junction_quality`` column from ``map --junction-quality`` and RAISES without it.
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
    if ec_mode not in EC_MODES:
        raise ValueError(f"ec_mode must be one of {sorted(EC_MODES)}, got {ec_mode!r}")
    if clonotype_key not in CLONOTYPE_KEYS:
        raise ValueError(f"clonotype_key must be one of {list(CLONOTYPE_KEYS)}, got {clonotype_key!r}")
    preset = EC_MODES[ec_mode]
    regime = preset.get("regime", "fast")
    if error_method is None:
        error_method = preset["error_method"]
    if min_junction_q is None:
        min_junction_q = preset["min_junction_q"]
    for name, val in (("error_rate", error_rate), ("indel_rate", indel_rate)):
        # p_err < 1 keeps counts strictly increasing along parent pointers (no cycles); p_err == 0
        # would make a single mismatch collapse anything, p_err >= 1 never collapses.
        if not 0.0 < val < 1.0:
            raise ValueError(f"{name} must be in (0, 1), got {val}")
    if error_method not in ("simple", "binom", "betabinom"):
        raise ValueError(f"error_method must be simple|binom|betabinom, got {error_method!r}")

    stage = Stage()
    output = Path(output)
    raw = _read_airr(airr_tsv)
    if extra_airr is not None:
        extra = _read_airr(extra_airr)
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
    # ⛔ The KEY is decided BEFORE any correction. The quality gate names a read's parent by its
    # clonotype key, so canonicalising afterwards leaves those names pointing at keys that no
    # longer exist -- 2 of Jurkat's 14,531 reads went missing that way, which is small enough to
    # have shipped unnoticed and is still a violation of the invariant.
    if clonotype_key == "junction":
        # Collapse CALL SPLITS by giving every read of a (locus, junction) that junction's majority
        # (v_call, j_call) before grouping. Rewriting the labels rather than shortening the group
        # key keeps everything downstream -- coverage's exact-key pass, the read map, the emitted
        # columns -- working on a 4-tuple, so this is a relabel, not a second code path.
        # ⛔ Ties break on count then lexicographically, never on row order: `correct`'s
        # reproducibility rests on exactly this kind of decision being total.
        votes: dict[tuple[str, str], Counter] = defaultdict(Counter)
        for loc, vv, jj, jn in zip(df["locus"].to_list(), df["v_call"].to_list(),
                                   df["j_call"].to_list(), df["junction"].to_list()):
            votes[(loc or "", jn or "")][(vv or "", jj or "")] += 1
        best = {k: min(sorted(c), key=lambda pair: (-c[pair], pair)) for k, c in votes.items()}
        # ⛔ `raw` must be relabelled TOO, not just `df`. Coverage assignment re-derives every
        # read's clonotype from the UNFILTERED frame by its own (locus, v, j, junction) key, so a
        # read whose label the canonicalisation changed would miss the exact-key pass and fall to
        # the aligner or go unassigned. Measured when this was missing: 128 of Jurkat's 14,531
        # reads vanished -- the same read-conservation leak the quality gate had, one flag later.
        rl, rj = raw["locus"].to_list(), raw["junction"].to_list()
        rv, rjc = raw["v_call"].to_list(), raw["j_call"].to_list()
        raw = raw.with_columns(
            pl.Series("v_call", [best.get((a or "", b or ""), (c, d))[0]
                                 for a, b, c, d in zip(rl, rj, rv, rjc)]),
            pl.Series("j_call", [best.get((a or "", b or ""), (c, d))[1]
                                 for a, b, c, d in zip(rl, rj, rv, rjc)]),
        )
        df = df.with_columns(
            pl.Series("v_call", [best[(loc or "", n or "")][0]
                                 for loc, n in zip(df["locus"].to_list(), df["junction"].to_list())]),
            pl.Series("j_call", [best[(loc or "", n or "")][1]
                                 for loc, n in zip(df["locus"].to_list(), df["junction"].to_list())]),
        )

    n_low_q = n_clono_low_q = 0
    moved_sid: dict[str, tuple[str, str, str, str]] = {}
    emptied_keys: dict[tuple[str, str, str, str], tuple[str, str, str, str]] = {}
    if min_junction_q > 0:
        df, n_low_q, n_clono_low_q, moved_sid, emptied_keys = _quality_gate(
            df, min_q=min_junction_q, max_subs=max_subs, require_vj=require_vj,
            error_rate=error_rate)

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
                           reads_incomplete=n_incomplete,
                           reads_low_quality=n_low_q,
                           clonotypes_low_quality=n_clono_low_q)
    parent = _parents(junctions, span_counts, v, j, max_subs=max_subs, max_indel=max_indel,
                      error_rate=error_rate, indel_rate=indel_rate, require_vj=require_vj)
    if error_method != "simple":
        parent = _error_pileup(raw, junctions, v, j, locus, parent, span_counts,
                               error_rate=error_rate, indel_rate=indel_rate, method=error_method)

    # QUALITY-DIRECTED RESCUE (arda.rnaseq.denoise). The abundance model above is only valid where
    # a ladder of observed intermediates exists -- measured at k <= 2, with 3 as headroom. Beyond
    # that there is no ladder (0 of 13 clonotypes at k=4 on Jurkat have any observed k-1 neighbour,
    # against 0.0019 predicted) and `--min-junction-q` cannot see the class either, because a
    # clonotype with no neighbour within `max_subs` has no discriminating base to judge. What IS
    # true of it is that its reads are bad. So: only clonotypes whose reads are measurably bad, only
    # onto a much more abundant parent, and a candidate with no parent KEEPS ITS READS.
    rescue_report = None
    if regime and REGIMES[regime].enabled() and "junction_quality" not in df.columns:
        # ⛔ Raise, never degrade -- the same rule `_quality_gate` enforces above, for the same
        # reason. Skipping the rescue silently produces a report indistinguishable from a rescue
        # that ran and found nothing (rescued/orphan counters all 0) and a table byte-identical to
        # `--ec-mode fast`. `rnaseq run` wires the column up itself, but the standalone `correct`
        # entry point cannot, so an AIRR mapped without --junction-quality reached here and the
        # whole point of --ec-mode amplicon|rnaseq was quietly dropped.
        raise ValueError(
            f"--ec-mode {regime} needs a `junction_quality` column, which this AIRR does not have. "
            "Re-run Stage 1 with `arda rnaseq map --junction-quality` (Stage 1 is the only place "
            "the FASTQ quality is still in hand), or use --ec-mode fast.")
    if regime and REGIMES[regime].enabled() and "junction_quality" in df.columns:
        rq = read_quality(df["junction"].to_list(),
                          [x or "" for x in df["junction_quality"].to_list()])
        cq = clonotype_quality(list(zip([x or "" for x in df["locus"].to_list()],
                                        [x or "" for x in df["v_call"].to_list()],
                                        [x or "" for x in df["j_call"].to_list()],
                                        df["junction"].to_list())), rq)
        clono_q = [cq.get((locus[i], v[i], j[i], junctions[i]), -1.0)
                   for i in range(len(junctions))]
        # Only clonotypes the abundance model did NOT already claim: a read routed twice is a bug,
        # and re-deciding a parent the ladder already justified would be strictly worse evidence.
        free = set(i for i in range(len(junctions)) if parent[i] is None)
        gated_q = [clono_q[i] if i in free else -1.0 for i in range(len(junctions))]
        # ⛔ PER LOCUS. `quality_rescue` groups candidate parents by junction LENGTH and nothing
        # else -- no locus guard in either `_nearest_py` or the C++ `nearest_more_abundant` -- while
        # `amplicon` opens the radius to 12 substitutions. A rearrangement of a DIFFERENT locus is
        # not a sequencing error of this one, so those merges are wrong by construction. Measured on
        # SRR5233636 with `amplicon`: of 9,025 rescues, 3 crossed the locus -- 1-read TRB clonotypes
        # absorbed into abundant TRA clonotypes at 11-12 substitutions.
        #
        # Partitioning here rather than adding a group argument to the extension keeps the C++ ABI
        # (and the prebuilt wheels) untouched, and it also bounds the search: the scan is quadratic
        # within a length bin, and on a 397k-clonotype library the single-threaded whole-table
        # version ran over 39 minutes of CPU against ~5 for a 49k one.
        rescued: list[int | None] = [None] * len(junctions)
        rescue_report = DenoiseReport()
        by_locus: dict[str, list[int]] = defaultdict(list)
        for i in range(len(junctions)):
            by_locus[locus[i]].append(i)
        for _loc, idx in sorted(by_locus.items()):
            sub, rep = quality_rescue([junctions[i] for i in idx], [span_counts[i] for i in idx],
                                      [gated_q[i] for i in idx], REGIMES[regime])
            for k, pi in enumerate(sub):
                if pi is not None and pi >= 0:
                    rescued[idx[k]] = idx[pi]              # back to global clonotype indices
            rescue_report.lowq_clonotypes += rep.lowq_clonotypes
            rescue_report.rescued_clonotypes += rep.rescued_clonotypes
            rescue_report.rescued_reads += rep.rescued_reads
            rescue_report.orphan_clonotypes += rep.orphan_clonotypes
            rescue_report.orphan_reads += rep.orphan_reads
            rescue_report.quality_available |= rep.quality_available
            for d, c in rep.dist.items():
                rescue_report.dist[d] = rescue_report.dist.get(d, 0) + c
        for i, pi in enumerate(rescued):
            # Never overwrite an abundance parent, and never create a cycle: the rescue target is
            # strictly more abundant (enforced in `nearest_more_abundant`), which is the same
            # no-cycle argument `_parents` relies on.
            if pi is not None and parent[i] is None:
                parent[i] = pi
        report.rescued_clonotypes = rescue_report.rescued_clonotypes
        report.rescued_reads = rescue_report.rescued_reads
        report.orphan_clonotypes = rescue_report.orphan_clonotypes
        report.orphan_reads = rescue_report.orphan_reads

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
        # A clonotype the gate emptied is no longer in `exact`, so any OTHER read in `raw` still
        # carrying its key -- an INCOMPLETE-junction read, which never reached Stage 2 at all --
        # would lose its only anchor and go unassigned. (Measured: exactly one such read on Jurkat,
        # and it has no `sequence`, so the alignment pass cannot rescue it either.) Point the
        # vacated key at the parent instead, following chains of emptied clonotypes.
        for old_k, par_k in emptied_keys.items():
            seen = set()
            while par_k in emptied_keys and par_k not in seen:
                seen.add(par_k)
                par_k = emptied_keys[par_k]
            if par_k in exact:
                exact.setdefault(old_k, exact[par_k])
        # A quality-gated read is routed by its NEW clonotype, not by the key still sitting on it
        # in `raw` -- otherwise coverage hands it back to the error clonotype it was moved off.
        override = {sid: exact[k] for sid, k in moved_sid.items() if k in exact}
        # ...and the junctions the gate vacated stay in the ALIGNMENT index, pointing at the parent.
        # A partial read carrying no junction of its own can have covered only that sequence; with
        # the target gone it overlaps no surviving root by `min_ov` and is dropped. (Measured: 5 of
        # 9,208 on Ramos, every one of them a read with no junction at all.)
        aliases = [(old_k[3], old_k[0], exact[old_k]) for old_k in emptied_keys if old_k in exact]
        read_sets = _assign_coverage(raw, [junctions[i] for i in roots], [locus[i] for i in roots],
                                     exact, override=override, aliases=aliases,
                                     root_counts=[span_counts[i] for i in roots])
    else:
        read_sets = [agg_reads[i] for i in roots]
    dup = [len(rs) for rs in read_sets]                                 # AIRR duplicate_count: reads
    cons = [len({_strip_mate(x) for x in rs}) for rs in read_sets]      # AIRR consensus_count: fragments
    # ⛔ THE read-conservation invariant, reported so it is checkable on a real run. `reads` above
    # cannot move (it is counted before correction), so a benchmark comparing it across --ec-mode
    # sees 0 everywhere and reads that as conservation. Measured full-depth on the golden set with
    # only `reads` available: SRR5233636 fast 1,838,213 -> accurate 1,812,745 (-25,468) went unseen.
    report.reads_assigned = sum(dup)

    def _dominant_ccall(read_list: list[str]) -> str:
        # A clonotype's isotype = the dominant RESOLVED class over its fragments' constant mates.
        # `isotype_class` emits the generic `IGHC` only on cross-class ambiguity (IGHG1,IGHM), so
        # report IGHC only if NO read resolves -- a handful of ambiguous reads must not outvote it.
        # One vote per FRAGMENT, not per assigned mate. `read_list` holds sequence_ids, so a
        # fragment whose two mates were both assigned used to contribute its calls twice while a
        # fragment with one assigned mate contributed once -- weighting the vote by assigned
        # mates rather than by molecules, exactly as the first line of this comment says it must
        # not. A 1-fragment minority could then outvote a 2-fragment majority.
        calls: list[str] = []
        for frag in dict.fromkeys(_strip_mate(sid) for sid in read_list):
            calls.extend(frag_iso.get(frag, ()))
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
        out = out.with_columns(_clonotype_d(out, organism, d_max_evalue))
    out.write_csv(output, separator="\t", quote_style="never")

    if read_map is not None:
        rows = [(rid, junctions[roots[r]]) for r in range(len(roots)) for rid in read_sets[r]]
        pl.DataFrame(rows, schema=["sequence_id", "junction"], orient="row").write_csv(
            read_map, separator="\t", quote_style="never")
    stage.finish(report)
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report
