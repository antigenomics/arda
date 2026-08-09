"""Quality-aware clonotype denoising — the framework Stage 2's error model plugs into.

What this module is for
-----------------------
``correct``'s abundance model asks one question: could the parent have produced this many misreads
by chance? That question is **only answerable for near neighbours**, and the benchmark round 22
measured exactly where it stops being answerable.

On a monoclonal T line (Jurkat, dominant TRB clone 9,932 reads over 48 nt), clonotypes at Hamming
distance *k* from the dominant clone:

===  ===========  ==================================  =============  ==============================
 k   clonotypes   expected at ``error_rate`` 1e-3     median mean-Q  has an OBSERVED (k-1) neighbour
===  ===========  ==================================  =============  ==============================
 1          108                              476.74            31.4  --
 2           82                               11.20            25.2  **82 / 82**
 3           28                                0.17            24.2  3 / 28
 4           13                              0.0019            24.0  **0 / 13**
 5           18                              0.0000            20.1  **0 / 18**
>=6          14                              0.0000       16.5-18.2  **0 / 14**
===  ===========  ==================================  =============  ==============================

Two regimes, and they need different evidence:

* **k <= 2 (3 as headroom) is a LADDER.** Independent per-base errors accumulate through observed
  intermediates -- every one of the 82 two-substitution variants has an observed one-substitution
  neighbour on its path to the parent -- and the observed counts are within an order of magnitude
  of the binomial prediction. The abundance model is valid here and chain collapse walks it.
* **k >= 4 is a CLIFF.** Zero intermediates at every k, and the model predicts 0.0019 clonotypes at
  k = 4 where 13 are observed. These are not accumulated substitutions; they are single reads whose
  whole junction window is unreliable, and the quality says so independently (median mean junction
  Phred falls monotonically 31.4 -> 16.5; the k >= 5 class is 100 % sub-Q30 against the dominant
  clone's 5.9 %).

So widening ``--max-subs`` to reach the cliff class "works" on Jurkat (53 -> 11 clonotypes) **for
the wrong reason**: the abundance test it applies has probability 0 to every printed digit there.
This module reaches that class on the evidence that actually distinguishes it -- read quality --
and never on abundance alone.

⛔ The invariant
---------------
**Nothing here discards a read.** A read that reached a complete junction came off a real
rearrangement of that locus; deciding its junction carries a miscall is a statement about the
bases, not about whether the molecule existed. Every function returns *parent assignments*, and a
clonotype with no parent is KEPT, never emptied. The caller asserts that the sum of
``duplicate_count`` does not fall.

⛔ And the reason the cliff class cannot simply be deleted: measured on a polyclonal, hypermutated
repertoire (IGH_repertoire, 31,943 clonotypes), a whole-junction mean-Q floor at 30 would strand
**3.70 %** of all junction-bearing reads with no parent to inherit them (47 % of everything it
removes), against 0.148 % at a floor of 20. That is why the floor here is a *rescue radius* rather
than a filter: reads move to a parent when one exists and stay put when one does not.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median

__all__ = ["DenoiseParams", "DenoiseReport", "REGIMES", "read_quality", "clonotype_quality",
           "quality_rescue"]

try:                                              # optional: built by scikit-build-core
    from .. import _denoise as _cpp
except ImportError:                               # pragma: no cover - source checkout without ext
    _cpp = None


@dataclass(frozen=True)
class DenoiseParams:
    """One regime's settings. See :data:`REGIMES` for the shipped ones and why they differ."""

    #: Phred floor on the base that discriminates a clonotype from its abundance parent. 0 = off.
    #: The plateau is Q20-32 and it starts eating real variants by Q35 (MIGEC spike-ins).
    min_junction_q: int = 0
    #: A clonotype whose reads are worse than this (median of per-read mean junction Phred) is a
    #: candidate for the wide-radius rescue below. 0 = off. Chosen from the k-table above: the
    #: dominant clone and the call splits both sit at 40.5, the cliff class at 16.5-20.1.
    lowq_mean_q: float = 0.0
    #: Rescue radius, in substitutions, for those candidates ONLY. Never applied to a clonotype
    #: whose reads are fine -- that would be the abundance model reaching past its evidence.
    lowq_max_subs: int = 0
    #: A rescue parent must be at least this many times more abundant. The cliff class is 1-read
    #: clonotypes beside a 9,932-read clone, so a large ratio costs nothing there while protecting
    #: a genuine low-frequency variant that happens to sit near an abundant clone.
    lowq_min_ratio: float = 100.0

    def enabled(self) -> bool:
        return self.lowq_mean_q > 0 and self.lowq_max_subs > 0


#: Shipped regimes. ``fast`` and ``accurate`` are the historical ``--ec-mode`` values and keep
#: their exact behaviour; ``amplicon`` and ``rnaseq`` add the quality-directed rescue.
#:
#: ⛔ **The two regimes differ because their clonotype-size distributions differ, not by taste.**
#: An amplicon library is deep -- a real clonotype has many reads, so a 1-read neighbour of an
#: abundant clone is almost always error, and the rescue can be wide. Bulk RNA-seq is sparse:
#: singletons are the norm and most of them are real, the receptor fraction being 0.02-3 %. So
#: ``rnaseq`` keeps the rescue narrow and demands a much stronger abundance ratio, and ``amplicon``
#: opens it up. Both keep every read either way.
REGIMES: dict[str, DenoiseParams] = {
    "fast": DenoiseParams(),
    "accurate": DenoiseParams(min_junction_q=20),
    "amplicon": DenoiseParams(min_junction_q=20, lowq_mean_q=25.0, lowq_max_subs=12,
                              lowq_min_ratio=50.0),
    "rnaseq": DenoiseParams(min_junction_q=20, lowq_mean_q=20.0, lowq_max_subs=6,
                            lowq_min_ratio=200.0),
}


@dataclass
class DenoiseReport:
    """What the rescue did, in the terms the invariant is checked in."""

    lowq_clonotypes: int = 0        # clonotypes whose reads are below `lowq_mean_q`
    rescued_clonotypes: int = 0     # ...of which found a parent and were routed to it
    rescued_reads: int = 0          # reads those carried (MOVED to the parent, not dropped)
    orphan_clonotypes: int = 0      # low-quality but no parent found -- KEPT, deliberately
    orphan_reads: int = 0
    quality_available: bool = False
    dist: dict[int, int] = field(default_factory=dict)   # rescue distance -> clonotypes

    def as_dict(self) -> dict:
        return {**self.__dict__, "dist": dict(sorted(self.dist.items()))}


def read_quality(junctions: list[str], quals: list[str]) -> list[float]:
    """Per-read mean Phred over the junction; ``-1.0`` where there is no usable quality.

    ⛔ ``-1.0`` means ABSENT evidence and must never be read as bad evidence. A quality string of
    the right length taken from the wrong strand or offset is the one corruption nothing downstream
    can detect, so a length disagreement is refused here rather than averaged.
    """
    if _cpp is not None:
        return _cpp.mean_phred(junctions, quals)
    return _read_quality_py(junctions, quals)


def _read_quality_py(junctions: list[str], quals: list[str]) -> list[float]:
    """Reference implementation of :func:`read_quality` (asserted identical in the tests)."""
    if len(junctions) != len(quals):
        raise ValueError("junctions and quals must be the same length")
    out = []
    for jn, q in zip(junctions, quals):
        if not q or len(q) != len(jn):
            out.append(-1.0)
        else:
            out.append(sum(ord(c) - 33 for c in q) / len(q))
    return out


def clonotype_quality(keys: list, read_q: list[float]) -> dict:
    """Median of the per-read mean junction Phred, per clonotype key.

    The MEDIAN, not the mean: one catastrophic read must not drag an otherwise clean clonotype
    below the threshold, and a clonotype of genuinely bad reads is bad at its median too. Reads
    with no usable quality (``-1.0``) are excluded rather than counted as 0; a clonotype with no
    usable read at all is absent from the result and is therefore never a rescue candidate.
    """
    acc: dict = {}
    for k, q in zip(keys, read_q):
        if q >= 0:
            acc.setdefault(k, []).append(q)
    return {k: median(v) for k, v in acc.items()}


def quality_rescue(seqs: list[str], counts: list[int], clono_q: list[float],
                   params: DenoiseParams) -> tuple[list[int | None], DenoiseReport]:
    """Parent index for each clonotype the QUALITY evidence justifies collapsing, else ``None``.

    Only clonotypes whose reads are below ``params.lowq_mean_q`` are considered, and each may only
    join a neighbour at least ``params.lowq_min_ratio`` times more abundant within
    ``params.lowq_max_subs`` substitutions. A candidate with no such neighbour gets ``None`` and
    **keeps its reads** -- see the module docstring for why deleting it is not on the table.

    ``clono_q`` is parallel to ``seqs``; use ``-1.0`` for a clonotype with no usable quality, which
    excludes it (absent evidence, not bad evidence).

    ⛔ **V and J are deliberately IGNORED here, and that is a decision, not an omission.** The
    abundance model in ``correct._parents`` defaults to ``require_vj=True`` on the principle that a
    true sequencing error keeps the germline V/J call -- correct there, because it collapses
    1-3 substitution neighbours, which rarely move an alignment onto a different gene. This
    function targets the opposite class: the cliff, where the *whole junction window* is unreliable
    (median mean Phred 16.5-20.1, 6-12 substitutions). **A read that bad has an unreliable V/J call
    for exactly the same reason its junction is unreliable** -- the call came from aligning that
    sequence -- so requiring the calls to agree would filter on the corrupted evidence and defeat
    the rescue. Measured on SRR5233636 at full depth: of 9,025 rescues under ``amplicon``,
    **4,593 (50.9 %) cross the V call** and 908 the J; under ``rnaseq`` (6 subs, 200x) it is 24 and
    2 of 215. What protects a genuine clone here is not the call but the two gates that ARE
    trustworthy: its reads must be measurably bad (``lowq_mean_q``) and the parent must be
    ``lowq_min_ratio`` times more abundant.

    ⛔ The LOCUS is a different matter and is NOT ignored -- ``correct_airr`` partitions the search
    by it. A locus is fixed by the whole read (V, J and C genes together), not by junction bases, so
    a locus flip is not a plausible consequence of junction miscalls, and a rearrangement of another
    locus is not a sequencing error of this one. Measured before the partition: 3 of those 9,025
    were 1-read TRB clonotypes absorbed into abundant TRA clonotypes at 11-12 substitutions.
    """
    rep = DenoiseReport(quality_available=any(q >= 0 for q in clono_q))
    n = len(seqs)
    parent: list[int | None] = [None] * n
    if not params.enabled() or not rep.quality_available:
        return parent, rep

    cand = [0 <= clono_q[i] < params.lowq_mean_q for i in range(n)]
    rep.lowq_clonotypes = sum(cand)
    if not rep.lowq_clonotypes:
        return parent, rep

    if _cpp is not None:
        found = _cpp.nearest_more_abundant(seqs, [int(c) for c in counts], cand,
                                           params.lowq_max_subs, params.lowq_min_ratio)
    else:                                          # pragma: no cover - source checkout without ext
        found = _nearest_py(seqs, counts, cand, params.lowq_max_subs, params.lowq_min_ratio)

    for i in range(n):
        if not cand[i]:
            continue
        p = found[i]
        if p is None or p < 0:
            rep.orphan_clonotypes += 1
            rep.orphan_reads += counts[i]
            continue
        parent[i] = int(p)
        rep.rescued_clonotypes += 1
        rep.rescued_reads += counts[i]
        d = sum(1 for a, b in zip(seqs[i], seqs[int(p)]) if a != b)
        rep.dist[d] = rep.dist.get(d, 0) + 1
    return parent, rep


def _nearest_py(seqs, counts, cand, max_subs, min_ratio):
    """Reference implementation of the wide-radius search (asserted identical in the tests)."""
    by_len: dict[int, list[int]] = {}
    for i, s in enumerate(seqs):
        by_len.setdefault(len(s), []).append(i)
    out: list[int] = [-1] * len(seqs)
    for ci, want in enumerate(cand):
        if not want or not seqs[ci]:
            continue
        need = counts[ci] * min_ratio
        best, best_count = -1, -1
        for nj in by_len.get(len(seqs[ci]), ()):
            if nj == ci or counts[nj] <= counts[ci] or counts[nj] < need:
                continue
            mm = 0
            for a, b in zip(seqs[ci], seqs[nj]):
                if a != b:
                    mm += 1
                    if mm > max_subs:
                        break
            if mm == 0 or mm > max_subs:
                continue
            if counts[nj] > best_count or (counts[nj] == best_count and seqs[nj] < seqs[best]):
                best, best_count = nj, counts[nj]
        out[ci] = best
    return out
