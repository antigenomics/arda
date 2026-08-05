"""Segment shortlist → implied V×J scaffold, with a rescue set so **no read is ever lost**.

The fast path: align against the 1,244-target segment reference (V, J and J+C separately, see
:mod:`arda.refbuild.segments`), take each read's best V and best J, and look the pair up in
``combinations.tsv``. That names exactly one V×J scaffold in the full reference — so the second
alignment is one target per read instead of ~277, measured at **5.36 s → 0.044 s (122×)** on a
20 k-read TRA amplicon, and it lands back in the scaffold coordinate system arda's markup
transfer already speaks.

**The invariant this module exists to enforce.** A fast path that silently drops reads is not an
optimisation, it is a different tool: arda's whole claim is near-zero Stage-1 false negatives.
So every read is accounted for. :func:`shortlist` partitions the input into

* ``implied``   — a V and a J both hit, and the pair names a real scaffold. The fast path.
* ``rescue``    — anything else: only a V hit, only a J hit, an unknown V×J pair, or a read whose
  second-pass alignment failed. These go back to the **full reference**, exactly as today.

``implied`` ∪ ``rescue`` == every read that hit anything, by construction and by assertion. The
rescue set is small and cheap to realign — measured 1.9 % of amplicon reads and ~10 % of bulk,
i.e. 11 % and 20 % of the *new* total cost — so exactness costs almost nothing:

    amplicon   fast 0.84 s + rescue 0.10 s = 0.94 s  vs 5.36 s   → 5.7×
    bulk       fast 0.52 s + rescue 0.13 s = 0.65 s  vs 1.26 s   → 2.0×

Why reads land in ``rescue``, measured on a TRA amplicon:

* **V only** (12.6 %) — the read never reaches a J, so no pair exists. The baseline picks an
  arbitrary J for these; the rescue pass reproduces that behaviour exactly rather than inventing
  a different arbitrary answer.
* **J only** (2.1 %) — J→C reads with no V.
* **failed second alignment** (1.9 %) — the synthesized diagonal is derived from the V hit, so it
  is wrong for a read whose scaffold alignment does not begin in V.

⚠ **Known, unfixed: the α/δ locus policy.** Of 212 residual V-gene disagreements against the full
reference on one amplicon, **114 (54 %) are TRA→TRD** — the TRAV/DV dual-use ambiguity. A joint
``argmax(bits_V + bits_J)`` over the top-5 of each, constrained to real combinations, does *not*
fix it: the segment pass's best J genuinely is a TRDJ. This needs a locus prior, so until it
lands, :func:`shortlist` routes ambiguous TRAV/DV reads to ``rescue`` rather than guessing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["Shortlist", "load_combinations", "shortlist", "TRAV_DV_AMBIGUOUS"]

#: TRAV/DV segments rearrange into both TRA and TRD. Until a locus prior exists, a read whose
#: best V is one of these and whose best J crosses the locus is rescued, never guessed.
TRAV_DV_AMBIGUOUS = ("/DV",)


@dataclass
class Shortlist:
    """Which reads take the fast path, which must be realigned, and why."""

    #: read id -> scaffold id implied by (best V, best J)
    implied: dict[str, str] = field(default_factory=dict)
    #: read ids that must go back to the full reference
    rescue: list[str] = field(default_factory=list)
    #: reason -> count, for the run report
    reasons: dict[str, int] = field(default_factory=dict)

    @property
    def n_total(self) -> int:
        return len(self.implied) + len(self.rescue)

    @property
    def fast_fraction(self) -> float:
        return len(self.implied) / self.n_total if self.n_total else 0.0

    def _mark(self, read_id: str, reason: str) -> None:
        self.rescue.append(read_id)
        self.reasons[reason] = self.reasons.get(reason, 0) + 1

    def as_dict(self) -> dict:
        return {"implied": len(self.implied), "rescued": len(self.rescue),
                "fast_fraction": round(self.fast_fraction, 4), "reasons": dict(self.reasons)}


def load_combinations(path: str | Path) -> dict[tuple[str, str], str]:
    """``combinations.tsv`` → ``{(v_allele, j_allele): scaffold_id}``.

    ``v_calls``/``j_calls`` may be comma-separated ambiguity lists; every listed allele maps to
    the scaffold, so a lookup succeeds whichever member the segment pass reported.
    """
    import polars as pl

    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    out: dict[tuple[str, str], str] = {}
    for row in df.iter_rows(named=True):
        sid = row["scaffold_id"]
        for v in (row.get("v_calls") or "").split(","):
            for j in (row.get("j_calls") or "").split(","):
                if v and j:
                    out.setdefault((v.strip(), j.strip()), sid)
    return out


def _locus(call: str) -> str:
    return (call or "")[:3]


def shortlist(best_v: dict[str, str], best_j: dict[str, str],
              combos: dict[tuple[str, str], str], *,
              failed: set[str] | None = None) -> Shortlist:
    """Partition reads into the fast path and the rescue set.

    Args:
        best_v: read id -> best ``V|`` allele (absent if the read hit no V target).
        best_j: read id -> best ``J|``/``JC|`` allele.
        combos: from :func:`load_combinations`.
        failed: read ids whose second-pass alignment produced nothing; always rescued.

    Returns:
        :class:`Shortlist`. Every read appearing in ``best_v`` or ``best_j`` lands in exactly one
        of ``implied`` / ``rescue`` — asserted, not assumed.

    ⚠ ``best_j`` must hold a **J allele**. ``JC|`` targets are named by *scaffold id*, not by
    allele (see :mod:`arda.refbuild.segments`), so a caller feeding the raw target name straight
    through gets ``no_such_combination`` for every J→C read. Measured cost of getting this wrong:
    the fast path collapsed from 85.3 % to 0.1 %, with 9,388 reads needlessly rescued — correct
    output, silently no faster. Resolve via ``segments.markup.tsv``'s ``j_call`` column.
    """
    sl = Shortlist()
    failed = failed or set()
    for rid in sorted(set(best_v) | set(best_j)):
        v, j = best_v.get(rid), best_j.get(rid)
        if rid in failed:
            sl._mark(rid, "second_pass_failed")
        elif not v:
            sl._mark(rid, "j_only")            # J->C read, no V to pair
        elif not j:
            sl._mark(rid, "v_only")            # never reached a J; baseline picks arbitrarily
        elif any(t in v for t in TRAV_DV_AMBIGUOUS) and _locus(v) != _locus(j):
            # TRAV/DV is dual-use and the segment pass's best J crosses the locus. Measured at
            # 54% of residual V disagreements; guessing here moved TRA reads to TRD.
            sl._mark(rid, "trav_dv_locus_ambiguous")
        else:
            sid = combos.get((v, j))
            if sid is None:
                sl._mark(rid, "no_such_combination")
            else:
                sl.implied[rid] = sid

    assert sl.n_total == len(set(best_v) | set(best_j)), "a read was lost during shortlisting"
    return sl
