"""Two AIRR fields were being reported without evidence, and both are expressible offline.

1. ``junction`` — a rearrangement can trim V back **past Cys104**. The scaffold projection still
   lands somewhere and emits a junction opening on bases the V germline never templated. On a TRA
   amplicon (results/round18) 1,396 of 46,785 reads disagreed with IgBLAST and **every single one
   was a pure 5' over-extension** — the 3' [FW]118 end was right in 100 % of them. Requiring the
   called V's own junction germline to explain the opening codon kills 1,360 of those for 92
   correct junctions.

2. ``j_call`` — it was copied from the scaffold unconditionally, so a V-only read inherited the J
   of whichever V×J scaffold it landed on. On bulk RNA-seq that was 1,823 of 2,737 mapped reads,
   and j_call precision against IgBLAST was .1129.

Both are expressed here as a hand-built ``RefEntry`` + hit dict: no reference DB, no MMseqs2,
nothing to skip. CLAUDE.md: gating a regression test on ``requires_human_db`` is how the last
reference defect regressed.
"""

from __future__ import annotations

from types import SimpleNamespace

from arda.annotate.reference import REGIONS, RefEntry
from arda.annotate.transfer import (
    _ANCHOR_MIN_MATCH,
    _ANCHOR_WINDOW,
    MIN_V_ANCHOR_PREFIX,
    transfer_hit,
    v_anchor_match,
    v_anchor_ok,
    v_anchor_prefix,
)

# TRAV25*01's templated junction germline: Cys104 codon + AGG GGG. Its Cys sits at germline 264,
# and the read that produced the +9 nt cluster had V trimmed to 258 — six nt short of the anchor.
TRAV25_GERMLINE_NT = "TGTGCAGGG"


def _anchors(germline_nt=TRAV25_GERMLINE_NT):
    return {("V", "TRAV25*01"): SimpleNamespace(status="ok", germline_nt=germline_nt)}


def test_a_junction_trimmed_past_cys104_is_not_explained_by_its_own_v():
    """The real failing read: arda's window opened 9 nt early, on bases V never templated."""
    over_extended = "ACCATGAACCAGGGAGGAAAGCTTATCTTC"      # what arda emitted
    assert v_anchor_prefix(over_extended, "TRAV25*01", _anchors()) < MIN_V_ANCHOR_PREFIX


def test_a_correctly_anchored_junction_is_explained_to_the_end_of_the_templated_stretch():
    assert v_anchor_prefix("TGTGCAGGGAAAGCTTATCTTC", "TRAV25*01",
                           _anchors()) == len(TRAV25_GERMLINE_NT)


def test_the_cut_survives_the_one_shm_event_that_hits_the_conserved_codon():
    """A synonymous TGT->TGC leaves two bases, which is exactly why the cut is 2 and not 3."""
    assert MIN_V_ANCHOR_PREFIX == 2
    assert v_anchor_prefix("TGCGCAGGGAAAGCTTATCTTC", "TRAV25*01",
                           _anchors()) >= MIN_V_ANCHOR_PREFIX


def test_an_ambiguous_v_call_is_explained_by_whichever_allele_explains_most():
    anchors = _anchors()
    anchors[("V", "TRAV25*04")] = SimpleNamespace(status="ok", germline_nt="TGTGCTGGG")
    assert v_anchor_prefix("TGTGCTGGGAAAG", "TRAV25*01,TRAV25*04", anchors) == 9


def test_the_over_extension_stays_refused_by_the_mismatch_tolerant_rule_too():
    """Loosening the cut must not readmit what it exists to catch: 0 of 8 bases match here."""
    over_extended = "ACCATGAACCAGGGAGGAAAGCTTATCTTC"
    assert v_anchor_match(over_extended, "TRAV25*01", _anchors()) == 0
    assert not v_anchor_ok(over_extended, "TRAV25*01", _anchors())


def test_a_substitution_in_the_first_base_no_longer_throws_the_junction_away():
    """One mismatch takes a longest common prefix to 0; the window still sees seven of eight.

    That is the whole defect: an exact prefix is one sequencing error from zero, and on the
    held-out TRB amplicon it was discarding 116 correct junctions to catch 21 wrong ones
    (benchmark round 28).
    """
    mutated = "AGTGCAGGGAAAGCTTATCTTC"               # TRAV25*01's Cys104 with base 1 substituted
    assert v_anchor_prefix(mutated, "TRAV25*01", _anchors()) < MIN_V_ANCHOR_PREFIX
    assert v_anchor_match(mutated, "TRAV25*01", _anchors()) == _ANCHOR_WINDOW - 1
    assert v_anchor_ok(mutated, "TRAV25*01", _anchors())


def test_the_window_is_a_rescue_and_never_a_second_chance_to_refuse():
    """An exact prefix admits regardless of what the rest of the window does."""
    assert _ANCHOR_MIN_MATCH <= _ANCHOR_WINDOW
    short = {("V", "TRAV25*01"): SimpleNamespace(status="ok", germline_nt="TGT")}
    assert v_anchor_match("TGTAAAAAA", "TRAV25*01", short) == 3 < _ANCHOR_MIN_MATCH
    assert v_anchor_ok("TGTAAAAAA", "TRAV25*01", short)


def test_an_unusable_anchor_row_never_licenses_a_junction():
    """``status != ok`` means the anchor was not established; it must not read as evidence."""
    bad = {("V", "TRAV25*01"): SimpleNamespace(status="no_anchor", germline_nt=TRAV25_GERMLINE_NT)}
    assert v_anchor_prefix("TGTGCAGGGAAAG", "TRAV25*01", bad) == 0


# --- j_call evidence gate ------------------------------------------------------------------

#: A minimal `V + 9 nt N-pad + J` scaffold: V germline 1..30, pad 31..39, J germline 40..69.
_V_END, _J_START, _VJ_END = 30, 40, 69


def _ref():
    return RefEntry(locus="TRB", v_call="TRBV1*01", j_call="TRBJ1-1*01",
                    starts=[-1] * len(REGIONS), ends=[-1] * len(REGIONS),
                    v_sequence_end=_V_END, j_sequence_start=_J_START, vj_end=_VJ_END)


def _hit(tstart: int, tend: int) -> dict:
    n = tend - tstart + 1
    aln = "A" * n
    return {"qaln": aln, "taln": aln, "qstart": 1, "qend": n, "qlen": n,
            "tstart": tstart, "tend": tend, "tlen": _VJ_END}


def test_a_read_that_never_reaches_j_germline_gets_no_j_call():
    rec = transfer_hit("r1", "A" * 25, _hit(1, 25), _ref())
    assert rec["v_call"] == "TRBV1*01"          # the V is real: the read is inside V germline
    assert rec["j_germline_start"] == ""
    assert rec["j_call"] == "", "the scaffold's J is not evidence that this read carries a J"


def test_a_read_that_reaches_j_germline_keeps_its_j_call():
    rec = transfer_hit("r2", "A" * _VJ_END, _hit(1, _VJ_END), _ref())
    assert rec["j_germline_start"] == 1
    assert rec["j_call"] == "TRBJ1-1*01"
