"""A junction whose read stopped before [FW]118, finished from the called J's germline.

The 3' boundary is the only one arda can rebuild without observing it: the J's 5' chew-back and
the N/P additions all sit UPSTREAM of the read's last aligned J base, so what is missing from there
to [FW]118 is J-germline templated. The V side has no counterpart, which is why
``v_anchor_prefix`` refuses rather than extrapolates.

Everything here is about what the completion must REFUSE. Measured on a 100 k-read TRA amplicon:
of 266 reads that reached Cys104 without reaching [FW]118, **236 carry more than a partial codon of
unaligned 3' sequence, and that sequence is TRAC** -- V straight into the constant region with no J
at all (the aligner still names a J off a few coincidental bases). Completing those would have
manufactured 236 junctions for reads that carry no rearrangement to Cys104 and beyond.
"""
from __future__ import annotations

from arda.annotate.transfer import _MAX_UNALIGNED_TAIL_NT, _germline_completed_junction
from arda.cdr3fix import Anchor


# TRAJ-shaped: `germline_nt` runs the J's 5' end through the [FW]118 codon's last base, so the
# anchor codon is the final three nt -- the invariant the completion grafts onto.
J_GERMLINE = "AACTATGGTCAGAATTTTGTCTTT"          # ...N Y G Q N F V F
ANCHORS = {("J", "TRAJ26*01"): Anchor(
    locus="TRA", segment="J", templated_aa="NYGQNFVF", functionality="F", status="ok",
    source="aux", anchor_nt=len(J_GERMLINE) - 3, germline_nt=J_GERMLINE)}

# A read: V framework .. Cys104 .. N/P .. the first 10 nt of the J, and then it stops.
V_TAIL = "GCTGTGTACTATTGC"                       # ends on the Cys104 codon (TGC)
NP = "ATCGTCAGAGTC"
READ = V_TAIL + NP + J_GERMLINE[:10]
CS = len(V_TAIL) + 1                             # CDR3 starts just after the Cys codon (1-based)


def _complete(read=READ, *, j_germline_end=10, q_aln_end=None, max_nt=40, j_call="TRAJ26*01"):
    return _germline_completed_junction(
        read, CS, coding_start=1, v_end_q=len(V_TAIL),
        j_call=j_call, j_germline_end=j_germline_end,
        q_aln_end=len(read) if q_aln_end is None else q_aln_end,
        anchors=ANCHORS, max_nt=max_nt)


def test_it_finishes_the_junction_at_the_anchor():
    """The completed junction is the observed bases plus exactly the untranscribed germline tail,
    and it closes on the [FW]118 codon -- which is what made it a junction rather than a fragment."""
    got = _complete()
    assert got is not None, "a read 14 nt short of [FW]118 was not completed"
    jnt, jaa, _cdr3, _phase, n_done = got

    assert n_done == len(J_GERMLINE) - 10, "wrong number of germline bases imputed"
    assert jnt == READ[CS - 4:] + J_GERMLINE[10:], "the graft is not observed-bases + germline-tail"
    assert jnt.startswith("TGC"), "the junction must still open on the V's Cys104 codon"
    assert jaa.endswith("F"), f"junction {jaa!r} does not close on [FW]118"


def test_it_refuses_when_the_read_carries_UNALIGNED_sequence_past_the_alignment():
    """⛔ The one that matters. A read whose alignment stops well short of its own 3' end stopped
    for a reason, and on real data that reason is usually that there is no J: 236 of 266 candidates
    on a TRA amplicon run straight from the V into TRAC, so the bases after the alignment are
    CONSTANT REGION. Grafting J germline over them invents a junction for a V-to-C chimera.
    """
    chimeric = READ[:len(V_TAIL) + len(NP)] + "ATCCAGAACCCTGACCCCCTTGCTTT"   # TRAC, not J
    aln_end = len(V_TAIL) + len(NP)

    assert _complete(chimeric, j_germline_end=6, q_aln_end=aln_end) is None, (
        "a V-to-C chimera was completed into a junction it has no J for"
    )
    # The boundary itself: a final partial codon that an ungapped extension would not carry is
    # still completed, one base more is not.
    assert _complete(READ + "A" * _MAX_UNALIGNED_TAIL_NT, q_aln_end=len(READ)) is not None
    assert _complete(READ + "A" * (_MAX_UNALIGNED_TAIL_NT + 1), q_aln_end=len(READ)) is None


def test_it_refuses_to_impute_more_than_asked():
    """`--complete-junctions N` is a budget, not a switch: the risk scales with how much is
    imputed, so a read needing more than N nt is declined outright rather than partly filled."""
    need = len(J_GERMLINE) - 10
    assert _complete(max_nt=need) is not None
    assert _complete(max_nt=need - 1) is None
    assert _complete(max_nt=0) is None, "max_nt=0 is the shipped default and must complete nothing"


def test_it_refuses_without_a_usable_J_anchor():
    """An allele with no anchor cannot say where [FW]118 is, so there is no tail to graft. Left
    unscoped rather than guessed -- the same rule the SHM scoping follows."""
    assert _complete(j_call="TRAJ99*01") is None, "completed against an allele not in the anchors"
    assert _complete(j_call="") is None

    no_anchor = {("J", "TRAJ26*01"): Anchor(
        locus="TRA", segment="J", templated_aa="", functionality="F", status="no_anchor",
        source="no_anchor", anchor_nt=-1, germline_nt="")}
    assert _germline_completed_junction(
        READ, CS, coding_start=1, v_end_q=len(V_TAIL), j_call="TRAJ26*01",
        j_germline_end=10, q_aln_end=len(READ), anchors=no_anchor, max_nt=40) is None


def test_a_read_that_already_reached_the_anchor_has_nothing_to_complete():
    """No tail left means no completion -- and, in the caller, no `junction_completed_nt`. A row
    that says it imputed 0 nt would be indistinguishable from one that imputed nothing."""
    assert _complete(j_germline_end=len(J_GERMLINE)) is None
    assert _complete(j_germline_end=len(J_GERMLINE) + 5) is None
