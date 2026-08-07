"""The junction placed by arithmetic: the coordinate contract, and the refusals.

`project_junction` composes three coordinate systems that disagree about their origin -- `tstart` is
1-based on the forward target, `anchor_nt` is 0-based in the germline, and a minus-strand `qstart`
is in forward coordinates while the sequence it indexes is the reverse complement. Every one of
those is an off-by-one waiting to happen, and an off-by-one here still produces a junction that
looks plausible: right length, starts with a codon, ends with a codon. So the tests pin the
arithmetic against a hand-built case where the answer is known by construction, and then against the
SHIPPED reference, where the anchor codon must actually decode to Cys / Phe / Trp.
"""

from __future__ import annotations

import pytest

from arda.annotate.project import project_junction
from arda.cdr3fix import Anchor

from ..conftest import requires_human_db

COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def anchor(nt: int, status: str = "ok") -> Anchor:
    return Anchor(locus="TRB", segment="V", templated_aa="C", functionality="F", status=status,
                  source="test", anchor_nt=nt, partial_nt=0, germline_nt="TGT")


# A read built so the answer is known: 30 nt of V, then TGT (Cys104), then a 12 nt NDN, then TTT
# (Phe118), then 15 nt of FR4. The junction is therefore TGT + NDN + TTT = 18 nt.
V_LEAD, NDN, FR4 = "A" * 30, "GGGCCCAAATTT", "C" * 15
READ = V_LEAD + "TGT" + NDN + "TTT" + FR4
JUNCTION = "TGT" + NDN + "TTT"


@pytest.fixture
def scene():
    """V hit starting at read base 1 against germline base 1; J hit starting at the [FW] codon."""
    v_row = {"qstart": 1, "qend": 33, "tstart": 1, "split": False}
    j_row = {"qstart": 46, "qend": 61, "tstart": 1, "split": False}
    # V germline: the Cys codon begins at read offset 30 (0-based), and the hit is germline-aligned
    # from base 1, so anchor_nt is 30. J germline: the hit starts AT the [FW] codon, so anchor_nt 0.
    anchors = {("V", "TRBV1*01"): anchor(30), ("J", "TRBJ1*01"): anchor(0)}
    return dict(v_row=v_row, j_row=j_row, v_call="TRBV1*01", j_call="TRBJ1*01", anchors=anchors)


def test_the_projected_junction_is_the_one_that_is_actually_there(scene):
    p, why = project_junction(READ, len(READ), split_checked=True, **scene)
    assert why == ""
    assert p.junction == JUNCTION
    assert READ[p.start - 1:p.end] == JUNCTION


def test_it_returns_the_AIRR_junction_not_the_IMGT_cdr3(scene):
    # These differ by exactly the two anchor codons, and confusing them corrupts every downstream
    # comparison silently -- junction_aa is two residues longer than cdr3_aa.
    p, _ = project_junction(READ, len(READ), split_checked=True, **scene)
    assert p.junction.startswith("TGT") and p.junction.endswith("TTT")
    assert p.cdr3 == NDN
    assert len(p.junction) == len(p.cdr3) + 6


def test_a_reverse_complement_hit_lands_on_the_same_junction(scene):
    # This is the branch that matters most in practice: R2 of a paired amplicon is 100 % rc, so a
    # sign error here would corrupt half of every paired run while the forward half looked perfect.
    #
    # The INPUT read is rc(READ); segmap therefore finds V and J while scanning its reverse
    # complement, which is READ. `segmap.cpp:383` converts those back to FORWARD coordinates on the
    # input read, giving qstart > qend. The caller passes the strand the hits were MEASURED on --
    # revcomp(input) == READ -- because that is the frame the offsets are expressed in.
    n = len(READ)
    # V occupies READ 0-based [0, 32] -> qstart = n - 0 = 63, qend = n - 32 = 31.
    v_row = {"qstart": n - 0, "qend": n - 32, "tstart": 1, "split": False}
    # The [FW] codon begins at READ 0-based 45; the J hit runs to the end of FR4 at 0-based 62.
    j_row = {"qstart": n - 45, "qend": n - 62, "tstart": 1, "split": False}
    p, why = project_junction(READ, n, split_checked=True, **{**scene, "v_row": v_row, "j_row": j_row})
    assert why == ""
    assert p.rev_comp is True
    assert p.junction == JUNCTION
    assert READ[p.start - 1:p.end] == JUNCTION


def test_a_hit_starting_mid_germline_still_projects(scene):
    # A read that does not reach the start of V: the hit begins at germline base 11, so the anchor
    # is 10 bases closer. If `tstart` were ignored the junction would be off by exactly 10.
    trimmed = READ[10:]
    v_row = {"qstart": 1, "qend": 23, "tstart": 11, "split": False}
    j_row = {"qstart": 36, "qend": 51, "tstart": 1, "split": False}   # shifted with the read
    p, why = project_junction(trimmed, len(trimmed), split_checked=True,
                              **{**scene, "v_row": v_row, "j_row": j_row})
    assert why == ""
    assert p.junction == JUNCTION


@pytest.mark.parametrize("mutate, expected", [
    (lambda s: {**s, "anchors": {("V", "TRBV1*01"): anchor(30, status="no_anchor"),
                                 ("J", "TRBJ1*01"): anchor(0)}}, "no_anchor"),
    (lambda s: {**s, "v_row": {**s["v_row"], "split": True}}, "indel_split"),
    (lambda s: {**s, "j_row": {"qstart": 61, "qend": 46, "tstart": 1, "split": False}},
     "strand_mismatch"),
    (lambda s: {**s, "j_row": {**s["j_row"], "tstart": 400}}, "order"),
])
def test_every_refusal_returns_no_junction_and_names_itself(scene, mutate, expected):
    p, why = project_junction(READ, len(READ), split_checked=True, **mutate(scene))
    assert p is None and why == expected


def test_an_anchor_projected_past_the_end_of_the_read_is_refused(scene):
    # The read stops before FR4, so the J anchor projects off the end. Returning a truncated
    # junction here is exactly the "well-formed but systematically short" failure that shipped once.
    short = READ[:40]        # junction ends at 48, so the J anchor projects past the end
    p, why = project_junction(short, len(short), split_checked=True, **scene)
    assert p is None and why == "off_read"


def test_the_composite_allele_names_a_segment_target_can_carry_resolve(scene):
    # 23 of 775 human V targets are comma-joined duplicate sequences. A bare dict lookup on the
    # composite string misses them, which would refuse a read that is perfectly placeable.
    p, why = project_junction(READ, len(READ), split_checked=True,
                              **{**scene, "v_call": "TRBV1*01,TRBV1D*01"})
    assert why == "" and p.junction == JUNCTION


@requires_human_db
def test_the_anchor_offsets_decode_to_real_anchor_codons_in_the_shipped_reference():
    """The arithmetic is only meaningful if `anchor_nt` indexes the segment targets it is used with.

    ⛔ This is the assertion that catches a reference rebuild silently moving the coordinate system
    out from under the projection -- the class of failure this repo has shipped as "a reference swap
    can silently be a no-op". It reads the REAL segments.fasta, not a fixture.
    """
    from arda.annotate.reference import load_reference
    from arda.refbuild.segments import build_segment_reference

    ref = load_reference("human", "nt")
    seg = ref.target_fasta.parent / "segments.fasta"
    if not seg.exists():
        build_segment_reference("human", out_dir=ref.target_fasta.parent)

    seqs, name = {}, None
    for line in seg.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()

    ok = {"V": 0, "J": 0}
    bad = {"V": [], "J": []}
    for target, s in seqs.items():
        kind, _, call = target.partition("|")
        if kind not in ("V", "J"):
            continue
        for allele in call.split(","):
            a = ref.anchors.get((kind, allele.strip()))
            if a is None or a.status != "ok" or a.anchor_nt < 0:
                continue
            codon = s[a.anchor_nt:a.anchor_nt + 3].upper()
            want = {"TGT", "TGC"} if kind == "V" else {"TTT", "TTC", "TGG"}
            if codon in want:
                ok[kind] += 1
            else:
                bad[kind].append((target, codon))
            break

    # Not 100 %: TRBJ2-2P*01 (GGG) and TRBJ2-7*02 (GTC) are genuine germline exceptions, and a
    # handful of V alleles carry a non-canonical anchor. The bar is that the coordinate system is
    # right, which a systematic off-by-one would blow straight through.
    assert ok["V"] >= 700, f"V anchors decoding to Cys: {ok['V']}, misses {bad['V'][:5]}"
    assert ok["J"] >= 115, f"J anchors decoding to F/W: {ok['J']}, misses {bad['J'][:5]}"


def test_an_unvalidated_locus_is_declined_rather_than_guessed(scene):
    """TRD is refused because it has ZERO coverage, not because it is known wrong.

    The pre-registered bar was "byte-exact >= 0.99 at n >= 2,000 per locus, or the locus goes on
    the refusal list". Across two TR amplicons the segment pass never handed a single TRD read both
    anchors, so all 767 TRD junctions in the IgBLAST truth went to the aligner untouched -- n = 0,
    which is not a pass. Declining costs nothing today and stops a future reference change from
    routing TRD through arithmetic nobody has measured.
    """
    trd = Anchor(locus="TRD", segment="J", templated_aa="F", functionality="F", status="ok",
                 source="test", anchor_nt=0, partial_nt=0, germline_nt="TTT")
    anchors = {("V", "TRBV1*01"): anchor(30), ("J", "TRBJ1*01"): trd}
    p, why = project_junction(READ, len(READ), split_checked=True, **{**scene, "anchors": anchors})
    assert p is None and why == "unvalidated_locus"


def test_the_locus_is_read_off_the_J_anchor_not_the_V(scene):
    # TRAV/DV alleles rearrange to either TRAJ or TRDJ and the J decides the locus. A V-derived
    # locus would let a TRD read through under a TRA label -- the exact asymmetry refbuild/loci.py
    # encodes and that this project has already built, measured and reverted once.
    v_trd = Anchor(locus="TRD", segment="V", templated_aa="C", functionality="F", status="ok",
                   source="test", anchor_nt=30, partial_nt=0, germline_nt="TGT")
    anchors = {("V", "TRBV1*01"): v_trd, ("J", "TRBJ1*01"): anchor(0)}
    p, why = project_junction(READ, len(READ), split_checked=True, **{**scene, "anchors": anchors})
    assert why == "" and p.junction == JUNCTION


def test_an_unchecked_indel_gate_refuses_instead_of_passing_silently(scene):
    """`split` is 0 both when there is no indel and when nobody looked. Those must not be the same.

    `segment_rows` only populates `split` when `max_indel > 0`, i.e. when `--indel-rescue` is on.
    Without it every `split` is 0 and a projection would run with its indel protection INERT --
    this codebase's most repeated failure shape, and the reason `split_checked` has no default.

    It is worth a refusal because the gate is worth something: on IGH_repertoire (91.77 % median V
    identity) it takes byte-exact accuracy from .99634 to .99915 -- 332 wrong junctions down to 74 --
    for 3.97 % of fast-path yield.
    """
    p, why = project_junction(READ, len(READ), split_checked=False, **scene)
    assert p is None and why == "indel_unchecked"
