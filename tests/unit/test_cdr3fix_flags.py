"""The per-side flag vocabulary, and the four defects of issue #141.

Every case here is a real `(junction, V, J)` from VDJdb's curation corpus, and every expected
answer is the one the authoritative 2026-06-03 release ships -- so these are regressions against
measured behaviour, not against taste.

The corpus-scale numbers each case stands for are in `scripts/audit_cdr3fix.py`'s output and the
CHANGELOG; what is pinned here is the decision, per record.
"""
from __future__ import annotations

import pytest

from arda.cdr3fix import FLAGS, markup_cdr3

HS = "HomoSapiens"


def test_the_flag_vocabulary_is_closed():
    """Every flag any record can carry is one of the seven names, on both sides."""
    assert set(FLAGS) == {"ok", "allele", "sub", "add", "trim", "mismatch", "shallow",
                          "impossible"}
    for cdr3, v, j in (("CASSARSGELFF", "TRBV9*01", "TRBJ2-2*01"),
                       ("YFCASSQSPGGVAFFGQG", "TRBV14", "TRBJ1-1"),
                       ("CAISGEFGSGA", "TRBV10-1", "TRBJ2-6"),
                       ("CASSARSGELFF", "TRBV6", "TRBJ2-2"),
                       ("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01")):
        m = markup_cdr3(cdr3, v, j, HS)
        assert set(m.v_flags) <= set(FLAGS), (cdr3, m.v_flags)
        assert set(m.j_flags) <= set(FLAGS), (cdr3, m.j_flags)
        assert m.v_flags and m.j_flags, "a side always carries at least one flag"


def test_ok_and_impossible_are_never_on_the_same_side():
    """`ok` means nothing needed doing. It cannot coexist with a reason something does."""
    for cdr3, v, j in (("CASSARSGELFF", "TRBV9*01", "TRBJ2-2*01"),
                       ("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01"),
                       ("CAMREGDSSYKLIF", "TRAV14", "TRAJ12")):
        m = markup_cdr3(cdr3, v, j, HS)
        for flags in (m.v_flags, m.j_flags):
            assert not ("ok" in flags and len(flags) > 1), flags
            assert not ("ok" in flags and "shallow" in flags), flags


def test_good_means_neither_side_is_impossible():
    m = markup_cdr3("CASSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert m.v_flags == m.j_flags == ("ok",) and m.good
    bad = markup_cdr3("CALRPA", "TRAV9-2", "TRAJ17", HS)     # no J alignment at all
    assert "impossible" in bad.j_flags and not bad.good


def test_a_disagreement_inside_the_run_is_flagged_and_never_changed():
    """The policy in one record: not certain, so flag it -- do not change it.

    `CASSLGGNNKLFF` on TRBJ1-4 (`NEKLFF`) disagrees four residues in from Phe118. A curation error
    and an allele IMGT does not record look identical from one junction, so the submission stands
    and `mismatch` says why a curator should look. Rewriting these collapsed distinct records onto
    one string: `CAAAETSYDKV{M,R,T,V}F` on TRAJ50 are four junctions that all became
    `CAAAETSYDKVIF`, which is four residues at one templated position -- a suspect call or an
    unrecorded allele, not four independent typos.

    `mismatch` is NOT a refusal: the junction is well formed, every residue is the curator's, and
    the record stays `good`. `v_fix`/`j_fix` keep legacy's name because legacy's substring scanner
    never saw inside the run either, so the `cdr3fix` JSON stays comparable.
    """
    m = markup_cdr3("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01", HS)
    assert m.cdr3_repaired == m.cdr3
    assert m.j_flags == ("mismatch",) and m.j_fix == "NoFixNeeded" and m.good
    (e,) = [e for e in m.errors if e.side == "J"]
    assert (e.kind, e.dist, e.applied) == ("sub", 4, False)
    # and no setting of `max_replace` buys it: the knob bounds the anchor, not the run.
    for rung in (0, 1, 2, 4):
        assert markup_cdr3("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01", HS,
                           max_replace=rung).cdr3_repaired == "CASSLGGNNKLFF"


@pytest.mark.parametrize("cdr3, v, j", [
    ("CAAAETSYDKVMF", "TRAV13-1*01", "TRAJ50*01"),
    ("CAAAETSYDKVRF", "TRAV13-1*01", "TRAJ50*01"),
    ("CAAAETSYDKVTF", "TRAV13-1*01", "TRAJ50*01"),
    ("CAAAETSYDKVVF", "TRAV13-1*01", "TRAJ50*01"),
])
def test_four_residues_at_one_templated_position_stay_four_junctions(cdr3, v, j):
    """The record class that made the rule. TRAJ50 templates I where these carry M/R/T/V."""
    m = markup_cdr3(cdr3, v, j, HS)
    assert m.cdr3_repaired == cdr3
    assert "mismatch" in m.j_flags and m.good


# --------------------------------------------------------------------------- defect 1
@pytest.mark.parametrize("cdr3, v, j, want", [
    # Same defect at both ends of two junctions. Before the fix the DECISION FLIPPED per end,
    # because `aligned` was read off the alignment score and a 3-residue flank costs exactly the
    # 3 germline matches it unlocks: TRBJ1-1 (`NTEAFF`) netted 0 and declined, TRBJ1-2 (`NYGYTF`)
    # netted +1 and repaired. Both of these are what the 2026-06-03 release ships.
    ("YFCASSQSPGGVAFFGQG", "TRBV14", "TRBJ1-1", "CASSQSPGGVAFF"),
    ("YFCAVVGTGLGYTFGSG", "TRBV9", "TRBJ1-2", "CAVVGTGLGYTF"),
])
def test_framework_past_an_anchor_is_trimmed_at_both_ends(cdr3, v, j, want):
    m = markup_cdr3(cdr3, v, j, HS)
    assert m.cdr3_repaired == want
    assert m.v_flags == ("trim",) and m.j_flags == ("trim",)
    assert m.good


def test_an_alignment_is_not_read_off_the_score():
    """The J side of `YFCASSQSPGGVAFFGQG` scores exactly 0 and is still an alignment.

    3 exact germline matches (`AFF`) minus a 3-residue flank at `_TRIM` each. Gating on
    `score > 0` called that "no alignment" while the same defect with 4 matches repaired.
    """
    m = markup_cdr3("YFCASSQSPGGVAFFGQG", "TRBV14", "TRBJ1-1", HS)
    assert m.j_fix == "FixTrim" and "impossible" not in m.j_flags


# --------------------------------------------------------------------------- defect 2
def test_a_junction_cut_inside_the_germline_is_completed_to_its_anchor():
    """`CAISGEFGSGA` ends on TRBJ2-6's first three templated residues; five were cut.

    `_MAX_FIX = 2` refused it at 5 invented residues. A deletion that restores the conserved
    anchor is exempt, because the alternative is not an ambiguous boundary -- it is a junction
    with no Phe118 at all, which `_canonicalise` refuses anyway. The release ships this as FixAdd.
    """
    m = markup_cdr3("CAISGEFGSGA", "TRBV10-1", "TRBJ2-6", HS)
    assert m.cdr3_repaired == "CAISGEFGSGANVLTF"
    assert m.j_flags == ("add",) and m.j_fix == "FixAdd"


def test_the_anchor_exemption_needs_the_anchor_to_be_genuinely_absent():
    """A junction that ALREADY ends in Phe118 must not be rewritten from germline.

    `_align` gives leading germline gaps away free, so `CAAGGGAWLARISF` on TRAJ33 matches two
    residues and reports `NYQLIW` deleted at dist 0. An exemption keyed on dist alone dropped the
    junction's own Phe and appended six invented residues -- `CAAGGGAWLARISNYQLIW`, exactly the
    failure `_MAX_FIX` exists to stop.
    """
    m = markup_cdr3("CAAGGGAWLARISF", "TRAV23/DV6*01", "TRAJ33*01", HS)
    assert m.cdr3_repaired == "CAAGGGAWLARISF", "an anchored junction is left alone"
    assert "impossible" in m.j_flags


def test_the_anchor_exemption_needs_min_hit_residues_of_evidence():
    """Restoring a germline run needs more than one coincidental residue to stand on.

    `CALRPA` on TRAJ17 has the anchor genuinely missing, so the exemption is eligible -- but the
    germline agreement behind it is below `_MIN_HIT`, and six residues invented off one residue of
    evidence is a guess. The junction is returned untouched and the side says so.
    """
    m = markup_cdr3("CALRPA", "TRAV9-2", "TRAJ17", HS)
    assert m.cdr3_repaired == "CALRPA"
    assert "impossible" in m.j_flags


# --------------------------------------------------------------------------- defect 3
def test_one_residue_of_agreement_is_shallow_not_clean():
    """`CASSQQQQQQQQQF` shares only its terminal Phe with TRBJ1-1's `NTEAFF`.

    That residue is the conserved anchor -- it is there by definition -- so agreeing on it alone
    is not germline corroboration of the call, and `ok` would claim support that does not exist.
    1,986 keys of VDJdb's corpus sit here (1.047 %), 1,425 V and 579 J.

    ⚠ The call also has to SURVIVE. `TRBJ1-5` (`NQPQHF`) collects a match-minus-mismatch score of 3
    here by agreeing on every other residue, which is this junction's Q-richness and not germline
    evidence -- which is why a re-call is decided on contiguous agreement from the anchor
    (`anchor_depth`), where TRBJ1-5 scores 1, exactly what the called TRBJ1-1 scores.
    """
    m = markup_cdr3("CASSQQQQQQQQQF", "TRBV14", "TRBJ1-1", HS)
    assert m.j_call == "TRBJ1-1*01", "an alternating match is not a better call"
    assert "shallow" in m.j_flags
    assert m.v_flags == ("ok",), "the V side genuinely matches CASSQ and must stay clean"


def test_shallow_is_a_caveat_and_does_not_cost_good():
    """`shallow` reports thin evidence; it is not a failure, and nothing disagrees.

    `CGGSARSGELFF` against TRBV9 (`CASSV`) agrees on the Cys and nothing else -- structurally the
    SAME record as the one above -- and there it is the right answer: the V is exonuclease-trimmed
    back to Cys104 and `GGS` is N region, which `test_ambiguous_boundary_is_not_an_error` has
    pinned since 2.16.0. No property of a junction separates the two, so both are `shallow`,
    both keep their boundary, and both stay `good`. Refusing either would delete a standing result.
    """
    m = markup_cdr3("CGGSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert m.v_end == 1 and [e for e in m.errors if e.side == "V"] == []
    assert m.v_flags == ("shallow",) and m.cdr3_repaired == m.cdr3
    assert m.good, "shallow is filterable, not a failure"


# --------------------------------------------------------------------------- defect 4
def test_a_contradicted_call_is_re_called_and_the_sequence_left_alone():
    """Defect 4's own record, answered by step one instead of by a substitution.

    `CAISGEFGSGA` was submitted as TRBV10-1 (`CASSE`), which its `CAIS` contradicts two residues in
    from Cys104. 2.16.0-2.31.0 reported that as `V sub@2 I>S` and 2.32.0 rewrote it. Neither is
    right: TRBV10-3 templates `CAISE`, the junction matches it exactly, and the defect was in the
    CALL. The allele moves, every residue stays, and `allele` is the flag a curator reads.
    """
    m = markup_cdr3("CAISGEFGSGA", "TRBV10-1", "TRBJ2-6", HS)
    assert m.v_call == "TRBV10-3*01" and m.v_flags == ("allele",)
    assert not [e for e in m.errors if e.side == "V"], "nothing to report once the call is right"
    assert m.cdr3_repaired == "CAISGEFGSGANVLTF" and m.j_fix == "FixAdd"


def test_a_side_that_repairs_part_and_declines_part_says_both():
    """A single worst-wins label cannot express this; the flag set can.

    `CAMAGHSGSSPLTL` on TRAJ11 (`NSGYSTLTF`) ends in L where the germline says F, with `LT`
    agreeing behind it -- so the anchor is substituted -- while two residues deeper in the run
    disagree and are only reported. One side, both facts: `sub` and `mismatch`.

    And the name must stay the APPLIED fix: `markup_cdr3` gates the repair on `j_fix in _GOOD`, so
    naming such a side `FailedReplace` does not merely mislabel it, it throws the repair away.
    """
    m = markup_cdr3("CAMAGHSGSSPLTL", "TRAV17*01", "TRAJ11*01", HS)
    assert m.cdr3_repaired == "CAMAGHSGSSPLTF", "the anchor substitution still applies"
    assert m.j_fix == "FixReplace"
    assert set(m.j_flags) == {"sub", "mismatch"}
    applied = [e for e in m.errors if e.side == "J" and e.applied]
    assert [(e.kind, e.dist) for e in applied] == [("sub", 0)], "only the anchor is ever rewritten"


# --------------------------------------------------------------------------- the vocabulary rungs
@pytest.mark.parametrize("call, segment, expect", [
    ("TRAV14", "V", "TRAV14/DV4*01"),        # IMGT files the dual-use genes under a slash
    ("TRAV15D-1-DV6D-1", "V", None),         # mouse; the slash written as a dash
    ("TRBV19-1*01", "V", "TRBV19*01"),       # a `-1` suffix IMGT does not use for this gene
    ("TRBV9*99", "V", "TRBV9*01"),           # an allele IMGT does not mint
    ("TRBJ2.1", "J", "TRBJ2-1*01"),          # VDJdb-era dot nomenclature
])
def test_a_call_is_checked_against_every_gene_name(call, segment, expect):
    from arda.cdr3fix import load_anchors, resolve_allele
    got = resolve_allele(call, segment, load_anchors("human"), "human")
    if expect is not None:
        assert got == expect


@pytest.mark.parametrize("call, segment", [("TRBV6", "V"), ("TRBV12", "V"), ("TRBJ2", "J")])
def test_an_ambiguous_family_is_still_refused(call, segment):
    """The uniqueness rule is what separates resolving from guessing.

    `TRBV6` names five functional genes and `TRBJ2` thirteen. Breaking that tie by taking the
    lowest-numbered candidate is how `TRAV6-7-DV9` became `TRAV6-1*01`.
    """
    from arda.cdr3fix import load_anchors, resolve_allele
    assert resolve_allele(call, segment, load_anchors("human"), "human") == ""


@pytest.mark.parametrize("submitted, expect_flag", [
    ("TRBV9*01", False),        # exact
    ("TRBV9", False),           # a bare gene taking its *01 is not a finding
    ("TRBV9*99", True),         # an allele the submission named and arda did not use
    ("TRBV19-1", True),         # a different gene
])
def test_the_allele_flag_is_a_proofreading_signal(submitted, expect_flag):
    """It has to fire on exactly what a curator must look at, and stay quiet otherwise."""
    m = markup_cdr3("CASSARSGELFF", submitted, "TRBJ2-2*01", HS)
    assert ("allele" in m.v_flags) is expect_flag, (submitted, m.v_flags, m.v_call)


def test_a_multi_candidate_call_takes_the_leading_one():
    """`A+B` behaves like `A,B`, which is what this module has always done with a comma."""
    comma = markup_cdr3("CASSARSGELFF", "TRBV9*01,TRBV10-1*01", "TRBJ2-2*01", HS)
    plus = markup_cdr3("CASSARSGELFF", "TRBV9*01+TRBV10-1*01", "TRBJ2-2*01", HS)
    assert comma.v_call == plus.v_call == "TRBV9*01"


def test_a_boundary_is_never_reported_past_the_end_of_the_junction():
    """`_supported` charges a leading flank the full `_GAP`, so a side carrying framework can
    score its whole germline run away and return 0 -- which as `len - 0` is an index past the
    end. 58 keys of VDJdb's corpus carried one. A side that cannot place a boundary says -1.
    """
    for cdr3, v, j in (("YFCASSQSPGGVAFFGQG", "TRBV14", "TRBJ1-1"),
                       ("CASSQQQQQQQQQF", "TRBV14", "TRBJ1-1"),
                       ("CALRPA", "TRAV9-2", "TRAJ17"),
                       ("CAAGGGAWLARISF", "TRAV23/DV6*01", "TRAJ33*01")):
        m = markup_cdr3(cdr3, v, j, HS)
        n = len(m.cdr3_repaired)
        assert m.v_end == -1 or 0 <= m.v_end <= n, (cdr3, m.v_end, n)
        assert m.j_start == -1 or 0 <= m.j_start < n, (cdr3, m.j_start, n)


def test_flags_reach_the_frame():
    import polars as pl
    from arda.cdr3fix import MARKUP_COLUMNS, markup_batch
    df = pl.DataFrame({"cdr3": ["CASSARSGELFF", "CASSLGGNNKLFF"],
                       "v": ["TRBV9*01", "TRBV9*01"], "j": ["TRBJ2-2*01", "TRBJ1-4*01"],
                       "species": [HS, HS]})
    out = markup_batch(df)
    assert "v_flags" in MARKUP_COLUMNS and "j_flags" in MARKUP_COLUMNS
    assert out["v_flags"].to_list() == ["ok", "ok"]
    assert out["j_flags"].to_list() == ["ok", "mismatch"]


# --------------------------------------------------------------------------- step one: the re-call
def test_a_re_call_needs_two_extra_residues_and_takes_the_allele_not_the_sequence():
    """`_RECALL_GAIN`, on the record the margin was measured against.

    `CGVGGQANTEAFF` was submitted as TRBV29-1 (`CSVGA`), whose `CS` explains two residues; TRBV20-1
    (`CSAR`) explains no more, so nothing clears the margin and 2.32.0's substitution of the G to
    the germline's S is not made. What the record gets is a report, not a rewrite.
    """
    m = markup_cdr3("CGVGGQANTEAFF", "TRBV29-1*01", "TRBJ2-3*01", HS)
    assert m.cdr3_repaired == "CGVGGQANTEAFF", "one residue of gain is not evidence"
    assert "mismatch" in m.v_flags


def test_a_re_call_never_breaks_an_anchor_the_submission_already_holds():
    """`CAARLGNNYKLIW` on TRAJ33, and the reason this needs a rule of its own.

    TRAJ33 (`DSNYQLIW`) explains `WIL` -- three residues CONTIGUOUS FROM ITS OWN ANCHOR, the W
    included. TRAJ12 (`MDSSYKLIF`) explains five, but only by skipping the anchor: its Phe is not
    the junction's Trp. Ranking on depth alone re-called the record and step three then turned its
    W into an F.

    ⚠ `_canonicalise` cannot catch that, and this is the important part: **both anchors are
    canonical**. W and F are each a legitimate [FW]118, so the rewritten junction passes every
    canonicality test there is. Only "do not break an anchor the submission already explains" sees
    it. The converse case is deliberately still allowed -- `CASSKRGGYEQYV` on TRBJ2-7*01 (`SYEQYF`)
    does NOT hold its anchor, so *02 (`SYEQYV`) is free to take the call.
    """
    from arda.cdr3fix import load_anchors
    anchors = load_anchors("human")
    assert anchors[("J", "TRAJ33*01")].templated_aa[-1] == "W"
    assert anchors[("J", "TRAJ12*01")].templated_aa[-1] == "F"
    m = markup_cdr3("CAARLGNNYKLIW", "TRAV13-1*01", "TRAJ33*01", HS)
    assert m.j_call == "TRAJ33*01" and m.cdr3_repaired == "CAARLGNNYKLIW"
    assert "mismatch" in m.j_flags, "the ambiguity is reported, and nothing is changed"

    keep = markup_cdr3("CASSKRGGYEQYV", "TRBV9*01", "TRBJ2-7*01", HS)
    assert keep.j_call == "TRBJ2-7*02" and keep.cdr3_repaired == "CASSKRGGYEQYV"


def test_a_partial_reference_record_neither_wins_nor_loses_a_re_call():
    """IMGT ships `TRBV5-1*02` as a partial sequence: `CAS` against `*01`'s `CASSL`.

    That is a fact about the reference, not about the rearrangement, so comparing their lengths
    would re-call every `*02` record to `*01`. Excluded from the candidate set, and skipped as a
    base, so the submitted allele stands and its boundary is still placed (`TruncatedGermline`).
    """
    m = markup_cdr3("CASSLASGQETQYF", "TRBV5-1*02", "TRBJ2-5*01", HS)
    assert m.v_call == "TRBV5-1*02" and m.v_fix == "TruncatedGermline"


def test_a_pseudogene_is_never_a_re_call_candidate():
    """A junction in hand came from a rearrangement that produced a chain.

    `TRBJ2-2P` (`GRLGG`) won `CASRPGAAGGRPELYF` on depth and trimmed it to `CASRPGAAG` -- canonical
    for that pseudogene, a junction for nothing. ORFs stay candidates: `TRBJ2-7*02` is one.
    """
    m = markup_cdr3("CASRPGAAGGRPELYF", "TRBV6-5*01", "TRBJ2-7*01", HS)
    assert not m.j_call.endswith("P*01") and m.cdr3_repaired != "CASRPGAAG"


def test_a_junction_that_lost_only_its_anchor_gets_it_back():
    """Neither legacy nor 2.32.0 repaired this, and the 2026-06-03 release ships it broken.

    `AQGLLTGGGNKLTF` on TRAV29/DV5 (`CAAS`) agrees with the germline from position 1 on, so the best
    placement is a one-residue coincidence below `_MIN_HIT` and both engines returned
    `FailedNoAlignment`. Prepending the anchor makes `CA` agree, which clears the floor, and the C
    is the one residue whose absence makes the string not a junction at all.
    """
    m = markup_cdr3("AQGLLTGGGNKLTF", "TRAV29/DV5*01", "TRAJ10*01", HS)
    assert m.cdr3_repaired == "CAQGLLTGGGNKLTF" and m.v_fix == "FixAdd"
    assert m.v_canonical and m.good
    (e,) = [e for e in m.errors if e.side == "V" and e.applied]
    assert (e.kind, e.to, e.dist) == ("del", "C", 0)


def test_the_anchor_is_not_restored_off_a_coincidence():
    """The same rule, where it must refuse: `CALRPA` on TRAJ17 stays refused.

    Prepending its Phe still leaves a one-residue hit, so six germline residues would be invented
    off a coincidence. The release refuses it too.
    """
    m = markup_cdr3("CALRPA", "TRAV9-2", "TRAJ17", HS)
    assert m.cdr3_repaired == "CALRPA" and m.j_fix == "FailedNoAlignment"


def test_max_replace_bounds_the_junction_side_and_nothing_inside_the_run():
    """Legacy's `max_replace_size`: junction residues standing where the anchor should be.

    `CAAPGAGSYQTF` on TRAJ28 needs two of them dropped before the germline run picks up, so it is
    refused at the default 1 and repaired at 2. Residues inside the run are unaffected at every
    setting -- that is what `test_a_disagreement_inside_the_run_is_flagged_and_never_changed` pins.
    """
    one = markup_cdr3("CAAPGAGSYQTF", "TRAV8-1*01", "TRAJ28*01", HS)
    assert one.cdr3_repaired == "CAAPGAGSYQTF" and one.j_fix == "FailedReplace"
    two = markup_cdr3("CAAPGAGSYQTF", "TRAV8-1*01", "TRAJ28*01", HS, max_replace=2)
    assert two.j_fix == "FixReplace" and two.cdr3_repaired != "CAAPGAGSYQTF"


def test_the_alleles_a_junction_cannot_separate_travel_with_the_answer():
    """`v_alts` / `j_alts`, and why one allele is not the honest answer.

    `CAISE` is the templated run of TRBV10-3*01, *02 AND *03 alike, so an amino-acid junction cannot
    tell them apart -- and resolving that by functionality and then by name binds a choice with no
    evidence behind it. The chosen call leads, every allele indistinguishable from it follows, and
    the consumer that CAN separate them is the nucleotide stage: `vdjtools.model.infer_nt_batch`
    scores a LIST of alleles per row, so the tie is settled by codon plausibility rather than by
    sort order.
    """
    m = markup_cdr3("CAISGEFGSGA", "TRBV10-1", "TRBJ2-6", HS)
    assert m.v_alts[0] == m.v_call, "the chosen allele leads, so alts[0] and v_call agree"
    assert set(m.v_alts) == {"TRBV10-3*01", "TRBV10-3*02", "TRBV10-3*03"}
    assert m.j_alts and m.j_alts[0] == m.j_call

    # A call the junction confirms still reports itself, so a consumer never has to special-case an
    # empty list.
    clean = markup_cdr3("CASSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert clean.v_alts[0] == "TRBV9*01" and clean.j_alts[0] == "TRBJ2-2*01"


def test_the_alts_reach_the_frame():
    import polars as pl
    from arda.cdr3fix import MARKUP_COLUMNS, markup_batch

    assert "v_alts" in MARKUP_COLUMNS and "j_alts" in MARKUP_COLUMNS
    out = markup_batch(pl.DataFrame({"cdr3": ["CAISGEFGSGA"], "v": ["TRBV10-1"],
                                     "j": ["TRBJ2-6"], "species": [HS]}))
    assert out["v_alts"][0].split(",") == ["TRBV10-3*01", "TRBV10-3*02", "TRBV10-3*03"]
