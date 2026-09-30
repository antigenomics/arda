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
    assert set(FLAGS) == {"ok", "allele", "sub", "add", "trim", "shallow", "impossible"}
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
    bad = markup_cdr3("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01", HS)
    assert "impossible" in bad.j_flags and not bad.good


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
    """
    m = markup_cdr3("CASSQQQQQQQQQF", "TRBV14", "TRBJ1-1", HS)
    assert m.j_flags == ("shallow",)
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
def test_a_declined_repair_never_reads_as_a_clean_record():
    """The whole of defect 4, on the record that exposed it.

    `CAISGEFGSGA`'s V side reports `sub@2 I>S` at dist 2, beyond `_MAX_REPLACE`, so it is not
    repaired. `_verdict` reading only the applied edits returned `NoFixNeeded` -- a side claiming
    its junction agrees with germline while carrying the disagreement in its own error list.
    """
    m = markup_cdr3("CAISGEFGSGA", "TRBV10-1", "TRBJ2-6", HS)
    (e,) = [e for e in m.errors if e.side == "V"]
    assert (e.kind, e.dist, e.applied) == ("sub", 2, False)
    assert m.v_flags == ("impossible",) and not m.good


def test_a_side_that_repairs_part_and_declines_part_says_both():
    """A single worst-wins label cannot express this; the flag set can.

    And the name must stay the APPLIED fix: `markup_cdr3` gates the repair on `v_fix in _GOOD`, so
    naming such a side `FailedReplace` does not merely mislabel it, it throws the repair away.
    733 junctions moved that way, `CAAFAGNLLAF` on TRAJ39 among them.
    """
    m = markup_cdr3("CAAFAGNLLAF", "TRAV29/DV5*01", "TRAJ39*01", HS)
    assert m.cdr3_repaired == "CAAFAGNLLTF", "the anchor-adjacent repair still applies"
    assert m.j_fix == "FixReplace"


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
    assert out["j_flags"].to_list() == ["ok", "impossible"]
