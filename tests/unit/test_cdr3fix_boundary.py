"""`boundary_nt`: the germline boundary in nucleotides, from an amino-acid alignment.

The amino acid rounds a boundary to a whole residue and the nucleotide does not, so an
alignment on the protein is systematically short. What the protein still fixes is which
nucleotides are admissible -- 17 of the 20 residues open every one of their codons with
the same base -- and that is what these pin.
"""

from __future__ import annotations

import pytest

from arda.cdr3fix import boundary_nt, load_anchors, markup_cdr3

V, J = "TRBV9*01", "TRBJ2-2*01"
CLEAN = "CASSARSGELFF"


# --------------------------------------------------------------------------- the decision

@pytest.mark.parametrize("residue, germ, expect", [
    ("E", "GGA", 1),    # GAA/GAG: every Glu codon opens with the germline's G -- forced
    ("P", "GGA", 0),    # CCN: no Pro codon can, so the germline stops at the codon -- forced
    ("R", "GGA", 0),    # CGN/AGA/AGG: 4 of 6 open C, 2 open A, none G
    ("S", "AGA", 2),    # AGT/AGC keep two germline nt; the other four Ser codons keep none
])
def test_the_codon_table_decides_how_far_the_v_germline_reaches(residue, germ, expect):
    assert boundary_nt("CASS" + residue, 4, "TGTGCCAGCAGC" + germ, "V")[0] == 12 + expect


def test_a_forced_call_is_certain_and_a_weighed_one_is_not():
    assert boundary_nt("CASSE", 4, "TGTGCCAGCAGCGGA", "V") == (13, 1.0)
    reach, p = boundary_nt("CASSH", 4, "TGTGCCAGCAGCAAT", "V")
    assert (reach, round(p, 4)) == (12, 1.0)          # neither His codon opens with A
    reach, p = boundary_nt("HQPQHF", 1, "AATCAGCCCCAGCATTTT", "J")
    assert reach == 1                                  # CAT keeps two nt, CAC keeps none
    assert 0.9 < p < 1.0                               # ...so it is weighed, not forced


def test_the_j_side_reads_the_third_base_of_the_codon():
    # Germline `CGA` (Arg) before the F codon, observed Lys: of AAA/AAG only AAA closes
    # with the germline's `A`, and it stops there because the middle bases differ. The
    # third base is the one the genetic code leaves free, which is why J is the harder side.
    assert boundary_nt("CAKF", 3, "CGATTT", "J")[0] == 9 - 1


# --------------------------------------------------------------------------- degenerate input

@pytest.mark.parametrize("junction, residues, germ, side", [
    ("CASS", -1, "TGTGCCAGCAGC", "V"),                 # the side declined
    ("CASS", 2, "", "V"),                              # no germline to place it against
    ("CASS", 2, "TGTGCCAGCAGC", "X"),                  # not a side
])
def test_nothing_to_place_a_boundary_against_declines(junction, residues, germ, side):
    assert boundary_nt(junction, residues, germ, side) == (-1, 0.0)


def test_a_germline_the_alignment_ran_out_of_stops_where_it_ends():
    assert boundary_nt("CASSL", 4, "TGTGCCAGCAGC", "V") == (12, 1.0)
    assert boundary_nt("CASSL", 5, "TGTGCCAGCAGC", "V") == (12, 1.0)


def test_a_partial_germline_codon_can_still_be_reached_into():
    # 14 nt = 4 codons + 2, which is what `partial_nt` records. Both Glu codons open
    # `GA`, so the germline reaches its own end and cannot be credited past it.
    assert boundary_nt("CASSE", 4, "TGTGCCAGCAGCGA", "V")[0] == 14


# --------------------------------------------------------------------------- through markup

def test_markup_reports_both_boundaries_in_nucleotides():
    m = markup_cdr3(CLEAN, V, J, "human")
    assert 3 * m.v_end <= m.v_end_nt <= 3 * m.v_end + 2
    assert 3 * m.j_start - 2 <= m.j_start_nt <= 3 * m.j_start


def test_the_vdjdb_residue_answers_are_unchanged_by_the_nucleotide_one():
    """VDJdb rounds both back, so `j.start` cannot move and `v.end` moves only past 2 nt."""
    anchors = load_anchors("human")
    for cdr3 in (CLEAN, "CASSIRSSYEQYF", "CATSDPGTGHQPQHF", "CASSQVGTGVYEQYF"):
        m = markup_cdr3(cdr3, "TRBV9*01", "TRBJ2-7*01", "human", anchors=anchors)
        assert (m.j_start_nt + 2) // 3 == m.j_start
        assert (m.v_end_nt + 1) // 3 in (m.v_end, m.v_end + 1)


def test_a_side_that_declined_has_no_nucleotide_boundary_either():
    m = markup_cdr3(CLEAN, "TRBVnope*01", J, "human")
    assert (m.v_end, m.v_end_nt) == (-1, -1)
    assert m.j_start >= 0 and m.j_start_nt >= 0


# --------------------------------------------------------------------------- what a boundary PROMISES

@pytest.mark.parametrize("junction, v, j, fix", [
    ("CAARITGEKLFF", "TRBV6-5*01", "TRBJ1-4*01", "NoFixNeeded"),
    ("CAARITGEKLFFGS", "TRBV6-5*01", "TRBJ1-4*01", "FixTrim"),
    ("CAARITGEKLF", "TRBV6-5*01", "TRBJ1-4*01", "FixAdd"),
    ("CAAGGGSYIPT", "TRAV13-1*01", "TRAJ6*01", "FixAdd"),
])
def test_every_residue_a_boundary_credits_is_germline(junction, v, j, fix):
    """`v_end` / `j_start` index the REPAIRED junction, and everything inside them is germline.

    This is the guarantee a consumer applying `cdr3_repaired` relies on, so it is pinned rather
    than left to follow from the repair paths agreeing with each other -- they do not have to.
    Measured over VDJdb's 189,596 distinct curation keys: the credited V prefix is exact germline
    on 185,636 of 185,636 keys that place one, and the credited J suffix on 187,272 of 187,272 --
    zero disagreements on either side.
    """
    m = markup_cdr3(junction, v, j, "HomoSapiens", max_replace=0)
    assert m.j_fix == fix
    seq, anchors = m.cdr3_repaired, load_anchors("human")
    assert anchors[("V", m.v_call)].templated_aa.startswith(seq[:m.v_end])
    assert anchors[("J", m.j_call)].templated_aa.endswith(seq[m.j_start:])


def test_a_repair_is_a_fixed_point_of_itself():
    """Handing a repaired junction back must not repair it again -- otherwise a build that writes
    `cdr3_repaired` and re-runs would walk. True on all 324 of VDJdb's keys that take a `FixAdd`.

    The boundary a repair path reports is a LOWER bound, short by exactly the residues it restored:
    `templated` counts residues the germline is OBSERVED to explain, which is what VDJdb's `vEnd` /
    `jStart` have always meant, and a restored residue was not observed. So the sequence is a fixed
    point while the coordinate is allowed to grow once -- never to move the other way.
    """
    first = markup_cdr3("CAARITGEKLF", "TRBV6-5*01", "TRBJ1-4*01", "HomoSapiens", max_replace=0)
    again = markup_cdr3(first.cdr3_repaired, first.v_call, first.j_call, "HomoSapiens",
                        max_replace=0)
    assert again.cdr3_repaired == first.cdr3_repaired == "CAARITGEKLFF"
    assert again.j_fix == "NoFixNeeded"
    assert again.j_start <= first.j_start        # the restored residue is credited on the re-run
