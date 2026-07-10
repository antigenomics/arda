"""Junction markup + repair of bare (junction_aa, V, J) records.

Everything here is JUNCTION space: Cys104 .. [FW]118, both included -- the
convention VDJdb's ``cdr3`` column uses. These tests pin that convention, the
per-allele anchors, and the deliberate split between *detecting* a germline
mismatch and *repairing* it.
"""

from __future__ import annotations

import pytest

from arda.cdr3fix import (
    load_anchors, markup_cdr3, resolve_allele, resolve_locus, resolve_species,
)
from tests.conftest import requires_human_db

pytestmark = [requires_human_db]

HS = "HomoSapiens"


@pytest.fixture(scope="module")
def anchors():
    a = load_anchors("human")
    if not a:
        pytest.skip("cdr3_anchors.tsv not built (run `arda build-db`)")
    return a


# --------------------------------------------------------------------------
# Anchors
# --------------------------------------------------------------------------

@pytest.mark.parametrize("segment,allele,anchor_nt,templated", [
    # V: Cys104. TRBV1*01/IGHV1-18*01 also match OLGA's V_gene_CDR3_anchors.csv.
    ("V", "TRBV1*01", 267, "CTSSQ"),
    ("V", "IGHV1-18*01", 285, "CAR"),
    # OLGA says 267 here, but its germline is not IMGT's: 270 is the Cys codon.
    ("V", "TRBV3-1*01", 270, "CASSQ"),
    # Cross-check against arda's own scaffold TRB_0 (junction_aa = CASSXXXSTDTQYF).
    ("V", "TRBV4-3*04", 219, "CASS"),
    # J: [FW]118. TRBJ1-4 is a double-terminal-F gene.
    ("J", "TRBJ1-4*01", 20, "TNEKLFF"),
    ("J", "TRBJ2-3*01", 18, "STDTQYF"),
    # TRAJ31*01's aux frame column disagrees with its own cdr3_stop; frame = anchor % 3.
    ("J", "TRAJ31*01", 23, "NNNARLMF"),
])
def test_anchor_golden_values(anchors, segment, allele, anchor_nt, templated):
    a = anchors[(segment, allele)]
    assert a.status == "ok"
    assert a.templated_aa == templated
    from arda.paths import vdj_dir
    import polars as pl
    df = pl.read_csv(vdj_dir("human") / "cdr3_anchors.tsv", separator="\t", infer_schema_length=0)
    row = df.filter((pl.col("segment") == segment) & (pl.col("allele") == allele)).row(0, named=True)
    assert int(row["anchor_nt"]) == anchor_nt


def test_v_anchor_is_always_a_cysteine(anchors):
    """A V anchor off by one codon corrupts every coordinate. Refuse, never guess."""
    bad = [a for (seg, a), v in anchors.items()
           if seg == "V" and v.status == "ok" and not v.templated_aa.startswith("C")]
    assert bad == []


def test_j_anchor_is_conserved_for_functional_alleles(anchors):
    """[FW]118 holds for functional J alleles; TRAJ35*01 genuinely has Cys118 in IMGT."""
    bad = sorted(a for (seg, a), v in anchors.items()
                 if seg == "J" and v.status == "ok" and v.functionality == "F"
                 and not v.templated_aa.endswith(("F", "W")))
    assert bad == ["TRAJ35*01"]


def test_truncated_allele_is_flagged_not_guessed(anchors):
    """IGHV1-18*02 never reaches Cys104 -> no anchor, and markup fails cleanly."""
    assert anchors[("V", "IGHV1-18*02")].status == "no_anchor"
    m = markup_cdr3("CASSARSGELFF", "IGHV1-18*02", "TRBJ2-2*01", HS)
    assert m.v_fix == "FailedBadSegment"
    assert m.v_end == -1
    assert not m.good


def test_pseudogene_alleles_are_kept(anchors):
    """VDJdb cites ORF/pseudogene alleles, so anchors must not gate on functionality."""
    assert ("J", "TRBJ2-2P*01") in anchors


# --------------------------------------------------------------------------
# Resolution
# --------------------------------------------------------------------------

def test_resolve_species():
    assert resolve_species("HomoSapiens") == "human"
    assert resolve_species("MusMusculus") == "mouse"
    assert resolve_species("MacacaMulatta") == "rhesus_monkey"
    assert resolve_species("human") == "human"


def test_resolve_locus_handles_dual_gene():
    assert resolve_locus("TRBV6-1", "TRBJ2-2") == "TRB"
    assert resolve_locus("TRAV29/DV5", "TRAJ6") == "TRA"   # leading token wins, not TRD
    assert resolve_locus("", "IGHJ4") == "IGH"


def test_resolve_allele_gene_level_fallback(anchors):
    assert resolve_allele("TRBV9", "V", anchors) == "TRBV9*01"
    assert resolve_allele("TRBV9*01", "V", anchors) == "TRBV9*01"
    assert resolve_allele("NOSUCHV", "V", anchors) == ""


# --------------------------------------------------------------------------
# Markup + repair — the three motivating cases
# --------------------------------------------------------------------------

def test_clean_record_matches_vdjdb_boundaries():
    """VDJdb's own cdr3fix for this record: vEnd=4, jStart=7, NoFixNeeded."""
    m = markup_cdr3("CASSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert (m.v_end, m.j_start) == (4, 7)
    assert m.v_fix == m.j_fix == "NoFixNeeded"
    assert m.errors == [] and not m.fix_needed and m.good


def test_substitution_next_to_the_cys_is_repaired():
    """CASS -> CCSS: substitution at index 1, adjacent to Cys104, so it is repaired."""
    m = markup_cdr3("CCSSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert m.cdr3_repaired == "CASSARSGELFF"
    assert m.v_fix == "FixReplace"
    (e,) = [e for e in m.errors if e.side == "V"]
    assert (e.kind, e.pos, e.length, e.frm, e.to, e.dist, e.applied) == \
        ("sub", 1, 1, "C", "A", 1, True)


def test_missing_terminal_f_is_added_never_truncated():
    """NEKLFF -> NEKLF: the double-F J gene lost its [FW]118. Repair adds it back."""
    m = markup_cdr3("CASSLGGNEKLF", "TRBV9*01", "TRBJ1-4*01", HS)
    assert m.cdr3_repaired == "CASSLGGNEKLFF"
    assert m.j_fix == "FixAdd"
    (e,) = [e for e in m.errors if e.side == "J"]
    assert (e.kind, e.length, e.to, e.dist, e.applied) == ("del", 1, "F", 0, True)


def test_deep_substitution_is_reported_but_not_repaired():
    """NEKLFF -> NNKLFF: 4 residues from the anchor, so it is indistinguishable from
    the N region starting early. Report where it is; do not rewrite the record."""
    m = markup_cdr3("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01", HS)
    assert m.cdr3_repaired == m.cdr3          # untouched
    assert m.j_fix == "NoFixNeeded"
    (e,) = [e for e in m.errors if e.side == "J"]
    assert (e.kind, e.pos, e.frm, e.to, e.dist, e.applied) == ("sub", 8, "N", "E", 4, False)
    # ...but the caller can opt in.
    deep = markup_cdr3("CASSLGGNNKLFF", "TRBV9*01", "TRBJ1-4*01", HS, max_replace=4)
    assert deep.cdr3_repaired == "CASSLGGNEKLFF" and deep.j_fix == "FixReplace"


def test_ambiguous_boundary_is_not_an_error():
    """CASS vs CGGS: two mismatches never pay for themselves, so v_end stops at the
    Cys and GGS is correctly read as N region, not as two typos."""
    m = markup_cdr3("CGGSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert m.v_end == 1
    assert [e for e in m.errors if e.side == "V"] == []
    assert m.cdr3_repaired == m.cdr3


def test_repair_is_idempotent():
    m = markup_cdr3("CCSSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    again = markup_cdr3(m.cdr3_repaired, "TRBV9*01", "TRBJ2-2*01", HS)
    assert not again.fix_needed and again.cdr3_repaired == m.cdr3_repaired


def test_cdr3fix_json_matches_vdjdb_key_set():
    m = markup_cdr3("CASSARSGELFF", "TRBV9*01", "TRBJ2-2*01", HS)
    assert set(m.to_cdr3fix()) == {
        "cdr3", "cdr3_old", "fixNeeded", "good", "jCanonical", "jFixType", "jId",
        "jStart", "vCanonical", "vEnd", "vFixType", "vId"}
    assert m.to_cdr3fix()["vEnd"] == 4 and m.to_cdr3fix()["jStart"] == 7


def test_non_canonical_junction_is_reported_not_hidden():
    m = markup_cdr3("XASSARSGELFX", "TRBV9*01", "TRBJ2-2*01", HS)
    assert not m.v_canonical and not m.j_canonical
