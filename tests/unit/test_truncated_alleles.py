"""Germline alleles that are truncated INTO the junction must not build scaffolds.

DB-free by design. The pre-existing coverage for this class of bug
(``tests/synthetic/test_locus_disambiguation.py``) is gated on ``requires_human_db`` and is
therefore *skipped in CI* -- which is exactly how the shared-TRAV/DV bug regressed once already.
These tests must never need a built reference.

The bug: a V allele whose germline stops short of the canonical 3' end builds a ``V + pad + J``
scaffold that LACKS nucleotides a real rearrangement has, so every read projected onto it yields a
junction short by exactly that much. ``TRAV20*03`` templates 3 nt into the junction where its
siblings template 13 -- so every clonotype built on it came out 3 aa short, and arda's TRA CDR3
length distribution peaked 3 aa below every other tool's. The junction still starts with Cys and
ends with [FW], so it looks canonical and nothing downstream catches it.
"""

import logging

import pytest

from arda.refbuild import build

log = logging.getLogger("test")


def test_truncated_allele_is_dropped():
    # *02 templates one codon less into the junction than *01 -> it can only ever emit short
    # junctions, so it must not reach the scaffold builder.
    alleles = {"TRAV1*01": "AAA", "TRAV1*02": "AAA"}
    kept = build.drop_truncated(alleles, {"TRAV1*01": 13, "TRAV1*02": 10},
                                logger=log, locus="TRA", segment="V")
    assert set(kept) == {"TRAV1*01"}


def test_truncation_is_relative_to_the_gene_so_no_gene_is_ever_orphaned():
    # The predicate is "short of your OWN gene's longest allele", so that allele always survives.
    # This is what makes the filter orphan-safe by construction -- verified across all 5 shipped
    # organisms: 0 genes lose every allele.
    alleles = {f"TRBV5*0{i}": "AAA" for i in (1, 2, 3)}
    kept = build.drop_truncated(alleles, {"TRBV5*01": 4, "TRBV5*02": 4, "TRBV5*03": 4},
                                logger=log, locus="TRB", segment="V")
    assert set(kept) == set(alleles), "alleles that agree with their gene must all survive"

    kept = build.drop_truncated(alleles, {"TRBV5*01": 2, "TRBV5*02": 2, "TRBV5*03": 2},
                                logger=log, locus="TRB", segment="V")
    assert set(kept) == set(alleles), "a uniformly short gene keeps every allele (nothing to compare to)"


def test_a_gene_is_judged_only_against_itself():
    # Different genes have wildly different V-templated CDR3 lengths (1-17 aa across loci), so a
    # global length cutoff would be nonsense. Only same-gene alleles are comparable.
    alleles = {"TRAV1*01": "AAA", "TRAV2*01": "AAA"}
    kept = build.drop_truncated(alleles, {"TRAV1*01": 13, "TRAV2*01": 3},
                                logger=log, locus="TRA", segment="V")
    assert set(kept) == set(alleles), "a short gene is not a truncated allele"


def test_unanchored_alleles_are_kept_for_mapping():
    # tail = -1 means Cys104 is not findable at all. Dropping those would delete 41 *functional*
    # mouse IGHV/IGLV genes, so they keep their scaffold (the read still maps) and the junction is
    # suppressed downstream instead. See _process_locus.
    alleles = {"IGHV1S21*01": "AAA", "IGHV1S21*02": "AAA"}
    kept = build.drop_truncated(alleles, {"IGHV1S21*01": -1, "IGHV1S21*02": -1},
                                logger=log, locus="IGH", segment="V")
    assert set(kept) == set(alleles)


def test_one_allele_per_gene_prefers_star01():
    alleles = {"IGHV1-2*01": "A", "IGHV1-2*02": "A", "IGHV1-2*03": "A"}
    kept = build.select_one_allele_per_gene(alleles, logger=log, locus="IGH", segment="V")
    assert set(kept) == {"IGHV1-2*01"}


def test_one_allele_per_gene_keeps_genes_that_have_no_star01():
    # 19 human V genes (IGHV2-70D, IGHV3-25, IGHV3-43D, ...) have no *01 record. A literal "*01 only"
    # filter would delete them silently; we fall back to the lowest-numbered allele instead.
    alleles = {"IGHV3-25*02": "A", "IGHV3-25*03": "A", "IGHV1-2*01": "A"}
    kept = build.select_one_allele_per_gene(alleles, logger=log, locus="IGH", segment="V")
    assert set(kept) == {"IGHV3-25*02", "IGHV1-2*01"}


@pytest.mark.parametrize("tail_nt,expect", [(13, True), (3, False)])
def test_v_cys_tail_measures_into_the_junction(tail_nt, expect):
    # v_cys_tail is the nt the V templates AFTER Cys104 -- not the raw record length. Using the raw
    # length conflates 5'-truncation (harmless: the junction is downstream) with 3'-truncation (the
    # bug), and ~50 of a naive 108-allele shortlist were 5'-partial false positives.
    cys = "TGT"                       # Cys codon
    v = "GGG" * 20 + cys + "A" * (tail_nt - 3)
    fwr3_stop = len("GGG" * 20) + 3   # 1-based-ish stop that puts the anchor on the Cys
    assert build.v_cys_tail(v, fwr3_stop) == tail_nt
    assert (tail_nt >= 13) is expect
