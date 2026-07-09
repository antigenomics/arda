"""The refbuild locus-coverage manifest — makes a silently-absent locus loud.

``_locus_manifest`` is pure (no IMGT, no mmseqs): it turns the scaffold rows + D germlines a build
produced into one row per DEFINED locus, so a locus that produced nothing is an explicit ``EMPTY``
row rather than an absence. The scenario that motivated it: rat/rabbit/rhesus have no TCR V/J in
IMGT, so those loci were ``continue``-d past with no scaffolds -- yet rabbit/rhesus still ship TRB/TRD
D germlines, which are then unreachable (D lookup is keyed by a hit scaffold's locus).
"""

from __future__ import annotations

from arda.refbuild.build import _locus_manifest
from arda.refbuild.loci import LOCI


def _row(locus, c_call=""):
    return {"locus": locus, "c_call": c_call}


def test_manifest_has_one_row_per_defined_locus():
    m = _locus_manifest([], [])
    assert {r["locus"] for r in m} == {l.name for l in LOCI}
    assert all(r["status"] == "EMPTY" for r in m)      # nothing built -> every locus empty


def test_counts_split_vj_and_jc_and_flag_empty_and_unreachable_d():
    nt_all = [
        _row("IGH"), _row("IGH"), _row("IGH"),          # 3 V-J scaffolds
        _row("IGH", "IGHM"), _row("IGH", "IGHG"),        # 2 J+C scaffolds
        _row("TRB"),                                     # 1 TRB V-J scaffold
    ]
    # TRD has D germlines but NO scaffolds -> unreachable (the rabbit/rhesus failure mode).
    d_germ = [("IGH", "IGHD1-1*01", "GGGACAGGGGGC"),
              ("TRD", "TRDD1*01", "GAAATAGT"), ("TRD", "TRDD2*01", "CCTTCCTAC")]
    m = {r["locus"]: r for r in _locus_manifest(nt_all, d_germ)}

    assert m["IGH"]["n_vj_scaffolds"] == 3 and m["IGH"]["n_jc_scaffolds"] == 2
    assert m["IGH"]["n_d_germlines"] == 1 and m["IGH"]["status"] == "ok"
    assert m["IGH"]["unreachable_d_germlines"] == 0     # IGH has scaffolds, so its D is reachable

    assert m["TRB"]["status"] == "ok" and m["TRB"]["n_vj_scaffolds"] == 1

    assert m["TRD"]["status"] == "EMPTY"
    assert m["TRD"]["n_vj_scaffolds"] == 0 and m["TRD"]["n_jc_scaffolds"] == 0
    assert m["TRD"]["unreachable_d_germlines"] == 2     # 2 D germlines with nowhere to be hit

    # a locus that produced nothing at all and has no D germlines is just EMPTY, not "unreachable_d"
    assert m["IGK"]["status"] == "EMPTY" and m["IGK"]["unreachable_d_germlines"] == 0
