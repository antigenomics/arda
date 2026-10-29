"""Raw IMGT partial flags outrank same-gene length agreement."""
import logging
from types import SimpleNamespace

import pytest

from arda.refbuild import build, imgt


@pytest.mark.parametrize("segment,label,excluded", [
    ("V", "partial in 3'", True), ("V", "partial in 5'", False),
    ("J", "partial in 5'", True), ("J", "partial in 3'", False),
    ("D", "partial in 5'", True), ("D", "partial in 3'", True),
    ("V", "partial in 5' and in 3'", True), ("J", "partial", True),
    ("V", "", False),
])
def test_junction_side_completeness(segment, label, excluded):
    assert imgt.junction_incomplete(label, segment) == excluded


def _raw(allele, label):
    return "|".join(["accession", allele, "species", "F", "V-REGION"] + [""] * 8 + [label, "", ""])


def test_annotation_loader_retains_partial_alignment_targets(tmp_path, monkeypatch):
    raw = tmp_path / "raw.fasta"
    raw.write_text("".join(f">{_raw(a, label)}\nTGTAAA\n" for a, label in
                           [("TRAV1*01", "partial in 3'"), ("TRAV2*01", "partial in 5'"), ("TRAV3*01", "")]))
    ungapped = tmp_path / "ungapped.fasta"
    ungapped.write_text("".join(f">TRAV{i}*01\nTGTAAA\n" for i in (1, 2, 3)))
    monkeypatch.setattr(imgt, "gene_fasta_path", lambda *a: raw)
    monkeypatch.setattr(imgt, "ungap_gene", lambda *a: ungapped)
    assert set(imgt.load_functional_alleles("sp", "TR", "TRAV")) == {"TRAV1*01", "TRAV2*01", "TRAV3*01"}


def test_anchor_metadata_preserves_partial_even_when_lengths_equal(tmp_path, monkeypatch):
    v, j = tmp_path / "v.fasta", tmp_path / "j.fasta"
    v.write_text("".join(f">{_raw(a, label)}\nTGTAAA\n" for a, label in
                         [("TRAV1*01", ""), ("TRAV1*02", "partial in 3'"), ("TRAV2*01", "partial in 5'")]))
    j.write_text("".join(f">{_raw(a, label)}\nTTTGGG\n" for a, label in
                         [("TRAJ1*01", "partial in 5'"), ("TRAJ2*01", "partial in 3'")]))
    locus = SimpleNamespace(name="TRA", group="TR", v="TRAV", j="TRAJ")
    monkeypatch.setattr(build, "LOCI", [locus])
    monkeypatch.setattr(imgt, "gene_fasta_path", lambda sp, group, stem: v if stem == "TRAV" else j)
    monkeypatch.setattr(build.combinations, "load_v_fwr3_stops", lambda _: {})
    monkeypatch.setattr(build.combinations, "load_j_fr4_offsets", lambda _: {"TRAJ1*01": (-1, 0), "TRAJ2*01": (-1, 0)})
    monkeypatch.setattr(build, "_v_anchor", lambda *a: (0, "ndm"))
    rows = {r["allele"]: r for r in build._collect_cdr3_anchors("human", "sp", logging.getLogger("test"))}
    assert rows["TRAV1*02"]["status"] == rows["TRAJ1*01"]["status"] == "partial"
    assert rows["TRAV1*02"]["imgt_completeness"] == "partial in 3'"
    assert all(rows[a]["status"] == "ok" for a in ("TRAV1*01", "TRAV2*01", "TRAJ2*01"))


def test_d_collector_preserves_preliminary_alignment_records(tmp_path, monkeypatch):
    raw, ungapped = tmp_path / "raw.fasta", tmp_path / "d.fasta"
    partial_header = _raw("TRBD1*01", "partial in 3'")
    complete_header = _raw("TRBD2*01", "")
    raw.write_text(f">{partial_header}\nGGG\n>{complete_header}\nAAA\n")
    ungapped.write_text(">TRBD1*01\nGGG\n>TRBD2*01\nAAA\n")
    monkeypatch.setattr(build, "VDJ_LOCI", [SimpleNamespace(name="TRB", group="TR", d="TRBD")])
    monkeypatch.setattr(imgt, "gene_fasta_path", lambda *a: raw)
    monkeypatch.setattr(imgt, "ungap_gene", lambda *a: ungapped)
    assert build._collect_d_germlines("sp", logging.getLogger("test")) == [("TRB", "TRBD1*01", "GGG"), ("TRB", "TRBD2*01", "AAA")]
    rows = build._d_germline_metadata("sp", [("TRB", "TRBD1*01", "GGG"), ("TRB", "TRBD2*01", "AAA")])
    assert [r["status"] for r in rows] == ["partial", "ok"]
    assert rows[0]["imgt_completeness"] == "partial in 3'"


def test_complete_duplicate_is_not_discarded(monkeypatch):
    monkeypatch.setattr(build, "VDJ_LOCI", [SimpleNamespace(name="TRB", group="TR", d="TRBD")])
    monkeypatch.setattr(build, "_gapped_alleles", lambda *a: [
        ("D1", "F", "AAA", "partial in 3'"), ("D1", "F", "AAA", "")])
    rows = build._d_germline_metadata("sp", [("TRB", "D1", "AAA")])
    assert rows[0]["status"] == "ok"


def test_raw_partial_remains_preliminary_placeable_evidence():
    from arda.cdr3fix import Anchor, _placeable

    anchor = Anchor(locus="TRB", segment="V", templated_aa="CASS", functionality="F",
                    status="partial", source="ndm")
    assert anchor.status != "ok"
    assert _placeable(anchor)
