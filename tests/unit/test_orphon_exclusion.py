import pytest
from arda.refbuild import imgt


def test_orphons_are_excluded_from_the_reference(tmp_path, monkeypatch):
    """⛔ An orphon sits OUTSIDE its locus and cannot rearrange, so a call naming one on a
    rearranged read is wrong by construction. Measured: `TRBV20-1` -> `TRBV20/OR9-2` was 98.51 %
    of all V disagreements on a MIGEC IGH library, taking v_gene from .9963 down to .7500, and it
    was the ONLY call emitted on those reads -- unfixable downstream."""
    fasta = tmp_path / "g.fasta"
    fasta.write_text(">TRBV20-1*01\nACGTACGTAC\n>TRBV20/OR9-2*01\nACGTACGTAA\n")
    monkeypatch.setattr(imgt, "gene_fasta_path", lambda *a: fasta)
    monkeypatch.setattr(imgt, "ungap_gene", lambda *a: fasta)
    monkeypatch.setattr(imgt, "parse_functionality",
                        lambda p: {"TRBV20-1*01": "F", "TRBV20/OR9-2*01": "F"})
    got = imgt.load_functional_alleles("sp", "TR", "TRBV")
    assert set(got) == {"TRBV20-1*01"}, "the orphon must not reach the reference"
    kept = imgt.load_functional_alleles("sp", "TR", "TRBV", keep_orphons=True)
    assert set(kept) == {"TRBV20-1*01", "TRBV20/OR9-2*01"}
