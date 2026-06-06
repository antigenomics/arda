"""Unit tests for native translation and frame utilities."""

from arda.refbuild.translate import (
    translate,
    detect_coding_frame,
    reverse_complement,
    back_translate,
    aa_coords_from_nt,
    CODON_TABLE,
)


def test_codon_table_complete():
    assert len(CODON_TABLE) == 64
    assert CODON_TABLE["ATG"] == "M"
    assert CODON_TABLE["TGG"] == "W"
    assert CODON_TABLE["TAA"] == CODON_TABLE["TAG"] == CODON_TABLE["TGA"] == "*"


def test_translate_basic():
    assert translate("ATGTGGTAA") == "MW*"
    # frame offset
    assert translate("CATGTGG", 1) == "MW"


def test_translate_n_is_x():
    assert translate("ATGNNNTGG") == "MXW"


def test_detect_coding_frame_stop_free():
    # A clean ORF reads in frame 0.
    assert detect_coding_frame("ATGAAACCCGGGTGG") == 0
    # Frame 0 has a stop (TAA); detection must pick a stop-free frame instead.
    seq = "ATGTAAATGGGG"
    assert translate(seq, 0).count("*") > 0
    f = detect_coding_frame(seq)
    assert translate(seq, f).count("*") == 0


def test_aa_coords_from_nt():
    # coding starts at nt 1; region nt 4..9 -> aa 2..3
    assert aa_coords_from_nt(4, 9, 1) == (2, 3)
    # coding starts at nt 3 (frame origin shifted)
    assert aa_coords_from_nt(3, 8, 3) == (1, 2)


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAN") == "NTTT"
    # double revcomp is identity
    s = "ACGTTGCANNACGT"
    assert reverse_complement(reverse_complement(s)) == s


def test_back_translate_human_codons():
    # M->ATG, W->TGG, C->TGC (most-frequent human codons)
    assert back_translate("MWC") == "ATGTGGTGC"
    # unknown residue -> NNN
    assert back_translate("X") == "NNN"
    # back-translation then translation recovers the protein
    aa = "MKWLV"
    assert translate(back_translate(aa)) == aa
