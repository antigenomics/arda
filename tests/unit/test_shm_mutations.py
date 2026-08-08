"""``v_mutations`` / ``j_mutations``: the SHM record a lineage-tree tool consumes.

The information was already in the output -- ``sequence_alignment`` and ``germline_alignment`` carry
every column, and on a real bulk IG library the germline they report matches the shipped allele on
28,365 of 28,365 mapped reads. What these tests pin is the part that is easy to get wrong once it is
summarised: the COORDINATE FRAME (each segment's own germline, not the scaffold and not the read),
and the SCOPE.

⛔ The scope is the domain constraint. V(D)J recombination chews the segment ends back and adds
non-templated N/P bases, so a mismatch inside the V..J interior is not attributable to any germline
-- the V-end / NDN / J-start partition of a junction is frequently not identifiable from the
sequence at all. So an NDN position must be excluded BY CONSTRUCTION, not by a downstream filter
someone can forget: the pad is not a segment, so it has no germline coordinate to be recorded under.
``test_a_mutation_in_the_ndn_is_not_recorded`` is the test that says so.

Pure: no reference DB, no mmseqs. The alignments are written out as literal strings, which is also
what makes the germline coordinates checkable by eye.
"""

from __future__ import annotations

import pytest

from arda.annotate.cigar import _segment_cigars_py, segment_cigars

IMPLS = pytest.mark.parametrize("impl", [segment_cigars, _segment_cigars_py],
                                ids=["shipped", "python-reference"])


@IMPLS
def test_known_mutations_at_known_germline_positions(impl):
    """Germline base, 1-based germline position, read base -- in that order, in germline order."""
    germ = "ACGTACGTACGTACGTACGT"
    read = "ACGTACGTTCGTACGAACGT"
    #                ^8 A>T      ^16 T>A
    out = impl(read, germ, 1, 1, len(read), t_vend=20, t_jstart=0, t_vjend=0)
    assert out["v_mutations"] == "A9T,T16A"
    assert out["v_cigar"] == "20M"
    assert "j_mutations" not in out


@IMPLS
def test_positions_are_germline_not_alignment_columns(impl):
    """A 5'-truncated read starts at germline 101, so its first column is position 101 -- not 1.

    Getting this wrong is invisible on a full-length amplicon read and wrong on every bulk read,
    and two reads of the same clone would then carry mutation sets that cannot be compared.
    """
    germ = "ACGTACGTAC"
    read = "ACGTAGGTAC"
    out = impl(read, germ, 1, 101, len(read), t_vend=300, t_jstart=0, t_vjend=0)
    assert out["v_mutations"] == "C106G"


@IMPLS
def test_j_positions_are_in_the_j_allele_frame(impl):
    """The J germline starts at scaffold ``t_jstart``, so J position == t - t_jstart + 1."""
    # V germline 1..4, N-pad 5..8, J germline 9..14 (so t_jstart=9 is J position 1).
    germ = "ACGT" + "NNNN" + "TTCGGA"
    read = "ACGT" + "GACA" + "TTCTGA"
    #                            ^ J position 4, T>G
    out = impl(read, germ, 1, 1, len(read), t_vend=4, t_jstart=9, t_vjend=14)
    assert out["j_mutations"] == "G4T"
    assert "v_mutations" not in out                       # the V part is germline


@IMPLS
def test_a_mutation_in_the_ndn_is_not_recorded(impl):
    """⛔ The domain rule. Every one of the four pad columns differs from the scaffold's N, and the
    read is otherwise pure germline: nothing may be reported, for either segment."""
    germ = "ACGTACGT" + "NNNN" + "TTCGGATT"
    read = "ACGTACGT" + "GACA" + "TTCGGATT"
    out = impl(read, germ, 1, 1, len(read), t_vend=8, t_jstart=13, t_vjend=20)
    assert "v_mutations" not in out
    assert "j_mutations" not in out
    assert out["v_cigar"] == "8M12S"                       # the alignment itself is unaffected


@IMPLS
def test_an_n_in_the_read_is_a_no_call_not_a_mutation(impl):
    """A base the sequencer did not call is not evidence of a somatic mutation."""
    germ = "ACGTACGTAC"
    read = "ACGTANGTAC"
    out = impl(read, germ, 1, 1, len(read), t_vend=10, t_jstart=0, t_vjend=0)
    assert "v_mutations" not in out


@IMPLS
def test_constant_region_mismatches_are_not_shm(impl):
    """The C segment is real germline, but a CH1 mismatch is not V(D)J SHM and has no place in a
    lineage frame; only V and J get mutation lists."""
    germ = "TTCGGA" + "ACGTACGT"
    read = "TTCGGA" + "ACGTAAGT"
    out = impl(read, germ, 1, 1, len(read), t_vend=0, t_jstart=1, t_vjend=6)
    assert out["c_cigar"] == "6S8M"
    assert "c_mutations" not in out
    assert "j_mutations" not in out


@IMPLS
def test_a_deletion_shifts_the_germline_frame_for_everything_after_it(impl):
    """An SHM deletion of 2 germline bases. The substitution after it is at its GERMLINE position,
    which is 2 further along than its position in the read -- the failure mode a naive
    read-offset implementation has, and the reason indels get their own test."""
    germ = "ACGTACGTACGT"
    read = "ACGT--GTACGA"
    #             deletion of germline 5..6      ^ germline 12, T>A
    out = impl(read, germ, 1, 1, len(read.replace("-", "")), t_vend=12, t_jstart=0, t_vjend=0)
    assert out["v_cigar"] == "4M2D6M"
    assert out["v_mutations"] == "T12A"


@IMPLS
def test_an_insertion_consumes_no_germline_position(impl):
    """Inserted read bases have no germline coordinate, so they are in the CIGAR (``I``) and not in
    the mutation list, and they do not advance the germline frame."""
    germ = "ACGT--ACGTAC"
    read = "ACGTTTACGTAG"
    #           inserted                     ^ germline 10, C>G
    out = impl(read, germ, 1, 1, len(read), t_vend=10, t_jstart=0, t_vjend=0)
    assert out["v_cigar"] == "4M2I6M"
    assert out["v_mutations"] == "C10G"


@IMPLS
def test_mutations_replay_onto_the_germline_and_reproduce_the_read(impl):
    """The round-trip a tree builder actually performs: germline + mutation set == observed
    sequence over the covered span. This is the property; the string format is an encoding of it."""
    germ = "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCC"
    read = "CAGGTTCAGCTAGTGCAGTCTGGAGCTGAGGTCAAGAAGCCTGGCGCC"
    tstart = 20                                            # germline 20..67
    out = impl(read, germ, 1, tstart, len(read), t_vend=200, t_jstart=0, t_vjend=0)

    rebuilt = list(germ)
    for m in out["v_mutations"].split(","):
        pos = int(m[1:-1])
        assert germ[pos - tstart] == m[0]                  # the germline base the record claims
        rebuilt[pos - tstart] = m[-1]
    assert "".join(rebuilt) == read


@IMPLS
def test_no_mutations_emits_no_column_at_all(impl):
    """A germline-identical read must not carry an empty-string field that reads as a measurement."""
    germ = read = "ACGTACGTACGT"
    out = impl(read, germ, 1, 1, len(read), t_vend=12, t_jstart=0, t_vjend=0)
    assert "v_mutations" not in out and "j_mutations" not in out
