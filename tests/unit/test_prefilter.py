"""The k-mer prefilter's contract: it may pass junk, it must not drop a read that can align.

The prefilter is the one component allowed to lose reads, so every test here is about bounding
what it loses. A false positive costs a search that would have happened anyway; a false negative
is a read arda silently never reports.
"""
from __future__ import annotations

import random

import pytest

from arda import prefilter

pytest.importorskip("arda._prefilter", reason="native prefilter extension not built")

from arda._prefilter import Prefilter  # noqa: E402


def test_a_read_that_is_a_substring_of_a_target_always_passes():
    """The floor: an exact fragment of the reference can obviously align, so it must survive."""
    target = "".join(random.Random(1).choice("ACGT") for _ in range(300))
    idx = Prefilter([target], 16)
    for start in range(0, len(target) - 60, 37):
        frag = target[start:start + 60]
        assert idx.hits(frag, 1) >= 1, f"lost an exact fragment at {start}"


def test_the_reverse_complement_of_a_target_passes():
    """Reads come off both strands. Indexing one strand would drop ~half the library, and on a
    stranded paired library that is precisely the R2 mates."""
    target = "".join(random.Random(2).choice("ACGT") for _ in range(200))
    rc = target[::-1].translate(str.maketrans("ACGT", "TGCA"))
    idx = Prefilter([target], 16)
    assert idx.hits(rc[20:100], 1) >= 1


def test_unrelated_sequence_is_rejected():
    idx = Prefilter(["".join(random.Random(3).choice("ACGT") for _ in range(300))], 16)
    rnd = random.Random(4)
    junk = ["".join(rnd.choice("ACGT") for _ in range(100)) for _ in range(2000)]
    passed = sum(idx.mask(junk, 1, 1))
    assert passed / len(junk) < 0.05, f"{passed}/2000 random reads passed — no specificity"


def test_a_window_covering_an_N_is_dropped_not_guessed():
    """An N is unknown, not a wildcard. Guessing it would invent a hit; the design drops the
    window instead, which is why a read of pure Ns must score zero rather than match everything."""
    target = "ACGT" * 20
    idx = Prefilter([target], 16)
    assert idx.hits("N" * 100, 1) == 0
    # An N in the middle kills the k windows covering it, but the flanks still hit.
    assert idx.hits(target[:20] + "N" + target[20:40], 1) >= 1


def test_hits_stops_counting_at_min_hits():
    """`hits` is capped so the pass case exits early; it is a threshold test, not a census."""
    target = "".join(random.Random(5).choice("ACGT") for _ in range(500))
    idx = Prefilter([target], 16)
    assert idx.hits(target, 1) == 1
    assert idx.hits(target, 3) == 3


def test_threading_does_not_change_the_answer():
    """`mask` writes from several threads into one buffer. std::vector<bool> would be a bitfield
    and would race; this asserts the answer is thread-count invariant."""
    target = "".join(random.Random(6).choice("ACGT") for _ in range(400))
    idx = Prefilter([target], 16)
    rnd = random.Random(7)
    seqs = [target[i % 300:i % 300 + 80] if i % 3 == 0
            else "".join(rnd.choice("ACGT") for _ in range(80)) for i in range(5000)]
    one = idx.mask(seqs, 1, 1)
    assert one == idx.mask(seqs, 1, 8)
    assert one == idx.mask(seqs, 1, 32)


def test_empty_input_is_not_an_error():
    idx = Prefilter(["ACGT" * 10], 16)
    assert idx.mask([], 1, 4) == []


def test_k_outside_the_supported_range_is_rejected():
    """k < 12 would make the first-level bitset index wider than the k-mer itself."""
    with pytest.raises(ValueError):
        Prefilter(["ACGT" * 10], 8)
    with pytest.raises(ValueError):
        Prefilter(["ACGT" * 10], 33)


def test_keep_mask_indexes_the_fasta_mmseqs_searches(tmp_path):
    """The index must come from the search target. A prefilter built over V+J alone loses 16.29 %
    of real reads, 69.27 % of them J->C -- so what gets indexed is the whole ballgame, and this
    asserts the wiring reads the same FASTA rather than a hand-listed segment set."""
    fa = tmp_path / "alleles.fasta"
    target = "".join(random.Random(8).choice("ACGT") for _ in range(300))
    fa.write_text(f">t1\n{target[:150]}\n{target[150:]}\n")   # wrapped: the reader must join lines
    recs = [("hit", target[100:200]), ("miss", "A" * 100)]
    assert prefilter.keep_mask(recs, fa, threads=2) == [1, 0]
