"""The structure-aware segment mapper: seed, vote by diagonal, extend ungapped.

DB-free and deterministic — every case builds its own tiny reference, so these run in CI.

The component exists because the segment pass asks a general homology engine a question that has a
structural answer: which V and which J, with coordinates, against a fixed 236 kb germline reference.
Measured on 100,000 amplicon reads against the shipped 924-target reference, 8 threads:

    mmseqs search   2770 ms      segmap   74 ms      37x

with V allele agreement .9997 and J allele .9998 against mmseqs' own best-per-side.
"""
from __future__ import annotations

import pytest

_segmap = pytest.importorskip("arda._segmap", reason="native _segmap extension not built")


def rc(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


V1 = "ACGTTGCAACGTTGCAACGTTGCAGGGGTTTTCCCCAAAA"          # 40 nt
V2 = "ACGTTGCAACGTTGCAACGTTGCAGGGGTTTTCCCCAAAT"          # V1 with the last base changed
J1 = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG"          # 40 nt, unrelated to V1


def mapper(seqs, groups, k=16):
    return _segmap.SegmentMapper(seqs, groups, k)


def test_finds_the_target_a_read_came_from():
    sm = mapper([V1, J1], [0, 1])
    (hit,) = sm.map([V1], 8, 40, 1)
    qi, ti, score, qstart, qend, tstart, is_rc, split = hit
    assert (qi, ti, is_rc) == (0, 0, 0)
    assert (qstart, qend, tstart) == (1, len(V1), 1)
    assert score == 2 * len(V1)
    assert split == 0          # indel detection is opt-in; `max_indel` defaults to 0


def test_reverse_complement_keeps_the_mmseqs_strand_convention():
    """A minus-strand hit is signalled by ``qstart > qend``.

    `_align_implied` reads exactly that to decide which reads need reverse-complementing before the
    second alignment, so the convention is load-bearing, not cosmetic.
    """
    sm = mapper([V1], [0])
    (fwd,) = sm.map([V1], 8, 40, 1)
    (rev,) = sm.map([rc(V1)], 8, 40, 1)
    assert fwd[6] == 0 and rev[6] == 1
    assert fwd[3] < fwd[4], "forward hit must have qstart < qend"
    assert rev[3] > rev[4], "minus-strand hit must have qstart > qend"
    assert fwd[2] == rev[2], "the same alignment, read either way, scores the same"


def test_one_best_per_group_plus_exact_ties():
    """Best hit per (read, group), and exactly-tied targets up to the cap — the rule the consuming
    loop in `_segment_best_hits` applies, mirrored here so it need not re-derive it."""
    sm = mapper([V1, V2, J1], [0, 0, 1])
    rows = sm.map([V1], 8, 40, 1)
    by_group = {}
    for _, ti, score, *_ in rows:
        by_group.setdefault(0 if ti in (0, 1) else 1, []).append((ti, score))
    # V1 beats V2 (one mismatch), so no tie: exactly one V row.
    assert len(by_group[0]) == 1 and by_group[0][0][0] == 0

    # An exact tie must return both, in target order.
    sm2 = mapper([V1, V1, J1], [0, 0, 1])
    rows2 = [r for r in sm2.map([V1], 8, 40, 1) if r[1] in (0, 1)]
    assert [r[1] for r in rows2] == [0, 1], "tied targets are returned in reference order"
    assert rows2[0][2] == rows2[1][2]


def test_max_tied_caps_the_expansion():
    sm = mapper([V1] * 6 + [J1], [0] * 6 + [1])
    assert len([r for r in sm.map([V1], 2, 40, 1) if r[1] < 6]) == 2
    assert len([r for r in sm.map([V1], 4, 40, 1) if r[1] < 6]) == 4


def test_min_score_is_the_significance_floor():
    """Without one, every seeded diagonal with a positive extension is reported.

    Calibrated against mmseqs on real data: at no floor, 43,010 reads get a constant-region hit
    against mmseqs' 473, and half of those spurious hits score exactly 38 — a bare 16-mer seed
    (16 x MATCH = 32) plus a couple of flanking matches. 40 is the first value above that, and it
    reproduces mmseqs to within 0.25 % on V and 0.09 % on J.
    """
    # The target and the read share exactly one 16-mer and then diverge into runs that cannot
    # match each other anywhere — poly-C against poly-A. A fixture reusing V1 does NOT work: V1
    # carries a GGGG run, so a G-tail extends the seed by four and scores 40.
    seed = "AAACCCGGGTTTACGT"                      # not its own reverse complement
    target = seed + "C" * 20
    seed_only = seed + "A" * 20
    sm = mapper([target], [0])
    (hit,) = sm.map([seed_only], 8, 0, 1)
    assert hit[2] == 2 * len(seed), "a bare seed scores exactly its own length"
    assert not sm.map([seed_only], 8, 40, 1), "the floor must reject a seed with no extension"


def test_a_read_matching_nothing_returns_nothing():
    sm = mapper([V1], [0])
    assert sm.map(["A" * 60], 8, 40, 1) == []
    assert sm.map([""], 8, 40, 1) == []
    assert sm.map([], 8, 40, 1) == []


def test_n_bases_break_the_seed_and_are_never_guessed():
    sm = mapper([V1], [0])
    assert not sm.map(["N" * len(V1)], 8, 40, 1)
    # An N in the middle destroys every window covering it, but the flanks can still seed.
    holed = V1[:20] + "N" + V1[21:]
    hits = sm.map([holed], 8, 40, 1)
    assert hits and hits[0][1] == 0


def test_threading_does_not_change_the_answer():
    """`map` is const and every worker shares one mapper, so per-call state must not live on the
    object. Running the same input at 1 and 8 threads is the cheapest check that it does not."""
    reads = [V1, rc(V1), J1, "A" * 50, V2] * 40
    sm = mapper([V1, V2, J1], [0, 0, 1])
    assert sm.map(reads, 8, 40, 1) == sm.map(reads, 8, 40, 8)


def test_index_shape_is_reported():
    sm = mapper([V1, J1], [0, 1])
    assert sm.k == 16
    assert sm.n_groups == 2
    assert sm.size > 0 and sm.postings >= sm.size


@pytest.mark.parametrize("bad", [
    ([V1], [0, 1]),          # groups longer than seqs
    ([V1, J1], [0]),         # groups shorter than seqs
])
def test_mismatched_groups_raise(bad):
    with pytest.raises(Exception):
        _segmap.SegmentMapper(bad[0], bad[1], 16)


def test_k_is_bounded():
    for k in (4, 40):
        with pytest.raises(Exception):
            _segmap.SegmentMapper([V1], [0], k)


# ---------------------------------------------------------------- wiring into the mapper

def test_segmap_wiring_matches_the_mmseqs_path_it_replaces():
    """The two constants that decide whether this path is equivalent, pinned to their counterparts.

    Both were wrong at first and both were caught by measurement, not by review:

    * **k.** 16 was inherited from `prefilter`, where it is calibrated for REJECTION. Seed length
      sets sensitivity to mismatches, and at k=16 the mapper seeded 53,048 reads against mmseqs'
      53,121 — those 73 reads are never rescued, because a read with no segment hit is assumed
      non-receptor. At k=12, the `-k` arda actually passes MMseqs2, `no_segment_hit` matches it
      exactly and the end-to-end AIRR delta fell from 30 lost reads to 6.
    * **the side mapping.** `C` must be its own group: a constant-region hit says what the isotype
      is and nothing about which J the read carries, and folding it into the J side is precisely
      what the pre-2.8.0 `JC|` targets did wrong.
    """
    from arda import segmap
    from arda.annotate.mapper import _KMER, _MAX_TIED, _SEGMENT_SIDE

    assert segmap.K == _KMER["nt"], (
        "the segment mapper must seed at the same k arda passes MMseqs2, or it is less sensitive "
        "than the search it replaces and silently drops reads that are never rescued")
    assert segmap.SIDE_GROUP.keys() == _SEGMENT_SIDE.keys(), (
        "segmap and the consuming loop disagree about which target kinds exist")
    # Same partition, whatever the group labels are called on each side.
    pairs = {(segmap.SIDE_GROUP[k], _SEGMENT_SIDE[k]) for k in _SEGMENT_SIDE}
    assert len({g for g, _ in pairs}) == len({s for _, s in pairs}) == 3
    assert len(pairs) == 3, "the two side mappings do not induce the same partition"
    assert max(_MAX_TIED.values()) >= _MAX_TIED["V"]


# --- indel detection (`max_indel` > 0) -------------------------------------------------------
#
# An ungapped extension follows ONE diagonal, so a read carrying an indel relative to germline
# scores only up to the indel. That is not a corner case: measured on 341,294 real IGH mates
# (benchmark round 12), 3.18 % of reads carry a V indel and the rate tracks SHM load --
# 0.74 % at >=98 % V identity, 8.00 % below 90 % -- because AID makes indels, not only
# substitutions. Those reads are REROUTED to the gapped rescue path, never dropped.

# 120 nt, germline-like, and every one of its 105 16-mers is distinct -- so any second diagonal
# these tests see comes from the engineered indel and not from an internal repeat.
IND_T = ("CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCC"
         "TCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCT")[:120]


def test_the_16mers_of_the_indel_fixture_are_distinct():
    """Guards the fixture itself: a repeat would make every case below pass for the wrong reason."""
    kmers = [IND_T[i:i + 16] for i in range(len(IND_T) - 15)]
    assert len(set(kmers)) == len(kmers)


def test_an_indel_halves_the_ungapped_score_which_is_why_this_exists():
    """The motivating measurement, in one assertion.

    A clean read scores 2 x len. Delete a single base in the middle and the SAME read scores
    half, because one diagonal covers only one half of it. Nothing is wrong with the extension --
    it is the right answer to the wrong question, and it is why an indel-bearing read must not
    have its scaffold decided on the fast path.
    """
    sm = mapper([IND_T], [0])
    (clean,) = sm.map([IND_T], 8, 40, 1, 0)
    (deleted,) = sm.map([IND_T[:60] + IND_T[61:]], 8, 40, 1, 0)
    assert clean[2] == 2 * len(IND_T)
    assert deleted[2] == pytest.approx(clean[2] / 2, abs=4)


@pytest.mark.parametrize("read,label", [
    (IND_T[:60] + IND_T[61:], "1 nt deletion"),
    (IND_T[:60] + IND_T[62:], "2 nt deletion"),
    (IND_T[:60] + "GGGG" + IND_T[60:], "4 nt insertion"),
])
def test_indel_reads_are_flagged_when_detection_is_on(read, label):
    sm = mapper([IND_T], [0])
    assert all(h[7] == 1 for h in sm.map([read], 8, 40, 1, 30)), label


def test_a_clean_read_is_never_flagged():
    sm = mapper([IND_T], [0])
    assert all(h[7] == 0 for h in sm.map([IND_T], 8, 40, 1, 30))


def test_detection_is_opt_in_and_off_by_default():
    """`max_indel=0` must leave `split` at 0 even for a read that plainly carries an indel.

    The flag is a routing decision with a speed cost, so the default has to be inert -- and the
    default is expressed in the C++ signature, not only in the Python wrapper.
    """
    sm = mapper([IND_T], [0])
    read = IND_T[:60] + IND_T[61:]
    assert all(h[7] == 0 for h in sm.map([read], 8, 40, 1, 0))
    assert all(h[7] == 0 for h in sm.map([read], 8, 40, 1))       # C++ default for max_indel


def test_diagonals_further_apart_than_max_indel_are_not_an_indel():
    """Two well-supported diagonals 40 nt apart are a rearrangement or a repeat, not one indel.

    Routing those to a gapped realignment would cost time and tell us nothing, so the bound is
    load-bearing rather than cosmetic -- and it is a BOUND, so raising it must flag the same read.
    """
    sm = mapper([IND_T], [0])
    far = IND_T[:60] + IND_T[100:]
    assert all(h[7] == 0 for h in sm.map([far], 8, 40, 1, 30))
    assert any(h[7] == 1 for h in sm.map([far], 8, 40, 1, 50))


def test_a_minus_strand_indel_read_is_flagged_and_keeps_the_strand_convention():
    """Detection runs per strand, so it must survive the reverse-complement path.

    `_align_implied` reads `qstart > qend` to know a read needs reverse-complementing; flagging
    must not disturb that.
    """
    sm = mapper([IND_T], [0])
    (hit,) = sm.map([rc(IND_T[:60] + IND_T[61:])], 8, 40, 1, 30)
    assert hit[7] == 1
    assert hit[3] > hit[4]


def test_segment_rows_exposes_split_and_defaults_to_not_looking(tmp_path):
    """The Python wrapper must carry the flag through and keep the same default as the extension."""
    from arda import segmap
    fa = tmp_path / "segments.fasta"
    fa.write_text(f">V|IGHV1-1*01\n{IND_T}\n")
    reads = {"r1": IND_T[:60] + IND_T[61:]}
    assert all(r["split"] is False for r in segmap.segment_rows(fa, reads, max_tied=8))
    flagged = segmap.segment_rows(fa, reads, max_tied=8, max_indel=segmap.MAX_INDEL_NT)
    assert flagged and all(r["split"] is True for r in flagged)
