"""Chained (two-diagonal) scoring in `segmap`.

An ungapped extension follows ONE diagonal. On a ``V + pad + J`` scaffold a junction-spanning read
therefore lands on two diagonals -- the V half and the J half -- and `collect` already extends both.
The shipped rule then takes the *better* of the two, which scores the read ``max(V, J)`` rather than
``V + J``.

That is not a tie-break detail. Every scaffold sharing the read's J scores *identically* under the
max rule, so the winner among ~50-70 V partners of that J is decided by target order rather than by
evidence. Measured on the real 15,414-scaffold reference (arda-benchmark round 17, SRR5233635): the
true V's best scaffold ranks strictly below a wrong-V scaffold on 1,361 of 2,000 reads (68.05 %),
while the true scaffold is always present in the candidate pool (recall@250 = 1.0000). That is the
``v_gene`` .3430 this project recorded as a *structural* failure of the mapper; it is a scoring rule.

These tests need no reference DB and no MMseqs2 -- the geometry is expressible directly.
"""

from __future__ import annotations

import random

from arda._segmap import SegmentMapper

PAD = "N" * 9          # the N-pad a real V x J scaffold carries between its two segments
K = 12


def _rnd(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _scaffolds(rng: random.Random, n_decoy: int = 3):
    """``(targets, true_index, read)``: decoys sharing the true scaffold's J, true one LAST.

    The read carries a SHORT V tail (the weak, discriminating half) and a LONG J (the strong half,
    identical across every target). Placing the true scaffold last means the lowest-target-index
    tie-break cannot rescue it -- exactly the situation on the real reference.
    """
    j = _rnd(60, rng)
    vs = [_rnd(120, rng) for _ in range(n_decoy + 1)]
    targets = [v + PAD + j for v in vs]
    true_i = len(targets) - 1
    read = vs[true_i][-30:] + _rnd(20, rng) + j[:60]
    return targets, true_i, read


def _winner(mapper, read, chain_offset):
    hits = mapper.map([read], max_tied=8, min_score=40, threads=1, chain_offset=chain_offset)
    if not hits:
        return None, []
    scored = sorted(((h[1], h[2]) for h in hits), key=lambda x: -x[1])
    return scored[0][0], scored


def test_max_over_diagonals_picks_a_wrong_scaffold_that_shares_the_j():
    """The defect, pinned. Without chaining the true scaffold does not win."""
    rng = random.Random(1)
    targets, true_i, read = _scaffolds(rng)
    mapper = SegmentMapper(targets, [0] * len(targets), K)

    win, scored = _winner(mapper, read, chain_offset=0)

    assert win != true_i, (
        "the shipped max-over-diagonals rule unexpectedly picked the true scaffold; this fixture "
        "is supposed to reproduce the ranking failure that chaining fixes"
    )
    # Every J-sharing scaffold ties, which is why the choice is arbitrary.
    assert len({s for _, s in scored}) == 1, f"expected one tied score band, got {scored}"


def test_chaining_two_diagonals_recovers_the_true_scaffold():
    """Summing the two best compatible diagonals ranks the true scaffold first."""
    rng = random.Random(1)
    targets, true_i, read = _scaffolds(rng)
    mapper = SegmentMapper(targets, [0] * len(targets), K)

    win, scored = _winner(mapper, read, chain_offset=60)

    assert win == true_i, f"chained scoring picked {win}, expected the true scaffold {true_i}"
    # The true scaffold is the only one that can pair two diagonals, so it stands alone above the
    # decoys rather than merely winning a tie.
    assert scored[0][1] > (scored[1][1] if len(scored) > 1 else 0)


def test_chaining_sums_the_two_halves():
    """The chained score is the sum of both halves, not the better one."""
    rng = random.Random(0)
    v, j = _rnd(120, rng), _rnd(60, rng)
    mapper = SegmentMapper([v + PAD + j], [0], K)
    read = v[-60:] + _rnd(20, rng) + j[:40]

    (_, plain), = [(h[1], h[2]) for h in
                   mapper.map([read], max_tied=8, min_score=40, threads=1, chain_offset=0)]
    (_, chained), = [(h[1], h[2]) for h in
                     mapper.map([read], max_tied=8, min_score=40, threads=1, chain_offset=60)]

    assert chained > plain, f"chained {chained} should exceed the single-diagonal {plain}"


def test_chain_offset_zero_is_the_shipped_behaviour():
    """The flag is off by default, and off must be byte-identical to the old rule."""
    rng = random.Random(7)
    targets, _, read = _scaffolds(rng)
    mapper = SegmentMapper(targets, [0] * len(targets), K)

    default = mapper.map([read], max_tied=8, min_score=40, threads=1)
    explicit = mapper.map([read], max_tied=8, min_score=40, threads=1, chain_offset=0)

    assert default == explicit


def test_a_distant_second_diagonal_is_not_chained():
    """Two diagonals far apart are a repeat, not one rearrangement -- do not sum them."""
    rng = random.Random(3)
    targets, true_i, read = _scaffolds(rng)
    mapper = SegmentMapper(targets, [0] * len(targets), K)

    near, _ = _winner(mapper, read, chain_offset=60)
    far, _ = _winner(mapper, read, chain_offset=1)   # too tight to pair the two halves

    assert near == true_i
    assert far != true_i, "an offset of 1 nt must not be able to pair a 20 nt NDN insert"


def test_coordinates_stay_on_one_diagonal():
    """`qstart`/`tstart` must remain a matched pair.

    ``annotate.project`` places the junction as ``qstart + (anchor_nt + 1 - tstart)``. If chaining
    merged the span of two diagonals that arithmetic would break silently -- the class of defect
    that shipped a well-formed wrong junction in 2.5.6.
    """
    rng = random.Random(0)
    v, j = _rnd(120, rng), _rnd(60, rng)
    mapper = SegmentMapper([v + PAD + j], [0], K)
    read = v[-60:] + _rnd(20, rng) + j[:40]

    (plain,) = mapper.map([read], max_tied=8, min_score=40, threads=1, chain_offset=0)
    (chained,) = mapper.map([read], max_tied=8, min_score=40, threads=1, chain_offset=60)

    # index 3,4,5 = qstart, qend, tstart -- unchanged; only the score (index 2) moves.
    assert chained[3:6] == plain[3:6], (
        f"chaining moved the coordinates {plain[3:6]} -> {chained[3:6]}; the junction projection "
        "depends on qstart/tstart being a matched pair on ONE diagonal"
    )
    assert chained[2] > plain[2]
