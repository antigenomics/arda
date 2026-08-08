"""An equal-overlap tie in `_assign_coverage` must break on mismatches, not on arrival order.

A read is assigned to the clonotype whose junction it best covers. Overlap length alone does not
decide that: reads routinely span a junction end to end, so *ties are the normal case*, and the
losing root can be the one the read matches perfectly.

This is not hypothetical. In arda <= 2.9.0 a phantom clonotype -- a junction window slid 10 nt into
V framework by a target-inverted alignment -- competed with the true Jurkat clone, and ~47 % of the
5,758 reads it took were exact 48-vs-48 overlap ties on which the losing TRUE root matched with
**0** mismatches against the phantom's **1**. The mismatch count was computed one line below the
tie-break and thrown away. Measured on the delivered Jurkat run, fixing the tie-break alone returns
**2,510** reads to the true clone (15,448 -> 17,958) at an identical clonotype count.

No reference DB and no MMseqs2 -- `_assign_coverage` takes a DataFrame and a list of junctions.
"""

from __future__ import annotations

import polars as pl

from arda.rnaseq.correct import _assign_coverage

# A 48 nt junction and a decoy that differs from it at ONE base (position 30, A->C). A read equal to
# the true junction overlaps both roots by the full 48 nt -- an exact tie on overlap.
TRUE = "TGTGCCAGCAGTTTCTCGACCTGTTCGGCTAACTATGGCTACACCTTC"
DECOY = TRUE[:30] + ("C" if TRUE[30] != "C" else "G") + TRUE[31:]


def _reads(seqs: dict[str, str]) -> pl.DataFrame:
    n = len(seqs)
    return pl.DataFrame({
        "sequence_id": list(seqs),
        "sequence": list(seqs.values()),
        "locus": ["TRB"] * n,
        "v_call": ["TRBV12-3*01"] * n,
        "j_call": ["TRBJ1-2*01"] * n,
        "c_call": ["TRBC1"] * n,
        "rev_comp": ["F"] * n,
        "junction": [""] * n,
    })


def _assign(root_jn: list[str], seqs: dict[str, str]) -> dict[str, int]:
    got = _assign_coverage(_reads(seqs), root_jn, ["TRB"] * len(root_jn), {})
    return {jn: len(ids) for jn, ids in zip(root_jn, got)}


def test_a_perfect_match_beats_a_one_mismatch_root_at_equal_overlap():
    """The read matches TRUE exactly and DECOY with one mismatch; overlap is 48 for both."""
    read = {"r1": TRUE}

    # DECOY listed FIRST, so encounter order favours the wrong root.
    counts = _assign([DECOY, TRUE], read)

    assert counts[TRUE] == 1, (
        "the read matches this root with 0 mismatches and the other with 1, at identical 48 nt "
        "overlap; it was assigned by arrival order instead"
    )
    assert counts[DECOY] == 0


def test_the_result_does_not_depend_on_root_order():
    """Same tie, both orders -- the assignment must not move."""
    read = {"r1": TRUE}

    first = _assign([DECOY, TRUE], read)
    second = _assign([TRUE, DECOY], read)

    assert first[TRUE] == second[TRUE] == 1
    assert first[DECOY] == second[DECOY] == 0


def test_a_longer_overlap_still_wins_over_a_better_match():
    """Overlap remains the PRIMARY key -- mismatches only break ties.

    A root that the read covers longer is the better explanation even if a shorter one matches
    perfectly, which is what the coverage read-count is for.
    """
    short = TRUE[:24]          # read covers this perfectly, but only over 24 nt
    read = {"r1": TRUE}

    counts = _assign([short, DECOY], read)

    assert counts[DECOY] == 1, "the 48 nt overlap should beat a perfect 24 nt one"
    assert counts[short] == 0


def test_reads_are_assigned_once_each():
    """The invariant the coverage count depends on."""
    reads = {f"r{i}": TRUE for i in range(5)}

    counts = _assign([TRUE, DECOY], reads)

    assert sum(counts.values()) == 5
    assert counts[TRUE] == 5
