"""`_segment_rows` reduces in polars; the consuming loop must not be able to tell.

The segment search runs at the full `--max-seqs`, so it emits tens of rows per read and
`_segment_best_hits` skips nearly all of them. Reducing in polars instead of building a dict per
row is only safe if the reduction keeps EXACTLY the rows the loop can act on: the top row per
`(query, side)` plus up to `_MAX_TIED_V` exactly-tied V rows.

This is a differential test against a brute-force reimplementation of the loop's own filter, run
over a frame with the shapes that matter -- ties at the top, ties below the top, more ties than
the cap, unrecognised prefixes, and both sides for one read.
"""
from __future__ import annotations

import pytest

from arda.annotate.mapper import _MAX_TIED_V, _SEGMENT_FORMAT, _segment_rows


def _write(tmp_path, rows):
    p = tmp_path / "seg.tsv"
    p.write_text("".join("\t".join(str(c) for c in r) + "\n" for r in rows))
    return p


def _brute(rows):
    """Which rows the loop in `_segment_best_hits` can act on, decided row by row."""
    cols = _SEGMENT_FORMAT.split(",")
    recs = [dict(zip(cols, r)) for r in rows]
    recs = [dict(r, bits=float(r["bits"])) for r in recs]
    recs.sort(key=lambda r: (-r["bits"], r["target"]))
    top, tied, kept = {}, {}, []
    for r in recs:
        kind, sep, _name = r["target"].partition("|")
        if not sep or kind not in ("V", "J", "JC", "C"):
            continue
        side = kind if kind in ("V", "C") else "J"
        if (r["query"], side) in top:
            if side == "V" and r["bits"] == top[(r["query"], side)] \
                    and len(tied[r["query"]]) < _MAX_TIED_V:
                tied[r["query"]].append(r)
                kept.append(r)
            continue
        top[(r["query"], side)] = r["bits"]
        if side == "V":
            tied[r["query"]] = [r]
        kept.append(r)
    return kept


def _key(rows):
    return [(r["query"], r["target"], float(r["bits"])) for r in rows]


@pytest.mark.parametrize("rows", [
    # top row only, one side
    [("r1", "V|TRAV1*01", 100, 1, 50, 1)],
    # exact tie at the top: both kept
    [("r1", "V|TRAV1*01", 100, 1, 50, 1), ("r1", "V|TRAV1*02", 100, 1, 50, 1)],
    # a lower-scoring V is dropped, a tie is not
    [("r1", "V|TRAV1*01", 100, 1, 50, 1), ("r1", "V|TRAV1*02", 100, 1, 50, 1),
     ("r1", "V|TRAV2*01", 99, 1, 50, 1), ("r1", "V|TRAV3*01", 10, 1, 50, 1)],
    # more ties than the cap
    [("r1", f"V|TRAV{i}*01", 100, 1, 50, 1) for i in range(_MAX_TIED_V + 5)],
    # both sides for one read; J ties are NOT collected (only V ties are)
    [("r1", "V|TRAV1*01", 100, 1, 50, 1), ("r1", "J|TRAJ1*01", 80, 51, 90, 1),
     ("r1", "J|TRAJ2*01", 80, 51, 90, 1), ("r1", "JC|TRA_5", 70, 51, 99, 1)],
    # unrecognised prefix must be dropped, never treated as a J
    [("r1", "junk", 500, 1, 50, 1), ("r1", "V|TRAV1*01", 100, 1, 50, 1)],
    # `C|` is its OWN side, so a C row does not compete with the J row and neither hides the other.
    # This is the case the reduction silently dropped when `C|` targets were added: the filter kept
    # only V/J/JC, so `best_c` was always empty, no constant-only read was ever rescued, and 15
    # J->C reads vanished without `no_segment_hit` moving -- the rows were gone before any counter
    # saw them.
    [("r1", "C|IGHG1*01", 150, 20, 100, 1), ("r1", "J|IGHJ4*02", 60, 1, 19, 1)],
    # a C row that outscores everything must not displace the V or the J
    [("r1", "C|TRBC2*01", 300, 30, 100, 1), ("r1", "V|TRBV1*01", 100, 1, 29, 1),
     ("r1", "J|TRBJ1*01", 90, 25, 40, 1)],
    # only the best C survives; C ties are NOT collected (that is a V-only rule)
    [("r1", "C|IGHG1*01", 150, 1, 50, 1), ("r1", "C|IGHG3*01", 150, 1, 50, 1),
     ("r1", "C|IGHM*01", 80, 1, 50, 1)],
    # a constant-region-only read: C is its only evidence, and it must survive the reduction or
    # the read never enters `seen` and is never rescued
    [("r1", "C|IGKC*01", 170, 1, 90, 1)],
    # pre-split `JC|` and post-split `C|` in one frame -- a mapper of this vintage must read a
    # reference of either, and JC stays on the J side while C does not
    [("r1", "JC|IGH_12", 200, 1, 90, 1), ("r1", "C|IGHG1*01", 150, 20, 100, 1),
     ("r1", "V|IGHV1*01", 90, 1, 19, 1)],
    # two reads interleaved by score, so the global sort splits each read's rows apart
    [("r1", "V|TRAV1*01", 100, 1, 50, 1), ("r2", "V|TRBV1*01", 99, 1, 50, 1),
     ("r1", "V|TRAV2*01", 98, 1, 50, 1), ("r2", "V|TRBV1*02", 99, 1, 50, 1),
     ("r2", "J|TRBJ1*01", 60, 51, 90, 1), ("r1", "J|TRAJ1*01", 61, 51, 90, 1)],
])
def test_reduction_keeps_exactly_the_actionable_rows(tmp_path, rows):
    got = _segment_rows(_write(tmp_path, rows))
    assert _key(got) == _key(_brute(rows))


def test_reduction_is_invariant_to_input_row_order(tmp_path):
    """The determinism guarantee: the same rows in a different file order give the same answer."""
    rows = [("r1", "V|TRAV1*01", 100, 1, 50, 1), ("r1", "V|TRAV1*02", 100, 1, 50, 1),
            ("r1", "V|TRAV2*01", 100, 1, 50, 1), ("r1", "J|TRAJ1*01", 80, 51, 90, 1),
            ("r2", "V|TRBV1*01", 100, 1, 50, 1)]
    a = _key(_segment_rows(_write(tmp_path, rows)))
    b = _key(_segment_rows(_write(tmp_path, list(reversed(rows)))))
    assert a == b


def test_empty_and_missing(tmp_path):
    assert _segment_rows(tmp_path / "nope.tsv") == []
    assert _segment_rows(_write(tmp_path, [])) == []


def test_the_side_mapping_is_the_one_the_two_pass_depends_on():
    """`_SEGMENT_SIDE` is shared by the reduction and the loop, so it is the whole contract.

    Spelling this rule out twice -- once in polars, once in Python -- is what let `C|` targets be
    added to the loop but not the filter: the reduction discarded every C row, `best_c` stayed
    empty, no constant-only read was rescued, and 15 J->C reads vanished **without
    `no_segment_hit` moving**, because the rows were gone before any counter saw them.
    """
    from arda.annotate.mapper import _SEGMENT_SIDE

    assert _SEGMENT_SIDE["C"] == "C", (
        "a constant-region hit says what the isotype is and NOTHING about which J the read "
        "carries; folding it into the J side is what the old JC| targets did wrong")
    assert _SEGMENT_SIDE["JC"] == "J", (
        "JC| is the pre-split kind and must stay J-side, or this mapper mis-reads a reference "
        "built before the constant region became its own target")
    assert _SEGMENT_SIDE["V"] == "V" and _SEGMENT_SIDE["J"] == "J"
    assert set(_SEGMENT_SIDE.values()) == {"V", "J", "C"}
