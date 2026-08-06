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
        if not sep or kind not in ("V", "J", "JC"):
            continue
        side = "V" if kind == "V" else "J"
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
