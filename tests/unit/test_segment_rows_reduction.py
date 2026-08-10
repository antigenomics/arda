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

from arda.annotate.mapper import _MAX_TIED, _MAX_TIED_J, _MAX_TIED_V, _SEGMENT_FORMAT, _segment_rows


def _write(tmp_path, rows):
    p = tmp_path / "seg.tsv"
    # Rows are written in `_SEGMENT_FORMAT` order. A case that does not care about the target span
    # gives the first six fields and gets a FORWARD `tend`, so only the inverted-row cases have to
    # spell one out. `_brute` never reads `tend`, and `zip` stops at the shorter of the two, so it
    # takes the unpadded rows unchanged.
    rows = [r if len(r) == 7 else (*r, r[5] + 40) for r in rows]
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
            # A side absent from `_MAX_TIED` (C) collects no ties at all: its cap is effectively 1,
            # which the `_rank == 0` clause already covers.
            if side in _MAX_TIED and r["bits"] == top[(r["query"], side)] \
                    and len(tied[(r["query"], side)]) < _MAX_TIED[side]:
                tied[(r["query"], side)].append(r)
                kept.append(r)
            continue
        top[(r["query"], side)] = r["bits"]
        if side in _MAX_TIED:
            tied[(r["query"], side)] = [r]
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
    # both sides for one read; J ties ARE collected now -- an exact tie between two `J|` targets
    # used to be broken lexicographically, and that decided which V×J scaffold the read was
    # aligned against. Measured: SRR5233639.3589/2 ties IGLJ2*01,IGLJ3*01 and IGLJ2A*01 at 54 bits,
    # the comma sorts before `A`, and the read landed on a scaffold scoring 93 instead of 96.
    [("r1", "V|TRAV1*01", 100, 1, 50, 1), ("r1", "J|TRAJ1*01", 80, 51, 90, 1),
     ("r1", "J|TRAJ2*01", 80, 51, 90, 1), ("r1", "JC|TRA_5", 70, 51, 99, 1)],
    # more J ties than the J cap (which is smaller than the V cap: J targets are ~40-60 nt)
    [("r1", f"J|TRAJ{i}*01", 80, 51, 90, 1) for i in range(_MAX_TIED_J + 3)],
    # a J tie one bit below the top is NOT a tie and must be dropped
    [("r1", "J|TRAJ1*01", 80, 51, 90, 1), ("r1", "J|TRAJ2*01", 79, 51, 90, 1)],
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


def test_a_target_inverted_segment_row_is_refused(tmp_path):
    """`tstart > tend` must not survive the reduction, on this path as well as `_best_hits`.

    A segment row's `tstart` is read as a FORWARD target offset by both of its consumers, and
    neither can tell that it is not one. `_align_implied` builds its prefilter diagonal from
    `qstart - tstart`, so an inverted row demotes the read to the full-reference rescue (time, not
    reads); `_cannot_reach_cys104` walks `tstart + |qend - qstart|` forward against the segment's
    cdr3 start, so it can answer "unreachable" for a read that does reach Cys104 -- and that read
    is then aligned V-segment-only and emits no junction at all.
    """
    inverted = ("r1", "V|TRAV1*01", 100, 1, 50, 180, 1)
    assert _segment_rows(_write(tmp_path, [inverted])) == []


def test_an_inverted_segment_row_never_displaces_a_forward_one(tmp_path):
    """Refusal must survive the sort: outscoring a real hit is the case where keeping it does the
    damage, because the reduction keeps the top row per `(query, side)`."""
    good = ("r1", "V|TRAV1*01", 100, 1, 50, 1, 50)
    inverted = ("r1", "V|TRAV2*01", 243, 1, 50, 180, 1)

    for name, rows in (("inv_first", [inverted, good]), ("good_first", [good, inverted])):
        got = _segment_rows(_write(tmp_path, rows))
        assert [r["target"] for r in got] == ["V|TRAV1*01"], name


def test_the_reduction_returns_the_keys_segmap_builds(tmp_path):
    """`tend` is asked of mmseqs only to refuse the inverted row and must not leak downstream:
    `segmap.segment_rows` builds its dicts by hand and has no such key, and the two feed the same
    consuming loop."""
    got = _segment_rows(_write(tmp_path, [("r1", "V|TRAV1*01", 100, 1, 50, 1)]))
    assert "tend" not in got[0]
    assert set(got[0]) == {"query", "target", "bits", "qstart", "qend", "tstart"}


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


def test_a_pre_split_segments_fasta_is_detected_by_format(tmp_path):
    """`_has_jc_targets` is the self-healing gate for an upgraded install.

    Upgrading arda does not rewrite a generated artifact, and `segments.fasta` changed shape in
    2.8.0 (345 `JC|` scaffolds -> 25 `C|` targets). The mapper still READS `JC|`, so a stale file
    produces correct output with no error and none of the speedup -- invisible forever. Detecting
    it by format rather than by mtime or version is what makes the repair automatic.
    """
    from arda.annotate.mapper import _has_jc_targets

    old = tmp_path / "old.fasta"
    old.write_text(">V|TRAV1*01\nACGT\n>J|TRAJ1*01\nACGT\n>JC|TRA_JC_0\nACGTACGT\n")
    assert _has_jc_targets(old)

    new = tmp_path / "new.fasta"
    new.write_text(">V|TRAV1*01\nACGT\n>J|TRAJ1*01\nACGT\n>C|TRAC*01\nACGTACGT\n")
    assert not _has_jc_targets(new)

    # `C|` must not be mistaken for `JC|` in either direction, and a sequence line that happens to
    # start with the marker text is not a header.
    tricky = tmp_path / "tricky.fasta"
    tricky.write_text(">C|IGHG1*01\nACGT\n>V|IGHV1*01\nACGT\n")
    assert not _has_jc_targets(tricky)
