"""Unit tests for the C++ markup-transfer hot path.

These encode the hand-verified projection cases (full / insertion / deletion /
truncation) so regressions in the alignment walk are caught immediately.
"""

from arda import _markup


def tr(qaln, taln, qs, ts, starts, ends):
    return _markup.transfer_regions(qaln, taln, qs, ts, starts, ends)


def test_full_identity():
    qaln = taln = "ABCDEFGHIJ"
    assert tr(qaln, taln, 1, 1, [1, 4, 7], [3, 6, 10]) == [(1, 3), (4, 6), (7, 10)]


def test_insertion_in_query_shifts_downstream():
    # 2-base insertion (gap in target) before region 2
    qaln = "ABCxxDEFGHIJ"
    taln = "ABC--DEFGHIJ"
    assert tr(qaln, taln, 1, 1, [1, 4, 7], [3, 6, 10]) == [(1, 3), (6, 8), (9, 12)]


def test_deletion_in_query():
    # target D,E deleted in query
    qaln = "ABC--FGHIJ"
    taln = "ABCDEFGHIJ"
    assert tr(qaln, taln, 1, 1, [1, 4, 7], [3, 6, 10]) == [(1, 3), (4, 4), (5, 8)]


def test_truncated_query_region_uncovered():
    # alignment starts at target position 4 -> region 1 (target 1..3) uncovered
    qaln = taln = "DEFGHIJ"
    assert tr(qaln, taln, 1, 4, [1, 4, 7], [3, 6, 10]) == [(-1, -1), (1, 3), (4, 7)]


def test_project_region_primitive():
    # 0-based inclusive primitive used internally
    taln = "ABCDEFG"
    assert _markup.project_region(taln, taln, 0, 0, 2, 4) == (2, 4)
    qaln = "ABCD-FG"
    assert _markup.project_region(qaln, taln, 0, 0, 2, 4) == (2, 3)
