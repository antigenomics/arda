"""The coverage assignment's inner loop, bounded — and why bounding it cannot move a read.

`correct._assign_coverage` decides which clonotype every partial read belongs to. Its inner scan
compared an integer mismatch count against a FLOAT budget (`max_mm * ov`) one base at a time in
Python, over a measured 3.55 M candidate diagonals per 20,000 reads.

⛔ Unlike the assembler's overlap test, this caller NEEDS the count and not a verdict: a read joins
the root with the longest overlap and, **on a tie, the fewer mismatches**. That tie-break is
load-bearing — when a phantom clonotype competed with the true Jurkat clone, ~47 % of the 5,758
reads it stole were exact 48-vs-48 overlap ties on which the losing TRUE root matched with 0
mismatches against the phantom's 1. So the bounded version must return the *true* count whenever
the row is accepted, and only collapse to "too many" when the row is rejected anyway.
"""
from __future__ import annotations

import pytest

from arda import _markup


def _reference(a: str, a_off: int, b: str, b_off: int, n: int, budget: float) -> tuple[int, bool]:
    """The Python this replaced: count one base at a time, break once past the float budget."""
    mm = 0
    for k in range(n):
        if a[a_off + k] != b[b_off + k]:
            mm += 1
            if mm > budget:
                break
    return mm, mm <= budget


@pytest.mark.parametrize("a,b,n", [
    ("ACGTACGTACGT", "ACGTACGTACGT", 12),
    ("ACGTACGTACGT", "ACGTACGTACGA", 12),
    ("ACGTACGTACGT", "TGCATGCATGCA", 12),
    ("ACGTACGTACGT", "ACGTACGTACGT", 0),
    ("A" * 48, "A" * 48, 48),
    ("A" * 48, "A" * 24 + "C" * 24, 48),
])
@pytest.mark.parametrize("max_mm_frac", [0.0, 0.02, 0.12, 0.5, 1.0])
def test_it_accepts_exactly_what_the_python_scan_accepted(a, b, n, max_mm_frac):
    """Same accept/reject, and the same COUNT on every accepted row."""
    budget = max_mm_frac * n
    want_mm, want_ok = _reference(a, 0, b, 0, n, budget)
    got = _markup.count_mismatches(a, 0, b, 0, n, int(budget))

    assert (got <= budget) is want_ok, f"accept/reject moved: got {got} vs budget {budget}"
    if want_ok:
        assert got == want_mm, "an ACCEPTED row reported the wrong count — the tie-break reads this"


def test_the_float_budget_collapses_to_floor():
    """⛔ The equivalence the speedup rests on. `max_mm * ov` is a float; the count is an integer,
    so `mm > budget` is exactly `mm > floor(budget)`. This is the real operating point:
    `max_mm=0.12` over a 48 nt junction is 5.76, i.e. 5 mismatches pass and 6 do not."""
    ov, budget = 48, 0.12 * 48                      # 5.76
    assert int(budget) == 5

    for count in range(0, 9):
        a = "C" * count + "A" * (ov - count)
        got = _markup.count_mismatches(a, 0, "A" * ov, 0, ov, int(budget))
        assert (got <= budget) is (count <= budget), f"{count} mismatches vs budget {budget}"
        if count <= budget:
            assert got == count


def test_offsets_address_the_diagonal():
    """The caller compares `s[lo:hi]` against `jr[lo+d:hi+d]` — a read laid on a diagonal against a
    junction — so both offsets are independent and neither string starts at 0."""
    read = "TTTT" + "ACGTACGT"
    junc = "GG" + "ACGTACGT" + "CCCC"
    assert _markup.count_mismatches(read, 4, junc, 2, 8, 0) == 0
    assert _markup.count_mismatches(read, 0, junc, 2, 8, 0) == 1   # exceeds -> max_mm+1


def test_it_never_reads_past_either_string():
    """`n` is clamped to what both strings actually have. A junction shorter than the projected
    overlap must not read out of bounds — in Python that raised IndexError, which would have
    surfaced; in C++ it would silently read adjacent memory."""
    assert _markup.count_mismatches("ACGT", 0, "ACGT", 0, 999, 0) == 0
    assert _markup.count_mismatches("ACGT", 0, "AC", 0, 4, 0) == 0
    assert _markup.count_mismatches("ACGT", 9, "ACGT", 0, 4, 0) == 1     # offset past the end
    assert _markup.count_mismatches("", 0, "", 0, 0, 0) == 0


def test_a_rejected_row_reports_more_than_the_budget():
    """The contract the caller relies on: when the scan gives up it must return something the
    `mm <= budget` test rejects, never a truncated count that could sneak past."""
    for max_mm in (0, 1, 5, 20):
        got = _markup.count_mismatches("A" * 100, 0, "C" * 100, 0, 100, max_mm)
        assert got == max_mm + 1
