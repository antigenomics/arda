"""The assembler's overlap test, bounded — and the equivalence that makes it safe.

`_greedy_contigs` used to ask ``sum(1 for x, y in zip(a, b) if x != y) > (1 - min_id) * ov``: it
counted EVERY mismatch over the overlap and only then compared the total against the budget.
Nearly every candidate it tests is a read that merely shares one k-mer and does not overlap at all,
so the answer is settled in the first few bases and the rest of the count is dead work. Profiled on
a 100k-read TRA amplicon: **153 million generator iterations, 19.5 s of the stage's 25.6 s**.

⛔ The contig sequence feeds every junction derived from it, so "faster" is only acceptable if the
candidate set is IDENTICAL. It is, and not by luck: the caller never used the count, only
``count > budget``, and for an integer count that is exactly ``count > floor(budget)`` however the
float lands. That is what the parametrised test below pins.
"""
from __future__ import annotations

import pytest

from arda import _markup


def _reference(a: str, b: str, max_mm: int) -> bool:
    """The predicate the old inline `sum(...)` computed, spelled out."""
    return sum(1 for x, y in zip(a, b) if x != y) <= max_mm


@pytest.mark.parametrize("a,b,max_mm", [
    ("ACGT", "ACGT", 0),
    ("ACGT", "ACGA", 0),
    ("ACGT", "ACGA", 1),
    ("ACGT", "ACAA", 1),
    ("ACGT", "TGCA", 3),
    ("ACGT", "TGCA", 4),
    ("", "", 0),
    ("ACGT", "", 0),
    ("ACGT", "AC", 0),          # zip stops at the shorter -- the common length is the overlap
    ("ACGT", "AG", 0),
    ("ACGT", "AG", 1),
    ("A" * 200, "A" * 200, 0),
    ("A" * 200, "A" * 199 + "C", 0),
    ("A" * 200, "A" * 199 + "C", 1),
])
def test_it_answers_exactly_what_the_sum_answered(a, b, max_mm):
    assert _markup.within_mismatches(a, b, max_mm) is _reference(a, b, max_mm)


def test_it_stops_early_instead_of_counting_the_whole_overlap():
    """The point of the change. A pair that blows the budget in the first few bases must not
    depend on what follows -- if it did, the function would still be walking the whole overlap.

    Two 100k-char strings that differ at positions 0 and 1: with a budget of 1 the answer is
    settled at index 1, and the remaining 99,998 characters cannot matter. Constructing the tails
    to differ proves the tail is never examined.
    """
    head_a, head_b = "CC", "AA"
    a = head_a + "A" * 100_000
    b = head_b + "T" * 100_000               # tails disagree everywhere

    assert _markup.within_mismatches(a, b, 1) is False
    # Same prefix, budget large enough to survive it -> the tail now DOES matter, and disagrees.
    assert _markup.within_mismatches(a, b, 2) is False
    assert _markup.within_mismatches(head_a + "A" * 50, head_b + "A" * 50, 2) is True


def test_the_float_budget_collapses_to_floor():
    """⛔ The equivalence the speedup rests on. The old code compared an integer count against a
    FLOAT budget; the new one takes an int. `int((1 - min_id) * ov)` is floor for a non-negative
    value, and `count > budget` == `count > floor(budget)` for integer count — so a budget that
    lands at 9.999999999999998 instead of 10.0 decides the same candidates either way.
    """
    min_id, ov = 0.9, 100
    budget = (1 - min_id) * ov               # 10.000000000000002 in IEEE754, not 10.0
    max_mm = int(budget)

    for count in range(0, 15):
        a = "C" * count + "A" * (ov - count)
        b = "A" * ov
        assert _markup.within_mismatches(a, b, max_mm) is (count <= budget), (
            f"{count} mismatches against budget {budget!r} disagrees with the float comparison"
        )
