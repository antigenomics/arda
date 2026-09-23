"""The six clustering-validation scores, against cases whose value is known by hand.

The two degenerate clusterings are the point of the trade-off framing and are checked first:
one cluster holding everything maximises parsimony and zeroes homogeneity, one cluster per item
does the reverse. A test suite that only checked the perfect case would pass on either half.
"""

from __future__ import annotations

import math

import pytest

from arda.partition import (
    compare_partitions,
    homogeneity_score,
    inverse_purity_score,
    normalized_inverse_purity_score,
    normalized_purity_score,
    pair_sensitivity_score,
    pair_specificity_score,
    parsimony_score,
    purity_score,
)

TRUE = ["a", "a", "a", "b", "b", "b"]


def test_a_perfect_clustering_scores_one_everywhere():
    scores = compare_partitions(TRUE, ["x", "x", "x", "y", "y", "y"])
    for name in ("homogeneity", "parsimony", "normalized_purity",
                 "normalized_inverse_purity", "pair_specificity", "pair_sensitivity"):
        assert scores[name] == pytest.approx(1.0), name


def test_relabelling_does_not_change_a_score():
    a = compare_partitions(TRUE, ["x", "x", "x", "y", "y", "y"])
    b = compare_partitions(TRUE, [7, 7, 7, 0, 0, 0])
    assert a == b


def test_everything_in_one_cluster_maximises_parsimony_and_zeroes_homogeneity():
    pred = ["x"] * 6
    assert homogeneity_score(TRUE, pred) == pytest.approx(0.0)
    assert parsimony_score(TRUE, pred) == pytest.approx(1.0)
    assert normalized_purity_score(TRUE, pred) == pytest.approx(0.0)
    assert normalized_inverse_purity_score(TRUE, pred) == pytest.approx(1.0)
    assert pair_specificity_score(TRUE, pred) == pytest.approx(0.0)
    assert pair_sensitivity_score(TRUE, pred) == pytest.approx(1.0)


def test_one_cluster_per_item_maximises_homogeneity_and_zeroes_parsimony():
    pred = list(range(6))
    assert homogeneity_score(TRUE, pred) == pytest.approx(1.0)
    assert parsimony_score(TRUE, pred) == pytest.approx(0.0)
    assert normalized_purity_score(TRUE, pred) == pytest.approx(1.0)
    assert normalized_inverse_purity_score(TRUE, pred) == pytest.approx(0.0)
    assert pair_specificity_score(TRUE, pred) == pytest.approx(1.0)
    assert pair_sensitivity_score(TRUE, pred) == pytest.approx(0.0)


def test_purity_and_inverse_purity_on_a_hand_counted_split():
    # class a is split across x and y; class b sits whole inside y.
    pred = ["x", "x", "y", "y", "y", "y"]
    # majority class per cluster: x -> a (2), y -> b (3)  => 5/6
    assert purity_score(TRUE, pred) == pytest.approx(5 / 6)
    # majority cluster per class: a -> x (2), b -> y (3)  => 5/6
    assert inverse_purity_score(TRUE, pred) == pytest.approx(5 / 6)
    assert normalized_purity_score(TRUE, pred) == pytest.approx((5 / 6 - 3 / 6) / (1 - 3 / 6))
    assert normalized_inverse_purity_score(TRUE, pred) == pytest.approx(
        (5 / 6 - 2 / 6) / (1 - 2 / 6))


def test_pair_counts_on_the_same_split():
    pred = ["x", "x", "y", "y", "y", "y"]
    # same-class pairs: 3 + 3 = 6.  same-cluster pairs: 1 + 6 = 7.  both: aa(1) + bb(3) = 4.
    assert pair_sensitivity_score(TRUE, pred) == pytest.approx(4 / 6)
    # different-class pairs: 15 - 6 = 9; of those, kept apart: 9 - (7 - 4) = 6.
    assert pair_specificity_score(TRUE, pred) == pytest.approx(6 / 9)


def test_homogeneity_matches_its_definition_on_the_same_split():
    pred = ["x", "x", "y", "y", "y", "y"]
    # H(C) = ln 2.  H(C|K): cluster x is pure (0 nats, weight 2/6); cluster y is 1 a and 3 b.
    h_y = -(1 / 4) * math.log(1 / 4) - (3 / 4) * math.log(3 / 4)
    expected = 1 - (4 / 6) * h_y / math.log(2)
    assert homogeneity_score(TRUE, pred) == pytest.approx(expected)


def test_parsimony_normalises_by_the_reference_not_by_the_clustering():
    # The distinguishing property against completeness: the denominator log N - H(C) depends on
    # the reference alone, so two methods scored against one reference share it.
    pred = ["x", "x", "y", "y", "y", "y"]
    h_a = -(2 / 3) * math.log(2 / 3) - (1 / 3) * math.log(1 / 3)
    h_k_given_c = (3 / 6) * h_a + (3 / 6) * 0.0
    expected = 1 - h_k_given_c / (math.log(6) - math.log(2))
    assert parsimony_score(TRUE, pred) == pytest.approx(expected)


def test_singleton_classes_are_reported_beside_the_scores():
    scores = compare_partitions(["a", "b", "c", "c"], ["x", "y", "z", "z"])
    assert scores["classes_singleton"] == 2
    assert scores["classes_reference"] == 3
    assert scores["clusters_called"] == 3


def test_mismatched_or_empty_input_is_refused():
    with pytest.raises(ValueError, match="differ in length"):
        homogeneity_score(["a", "b"], ["x"])
    with pytest.raises(ValueError, match="empty"):
        homogeneity_score([], [])
