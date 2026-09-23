"""External clustering validation: does our partition of the cells match a reference's?

Recall and concordance answer "did we find this chain". They do not answer "did we put the same
cells in the same clone", which is the question a repertoire analysis actually runs on, and which
a per-chain metric cannot see: splitting one clone in two and merging two clones into one both
leave every individual chain call correct.

The scores implemented here are the homogeneity-parsimony trade-off proposed in
`arXiv:2607.20799 <https://arxiv.org/abs/2607.20799>`_, whose reference implementation is
`qimmuno/clustereval <https://github.com/qimmuno/clustereval>`_. Three families, each a pair of
one score that punishes merging and one that punishes splitting:

======================  ==================  ==========================  =================
family                  punishes merging    punishes splitting          entropy
======================  ==================  ==========================  =================
information-theoretic   homogeneity         parsimony                   Shannon
set-matching            normalized purity   normalized inverse purity   min-entropy
pair-counting           pair specificity    pair sensitivity            collision
======================  ==================  ==========================  =================

Note: the formulas are implemented here rather than imported, so arda gains neither a dependency
nor an optional-import branch for sixty lines of counting. They are transcribed from the
reference implementation and :mod:`tests.unit.test_partition` checks them against its published
worked values; the trade-off is that a change upstream will not follow automatically, so the
version transcribed is recorded here: clustereval as of commit ``9f6daba``.

Never: **parsimony is not completeness.** Completeness normalises ``H(K|C)`` by ``H(K)``, the
entropy of the clustering being judged, so a method can raise its own score by making its
clustering more uniform. Parsimony normalises by ``log N - H(C)``, which depends only on the
reference -- so the denominator is the same for every method being compared, which is the whole
point of comparing them.
"""

from __future__ import annotations

import math
from collections import Counter

__all__ = [
    "homogeneity_score", "parsimony_score", "purity_score", "normalized_purity_score",
    "inverse_purity_score", "normalized_inverse_purity_score", "pair_specificity_score",
    "pair_sensitivity_score", "compare_partitions",
]


def _check(labels_true, labels_pred) -> tuple[list, list]:
    true, pred = list(labels_true), list(labels_pred)
    if len(true) != len(pred):
        raise ValueError(f"label arrays differ in length: {len(true)} against {len(pred)}")
    if not true:
        raise ValueError("empty label arrays")
    return true, pred


def _entropy(labels: list) -> float:
    counts = Counter(labels)
    if len(counts) == 1:
        return 0.0
    total = len(labels)
    return -sum((n / total) * math.log(n / total) for n in counts.values())


def _conditional_entropies(true: list, pred: list) -> tuple[float, float, float, float]:
    """``(H(C|K), H(K|C), H(C), H(K))`` in nats."""
    entropy_c, entropy_k = _entropy(true), _entropy(pred)
    joint = Counter(zip(true, pred))
    total = len(true)
    entropy_joint = -sum(n * math.log(n) for n in joint.values()) / total + math.log(total)
    return (max(0.0, entropy_joint - entropy_k), max(0.0, entropy_joint - entropy_c),
            entropy_c, entropy_k)


def homogeneity_score(labels_true, labels_pred) -> float:
    """1 when no cluster mixes two reference classes. Punishes MERGING."""
    true, pred = _check(labels_true, labels_pred)
    h_c_given_k, _, entropy_c, _ = _conditional_entropies(true, pred)
    return 1.0 if entropy_c == 0 else 1 - h_c_given_k / entropy_c


def parsimony_score(labels_true, labels_pred) -> float:
    """1 when no reference class is split across clusters. Punishes SPLITTING."""
    true, pred = _check(labels_true, labels_pred)
    _, h_k_given_c, entropy_c, _ = _conditional_entropies(true, pred)
    denominator = math.log(len(true)) - entropy_c
    return 1.0 if denominator == 0 else 1 - h_k_given_c / denominator


def _contingency(true: list, pred: list) -> dict[tuple, int]:
    return Counter(zip(true, pred))


def purity_score(labels_true, labels_pred) -> float:
    """Share of items in the majority reference class of their own cluster."""
    true, pred = _check(labels_true, labels_pred)
    best: dict = {}
    for (c, k), n in _contingency(true, pred).items():
        best[k] = max(best.get(k, 0), n)
    return sum(best.values()) / len(true)


def normalized_purity_score(labels_true, labels_pred) -> float:
    """Purity rescaled so that the one-cluster clustering scores 0 rather than a class share."""
    true, pred = _check(labels_true, labels_pred)
    purity = purity_score(true, pred)
    floor = max(Counter(true).values()) / len(true)
    return 1.0 if floor == 1.0 else (purity - floor) / (1.0 - floor)


def inverse_purity_score(labels_true, labels_pred) -> float:
    """Share of items in the majority cluster of their own reference class."""
    true, pred = _check(labels_true, labels_pred)
    best: dict = {}
    for (c, k), n in _contingency(true, pred).items():
        best[c] = max(best.get(c, 0), n)
    return sum(best.values()) / len(true)


def normalized_inverse_purity_score(labels_true, labels_pred) -> float:
    """Inverse purity rescaled so that the singleton clustering scores 0."""
    true, pred = _check(labels_true, labels_pred)
    inverse = inverse_purity_score(true, pred)
    floor = len(set(true)) / len(true)
    return 1.0 if floor == 1.0 else (inverse - floor) / (1.0 - floor)


def _pair_confusion(true: list, pred: list) -> tuple[int, int, int, int]:
    """``(tn, fp, fn, tp)`` over unordered pairs."""
    def pairs(counts) -> int:
        return sum(n * (n - 1) // 2 for n in counts)

    total = len(true) * (len(true) - 1) // 2
    same_both = pairs(_contingency(true, pred).values())
    same_true = pairs(Counter(true).values())
    same_pred = pairs(Counter(pred).values())
    tp = same_both
    fp = same_pred - same_both
    fn = same_true - same_both
    return total - tp - fp - fn, fp, fn, tp


def pair_specificity_score(labels_true, labels_pred) -> float:
    """Share of reference-different pairs we also kept apart. Punishes MERGING."""
    true, pred = _check(labels_true, labels_pred)
    tn, fp, _, _ = _pair_confusion(true, pred)
    return tn / (tn + fp) if (tn + fp) > 0 else 1.0


def pair_sensitivity_score(labels_true, labels_pred) -> float:
    """Share of reference-same pairs we also kept together. Punishes SPLITTING."""
    true, pred = _check(labels_true, labels_pred)
    _, _, fn, tp = _pair_confusion(true, pred)
    return tp / (tp + fn) if (tp + fn) > 0 else 1.0


def compare_partitions(labels_true, labels_pred) -> dict[str, float]:
    """All six scores plus the sizes, as one flat dict ready for a TSV row.

    Never: report BOTH scores of a family or neither. Either alone is trivially maximised --
    homogeneity by putting every item in its own cluster, parsimony by putting them all in one --
    which is exactly what the trade-off framing exists to prevent.

    Never: read ``classes_singleton`` before reading the scores. A reference partition that is
    almost all singletons carries almost no clustering to agree about -- an unexpanded T-cell
    repertoire is 476 clonotypes over 479 cells -- and every score is then decided by a handful
    of items. It is reported beside them for that reason, not as decoration.
    """
    true, pred = _check(labels_true, labels_pred)
    return {
        "items": len(true),
        "classes_reference": len(set(true)),
        "clusters_called": len(set(pred)),
        "classes_singleton": sum(1 for n in Counter(true).values() if n == 1),
        "homogeneity": homogeneity_score(true, pred),
        "parsimony": parsimony_score(true, pred),
        "normalized_purity": normalized_purity_score(true, pred),
        "normalized_inverse_purity": normalized_inverse_purity_score(true, pred),
        "pair_specificity": pair_specificity_score(true, pred),
        "pair_sensitivity": pair_sensitivity_score(true, pred),
    }
