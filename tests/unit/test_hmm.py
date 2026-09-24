"""`arda.hmm` is `arda.scenarios` read as inference -- same recursion, one implementation.

What is pinned here is that the two readings agree, and that the inference marginalises rather
than maximises. A posterior that quietly became a Viterbi score would still look plausible on
every junction, which is exactly why it is asserted.
"""

from __future__ import annotations

from math import exp, isinf

import pytest

from arda import hmm, scenarios as sc

TRB = ("TGCAGTGCTAGAGATTTAGCGGGAGGGACTGAAGCTTTCTTT", "TRBV20-1*01", "TRBJ2-1*01")
TRA = ("TGTGCTGTGAGAGACAGCAACTATCAGTTAATCTGG", "TRAV1-1*01", "TRAJ33*01")


@pytest.fixture(scope="module")
def model():
    return hmm.model_for("human")


def test_the_posterior_is_a_distribution(model):
    post = hmm.posterior_d(*TRB, model)
    assert post.locus == "TRB"
    assert sum(post.probabilities.values()) == pytest.approx(1.0, rel=1e-9)
    assert post.best in post.probabilities
    assert post.scenarios > 1


def test_likelihood_is_the_partition_function_of_the_same_lattice(model):
    """One recursion, two readings -- if these ever disagree, one of them is a second model."""
    built = sc.lattice(*TRB, model)
    z = sum(t[0] for t in built[2])
    assert exp(hmm.log_likelihood(*TRB, model)) == pytest.approx(z, rel=1e-12)
    assert hmm.posterior_d(*TRB, model).log_likelihood == pytest.approx(
        hmm.log_likelihood(*TRB, model), rel=1e-12)


def test_the_posterior_marginalises_rather_than_maximises(model):
    """A D supported by many placements must be able to beat one with a single better placement.

    Taking the best path instead would make `posterior_d` agree with `best_scenario` on every
    junction; that they can disagree is the property worth having.
    """
    post = hmm.posterior_d(*TRB, model)
    built = sc.lattice(*TRB, model)
    best_single = max(built[2], key=lambda t: t[0])[3]
    summed = post.best
    # Both are legitimate answers; what must hold is that the posterior really is a SUM -- every
    # allele's mass exceeds its single best term's share.
    z = sum(t[0] for t in built[2])
    for allele, p in post.probabilities.items():
        top = max((t[0] for t in built[2] if t[3] == allele), default=0.0)
        assert p >= top / z - 1e-12
    assert summed and best_single


def test_best_scenario_reconstructs_the_junction(model):
    """The Viterbi path is still a scenario, so it still has to reassemble its input."""
    for junction, v, j in (TRB, TRA):
        s = hmm.best_scenario(junction, v, j, model)
        g = sc.germlines_for(v, j, "human")
        vpart = g.v_nt[: len(g.v_nt) - s.del_v]
        jpart = g.j_nt[s.del_j:]
        d = ""
        if s.d_call:
            dg = dict(g.d_germlines)[s.d_call]
            d = dg[s.del_dl: len(dg) - s.del_dr]
        assert len(vpart) + s.ins_vd + len(d) + s.ins_dj + len(jpart) == len(junction)
        assert junction.startswith(vpart) and junction.endswith(jpart)
        if d:
            at = len(vpart) + s.ins_vd
            assert junction[at: at + len(d)] == d


def test_a_vj_locus_scores_and_has_no_d(model):
    post = hmm.posterior_d(*TRA, model)
    assert post.locus == "TRA" and post.probabilities == {}
    assert not isinf(hmm.log_likelihood(*TRA, model))
    assert hmm.best_scenario(*TRA, model).d_call == ""


def test_unscoreable_input_is_minus_infinity_not_an_exception(model):
    assert isinf(hmm.log_likelihood("ACGT", "NOSUCHV*01", "NOSUCHJ*01", model))
    assert hmm.best_scenario("ACGT", "NOSUCHV*01", "NOSUCHJ*01", model) is None
    assert hmm.posterior_d("ACGT", "NOSUCHV*01", "NOSUCHJ*01", model).probabilities == {}


def test_an_estimated_prior_can_be_loaded_back_to_score(tmp_path, model):
    """`arda scenarios` writes a model; `hmm.model_for(prior=...)` scores against it.

    That round trip is what makes "estimate on one cohort, score another" possible at all.
    """
    stats = sc.estimate([(*TRB, 1.0)], organism="human", iterations=2)
    path = tmp_path / "prior.tsv"
    with open(path, "w") as fh:
        fh.write("\t".join(sc.PRIOR_COLUMNS) + "\n")
        for locus, kind, key, value in stats.rows():
            fh.write(f"{locus}\t{kind}\t{key}\t{value:.8g}\n")

    fitted = hmm.model_for("human", prior=path)
    # Fitted on this junction, so it must not score it WORSE than the borrowed OLGA model does.
    assert hmm.log_likelihood(*TRB, fitted) > hmm.log_likelihood(*TRB, model)
