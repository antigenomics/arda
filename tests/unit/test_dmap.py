"""Genomic order constrains which D can join which J -- assert the callers obey it.

The TRB locus runs TRBD1 - TRBJ1 cluster - TRBC1 - TRBD2 - TRBJ2 cluster - TRBC2. V(D)J
joining deletes the intervening DNA, so TRBD2 cannot reach any TRBJ1. Left unenforced,
TRBD2 (16 nt) simply outscores TRBD1 (12 nt) on noise: 17 % of real human TRB J1-cluster
records were assigned an impossible TRBD2.

The load-bearing case is ``test_planted_trbd2_is_refused_on_a_j1_junction``: a *verbatim,
untrimmed* TRBD2 germline is planted in the junction. The sequence evidence is perfect and
the caller must still refuse it, because the J says it cannot be there.
"""

from __future__ import annotations

import pytest

from arda.annotate.dmap import map_d_junction
from arda.annotate.reference import _load_d_germlines
from arda.annotate.transfer import _allowed_d
from arda.cdr3fix import load_anchors
from arda.dpost import _mask_forbidden, load_d_prior, posterior_d
from arda.paths import vdj_dir

J1, J2 = "TRBJ1-1*01", "TRBJ2-1*01"
V = {"human": "TRBV20-1*01", "mouse": "TRBV1*01"}
NP1, NP2 = "CTAAC", "GGATC"


def _gene(call: str) -> set[str]:
    return {a.strip().split("*")[0] for a in (call or "").split(",") if a.strip()}


def _junction(org: str, d_allele: str, j_call: str) -> str:
    """V germline + np1 + verbatim D + np2 + J germline, in junction space."""
    anc = load_anchors(org)
    d = dict(_load_d_germlines(vdj_dir(org))["TRB"])[d_allele]
    return (anc[("V", V[org])].germline_nt + NP1 + d + NP2
            + anc[("J", j_call)].germline_nt).upper()


@pytest.mark.parametrize("org", ["human", "mouse"])
def test_planted_trbd2_is_refused_on_a_j1_junction(org):
    """A perfect TRBD2 match must lose to genomic order, not to a score threshold."""
    call = map_d_junction(_junction(org, "TRBD2*01", J1), V[org], J1, org)
    assert "TRBD2" not in _gene(call.d_call), call.d_call
    assert "TRBD2" not in _gene(call.d2_call), call.d2_call


@pytest.mark.parametrize("org", ["human", "mouse"])
def test_the_same_planted_trbd2_is_called_on_a_j2_junction(org):
    """Negative control: the filter refuses D2 for J1, not D2 everywhere."""
    call = map_d_junction(_junction(org, "TRBD2*01", J2), V[org], J2, org)
    assert _gene(call.d_call) == {"TRBD2"}, call.d_call


@pytest.mark.parametrize("org", ["human", "mouse"])
def test_trbd1_still_maps_on_a_j1_junction(org):
    """The filter must not cost the D that genomic order does allow."""
    call = map_d_junction(_junction(org, "TRBD1*01", J1), V[org], J1, org)
    assert _gene(call.d_call) == {"TRBD1"}, call.d_call


def test_allowed_d_leaves_other_loci_and_ambiguous_j_alone():
    igh = _load_d_germlines(vdj_dir("human"))["IGH"]
    assert _allowed_d(igh, "IGHJ1*01") is igh          # every IGHD is 5' of every IGHJ
    trb = _load_d_germlines(vdj_dir("human"))["TRB"]
    assert _allowed_d(trb, "TRBJ2-1*01") is trb
    assert _allowed_d(trb, "TRBJ1-1*01,TRBJ2-1*01") is trb   # spans clusters: forbids nothing
    assert _allowed_d(trb, "") is trb                        # unknown J: forbids nothing
    assert not any(a.startswith("TRBD2") for a, _ in _allowed_d(trb, "TRBJ1-1*01"))


@pytest.mark.parametrize("org", ["human", "mouse"])
def test_shipped_prior_zeroes_trbd2_on_every_j1(org):
    prior = load_d_prior(org)["TRB"]
    for j_allele, row in prior.d_given_j.items():
        mass = sum(p for d, p in row.items() if d.startswith("TRBD2"))
        if j_allele.startswith("TRBJ1-"):
            assert mass == 0.0, f"{org} {j_allele}: P(TRBD2|J) = {mass}"
        else:
            assert mass > 0.0, f"{org} {j_allele}: TRBD2 wrongly excluded"
        assert abs(sum(row.values()) - 1.0) < 1e-6, f"{org} {j_allele} unnormalised"


def test_marginal_backoff_also_forbids_trbd2():
    """TRBJ1-6*01 has no row in the shipped human model, so it exercises the backoff path."""
    prior = load_d_prior("human")["TRB"]
    assert "TRBJ1-6*01" not in prior.d_given_j          # guards the premise of this test
    assert any(d.startswith("TRBD2") for d in prior.d_marginal)
    backed_off = _mask_forbidden(prior.d_marginal, "TRBJ1-6*01")
    assert not any(d.startswith("TRBD2") for d in backed_off)
    assert abs(sum(backed_off.values()) - 1.0) < 1e-6
    assert _mask_forbidden(prior.d_marginal, J2) is prior.d_marginal


@pytest.mark.parametrize("j_call", [J1, "TRBJ1-6*01"])
def test_posterior_never_calls_trbd2_on_a_j1(j_call):
    """A junction whose middle *is* TRBD2 still posteriors onto TRBD1, with no entropy left."""
    post = posterior_d("CSARDGTGGYSGANVLTF", "TRBV20-1*01", j_call, "human")
    assert post is not None
    assert post.d_call == "TRBD1" and post.by_gene.get("TRBD2", 0.0) == 0.0
    assert post.posterior == pytest.approx(1.0) and post.entropy == pytest.approx(0.0)
