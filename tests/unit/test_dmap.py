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
from arda.annotate.transfer import _allowed_d, _dd_orientation_ok
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


def test_allowed_d_excludes_nothing_for_J_POSITION_outside_the_trb_j1_case():
    """The J-position rule is TRB-J1-specific; everywhere else it must forbid nothing.

    ⚠ Compared by CONTENT, not identity: `_allowed_d` also drops `/OR` orphons unconditionally
    (see the next test), so it always returns a new list. An `is` check here would pass only by
    accident of that filter being absent.
    """
    def names(xs):
        return [a for a, _ in xs]

    igh = _load_d_germlines(vdj_dir("human"))["IGH"]
    igh_real = [x for x in igh if "/OR" not in x[0]]
    assert names(_allowed_d(igh, "IGHJ1*01")) == names(igh_real)   # every IGHD is 5' of every IGHJ
    trb = _load_d_germlines(vdj_dir("human"))["TRB"]
    assert names(_allowed_d(trb, "TRBJ2-1*01")) == names(trb)
    assert names(_allowed_d(trb, "TRBJ1-1*01,TRBJ2-1*01")) == names(trb)  # spans clusters
    assert names(_allowed_d(trb, "")) == names(trb)                       # unknown J
    assert not any(a.startswith("TRBD2") for a, _ in _allowed_d(trb, "TRBJ1-1*01"))


def test_orphan_d_genes_are_never_candidates():
    """`/OR` D genes sit OUTSIDE their locus and cannot rearrange at all.

    IMGT ships 10 of human IGH's 48 D alleles as `IGHD.../OR15-...`, which are on **chromosome
    15**. They are not a usage preference to down-weight -- they are not producible. Measured on a
    real bulk IGH library, **11 of 11 tandem D-D calls named `IGHD2/OR15-2a*01,IGHD2/OR15-2b*01`
    as their second D**, so the whole tandem-D signal there was this one vocabulary artifact.
    """
    igh = _load_d_germlines(vdj_dir("human"))["IGH"]
    assert any("/OR" in a for a, _ in igh), "fixture check: IMGT should ship orphon IGHD genes"
    for j in ("IGHJ1*01", "IGHJ4*02", ""):
        assert not any("/OR" in a for a, _ in _allowed_d(igh, j)), j
    # And the real genes all survive -- the filter must be surgical.
    assert len(_allowed_d(igh, "IGHJ4*02")) == sum(1 for a, _ in igh if "/OR" not in a)


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


def test_map_d_junction_partitions_the_interior_exactly():
    """``np1 + D + np2`` must reconstruct the V..J interior, byte for byte."""
    junction = _junction("human", "TRBD1*01", J2)
    call = map_d_junction(junction, V["human"], J2, "human")
    assert call.called and not call.is_dd
    d = junction[call.d_sequence_start - 1 : call.d_sequence_end]
    interior = junction[call.v_sequence_end : call.j_sequence_start - 1]
    assert call.np1 + d + call.np2 == interior
    assert call.np1 == NP1 and call.np2 == NP2
    assert float(call.d_support) < 0.2, "d_support is the E-value the call was gated on"
    # A single-allele call anchors to one germline, so the AIRR extras come with it.
    assert call.extra["d_germline_start"] == 1
    assert call.extra["d_cigar"].endswith("S")


def test_map_d_junction_is_empty_on_a_vj_locus():
    """TRA has no D gene. Not a miss -- there is nothing to find."""
    call = map_d_junction("TGTGCTGTGAGAGATAGCAACTATCAGTTAATCTGG",
                          "TRAV1-1*01", "TRAJ33*01", "human")
    assert call.locus == "TRA"
    assert not call.called and not call.is_dd
    assert call.d_call == "" and call.v_sequence_end == -1


def test_map_d_junction_refuses_an_unresolvable_allele():
    call = map_d_junction(_junction("human", "TRBD1*01", J2), "NOSUCHV*01", J2, "human")
    assert not call.called and call.v_sequence_end == -1


def test_map_d_junction_needs_an_interior_to_search():
    """V germline abutting J germline: no room for a D, so no call and no coordinates."""
    anc = load_anchors("human")
    junction = (anc[("V", V["human"])].germline_nt + anc[("J", J2)].germline_nt).upper()
    call = map_d_junction(junction, V["human"], J2, "human")
    assert not call.called and call.v_sequence_end == -1


@pytest.mark.parametrize("j_call", [J1, "TRBJ1-6*01"])
def test_posterior_never_calls_trbd2_on_a_j1(j_call):
    """A junction whose middle *is* TRBD2 still posteriors onto TRBD1, with no entropy left."""
    post = posterior_d("CSARDGTGGYSGANVLTF", "TRBV20-1*01", j_call, "human")
    assert post is not None
    assert post.d_call == "TRBD1" and post.by_gene.get("TRBD2", 0.0) == 0.0
    assert post.posterior == pytest.approx(1.0) and post.entropy == pytest.approx(0.0)


# ---------------------------------------------------------------------------------------------
# Tandem D-D: orientation, and the markup it has to hand a consumer.
#
# A D-D fusion is a rearrangement like any other -- the upstream D joins the downstream one and
# the DNA between them is deleted -- so the product carries them in GENOMIC order. The read's
# 5' D must therefore be the 5' gene. This is the same argument `_allowed_d` makes about D vs J,
# applied to the second D; before it, 10 of 15 tandem calls on a real TRB amplicon were the
# impossible TRBD2 -> TRBD1.
#
# ⛔ Every assertion here is about a CALL (which gene, is there a second one) or about the
# partition CLOSING (parts concatenate back to the junction), never about where a boundary
# inside the NDN falls. The plants are untrimmed germlines, so `np1`/`np2`/`np3` come back
# verbatim -- that is a property of the fixture, not a claim that arda can recover chew-back.
NP3 = "TTAGC"

_TANDEM = {
    # locus: (V, J, producible 5'->3' pair)
    "TRB": ("TRBV20-1*01", J2, ("TRBD1*01", "TRBD2*01")),
    "TRD": ("TRDV1*01", "TRDJ1*01", ("TRDD2*01", "TRDD3*01")),
    "IGH": ("IGHV3-23*01", "IGHJ4*02", ("IGHD1-1*01", "IGHD6-19*01")),
}


def _tandem_junction(locus: str, d5: str, d3: str) -> str:
    """V germline + np1 + D + np2 + D + np3 + J germline, in junction space."""
    v_call, j_call, _ = _TANDEM[locus]
    anc = load_anchors("human")
    d = dict(_load_d_germlines(vdj_dir("human"))[locus])
    return (anc[("V", v_call)].germline_nt + NP1 + d[d5] + NP2 + d[d3] + NP3
            + anc[("J", j_call)].germline_nt).upper()


@pytest.mark.parametrize("locus", ["TRB", "TRD", "IGH"])
def test_tandem_dd_in_genomic_order_is_called_and_partitions_exactly(locus):
    """The producible direction: both Ds named, and np1+D+np2+D+np3 == the interior."""
    v_call, j_call, (d5, d3) = _TANDEM[locus]
    junction = _tandem_junction(locus, d5, d3)
    call = map_d_junction(junction, v_call, j_call, "human")
    assert call.is_dd, f"{d5} -> {d3} is producible and must be called"
    assert _gene(call.d_call) == _gene(d5) and _gene(call.d2_call) == _gene(d3)

    d1 = junction[call.d_sequence_start - 1 : call.d_sequence_end]
    d2 = junction[call.d2_sequence_start - 1 : call.d2_sequence_end]
    interior = junction[call.v_sequence_end : call.j_sequence_start - 1]
    assert call.np1 + d1 + call.np2 + d2 + call.np3 == interior
    assert (call.np1, call.np2, call.np3) == (NP1, NP2, NP3)
    assert call.d_sequence_end < call.d2_sequence_start, "the two Ds must not overlap"


@pytest.mark.parametrize("locus", ["TRB", "TRD"])
def test_tandem_dd_against_genomic_order_is_refused(locus):
    """The same two germlines, swapped. Deletional joining cannot make this product.

    The sequence evidence is *perfect* -- both germlines are planted verbatim -- and the second
    call must still go, exactly as a planted TRBD2 must lose to a TRBJ1. What is left is the
    single, higher-scoring D; the refused one is absorbed into ``np2``.
    """
    v_call, j_call, (d5, d3) = _TANDEM[locus]
    junction = _tandem_junction(locus, d3, d5)          # swapped: 3' gene planted 5'
    call = map_d_junction(junction, v_call, j_call, "human")
    assert call.called, "refusing the pair must not cost the single D"
    assert not call.is_dd, f"{d3} -> {d5} contradicts genomic order: {call.d2_call}"
    assert call.d2_sequence_start == -1 and call.np3 == ""


@pytest.mark.parametrize("locus", ["TRB", "TRD", "IGH"])
def test_tandem_markup_closes_over_the_whole_junction(locus):
    """``DCall.markup`` is the D-D deliverable: labelled parts that rebuild the junction."""
    v_call, j_call, (d5, d3) = _TANDEM[locus]
    junction = _tandem_junction(locus, d5, d3)
    call = map_d_junction(junction, v_call, j_call, "human")
    parts = call.markup(junction)
    assert "".join(seq for _, seq in parts) == junction
    assert [name for name, _ in parts] == [
        "V", "np1", call.d_call, "np2", call.d2_call, "np3", "J"]


def test_markup_closes_for_a_single_d_and_for_no_d_at_all():
    """The same contract on the other two shapes, so a consumer needs no special-casing."""
    single = _junction("human", "TRBD1*01", J2)
    call = map_d_junction(single, V["human"], J2, "human")
    assert [n for n, _ in call.markup(single)] == ["V", "np1", call.d_call, "np2", "J"]
    assert "".join(s for _, s in call.markup(single)) == single

    # No D: the interior is one unlabelled N stretch, and the parts still close.
    nod = map_d_junction(single, V["human"], J2, "human", d_max_evalue=1e-300)
    assert not nod.called
    assert [n for n, _ in nod.markup(single)] == ["V", "N", "J"]
    assert "".join(s for _, s in nod.markup(single)) == single

    # Nothing to mark up at all -> empty, not a half-built partition.
    assert map_d_junction("", V["human"], J2, "human").markup("") == []


def test_igh_tandem_is_accepted_in_BOTH_directions_and_still_partitions():
    """⛔ IGH carries no orientation constraint, on purpose -- assert that, don't assume it.

    `IGHD<family>-<position>` encodes genomic position in HUMAN IMGT and a family-member index
    in MOUSE, and the two vocabularies collide on real gene names, while `_map_d` is handed
    sequences and not an organism. So IGH is left out of `_D_GENOMIC_ORDER` and both directions
    are called. The markup contract still has to hold on the one arda does not gate.
    """
    _, _, (d5, d3) = _TANDEM["IGH"]
    v_call, j_call = _TANDEM["IGH"][0], _TANDEM["IGH"][1]
    for a, b in ((d5, d3), (d3, d5)):
        junction = _tandem_junction("IGH", a, b)
        call = map_d_junction(junction, v_call, j_call, "human")
        assert call.is_dd, f"IGH {a} -> {b} must not be orientation-gated"
        assert _gene(call.d_call) == _gene(a) and _gene(call.d2_call) == _gene(b)
        assert "".join(seq for _, seq in call.markup(junction)) == junction


def test_dd_orientation_rule_runs_in_both_directions():
    """The rule itself, without a junction in the way."""
    assert _dd_orientation_ok(("TRBD1*01",), ("TRBD2*01",))
    assert not _dd_orientation_ok(("TRBD2*01",), ("TRBD1*01",))
    assert _dd_orientation_ok(("TRDD1*01",), ("TRDD3*01",))
    assert not _dd_orientation_ok(("TRDD3*01",), ("TRDD1*01",))
    # The same gene twice needs two germline copies: not producible either.
    assert not _dd_orientation_ok(("TRBD1*01",), ("TRBD1*02",))
    # An ambiguity list passes on ONE producible assignment -- the call never claimed which.
    assert _dd_orientation_ok(("TRBD1*01", "TRBD2*01"), ("TRBD2*01",))
    # ...but a list whose every member sits 3' of the other side is still refused.
    assert not _dd_orientation_ok(("TRBD2*01",), ("TRBD1*01", "TRBD2*01"))
    # IGH is deliberately unranked (mouse `IGHD<f>-<n>` is not genomic order), so no constraint.
    assert _dd_orientation_ok(("IGHD7-27*01",), ("IGHD1-1*01",))


def test_d_max_evalue_tightens_the_gate_and_the_default_is_unchanged():
    """The strict band is opt-in; passing None must reproduce the shipped 0.2 exactly."""
    junction = _junction("human", "TRBD1*01", J2)
    shipped = map_d_junction(junction, V["human"], J2, "human")
    assert shipped == map_d_junction(junction, V["human"], J2, "human", d_max_evalue=0.2)
    assert shipped.called and float(shipped.d_support) < 0.2

    # Strict enough to refuse a call whose own E-value it undercuts, loose enough to keep it.
    e = float(shipped.d_support)
    assert map_d_junction(junction, V["human"], J2, "human", d_max_evalue=e / 100).d_call == ""
    assert map_d_junction(junction, V["human"], J2, "human", d_max_evalue=e * 100).called
