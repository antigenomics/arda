"""The amino-acid D posterior: what it claims, and what it refuses to claim.

``posterior_d`` exists because a VDJdb record has no nucleotides and a D segment is trimmed
at both ends, so the D is usually invisible in the translated junction. Two sources survive:
the junction's nucleotide length (which pins ``insVD + |D surviving| + insDJ``, placing the D
even with zero sequence evidence) and the amino-acid match (which identifies it, but only
where enough D survives).

These tests pin the contract, not the accuracy -- accuracy is measured against OLGA ground
truth and against nucleotide D calls on real sequences, and reported in the module docstring.
The contract that matters is: **return nothing rather than guess.**
"""

from __future__ import annotations

import pytest

from arda.dpost import load_d_prior, posterior_d

HUMAN_TRB = ("CASSLAPGATNEKLFF", "TRBV5-1*01", "TRBJ1-4*01")


def test_priors_ship_only_where_a_model_exists():
    assert sorted(load_d_prior("human")) == ["IGH", "TRB", "TRD"]
    assert sorted(load_d_prior("mouse")) == ["TRB"]
    assert load_d_prior("rat") == {}, "no published model: ship nothing, not a human proxy"
    assert load_d_prior("rhesus_monkey") == {}


def test_organism_without_a_model_returns_none_rather_than_a_human_proxy():
    # A perfectly well-formed rhesus TRB junction. There is no rhesus model, so: nothing.
    assert posterior_d("CASSLGMSEPRWETQYF", "TRBV11-1*01", "TRBJ2-5*01", "rhesus_monkey") is None


def test_vj_locus_has_no_d_to_posterior_over():
    assert posterior_d("CAVRDSNYQLIW", "TRAV3*01", "TRAJ33*01", "human") is None


def test_posterior_is_a_distribution_over_genes():
    post = posterior_d(*HUMAN_TRB, "human")
    assert post is not None
    assert sum(post.by_gene.values()) == pytest.approx(1.0)
    assert post.posterior == pytest.approx(max(post.by_gene.values()))
    assert post.d_call in post.by_gene
    assert post.entropy >= 0.0


def test_the_junction_length_pins_the_middle_and_places_the_d_inside_it():
    junction, v_call, j_call = HUMAN_TRB
    post = posterior_d(junction, v_call, j_call, "human")
    assert post is not None
    assert post.n_middle_nt % 3 == 0, "the middle is a whole number of codons of the junction"
    assert 0 < post.n_middle_nt <= 3 * len(junction)
    lo, hi = post.d_start_ci90
    assert 0 <= lo <= post.d_start <= hi < 3 * len(junction)


def test_confident_is_the_posterior_crossing_nine_tenths():
    post = posterior_d(*HUMAN_TRB, "human")
    assert post is not None
    assert post.confident == (post.posterior >= 0.9)


def test_a_trbj1_junction_is_certain_of_trbd1_because_only_trbd1_can_reach_it():
    """Genomic order, not evidence: TRBD2 sits 3' of the whole TRBJ1 cluster."""
    post = posterior_d(*HUMAN_TRB, "human")
    assert post is not None
    assert post.d_call == "TRBD1"
    assert post.by_gene.get("TRBD2", 0.0) == 0.0
    assert post.posterior == pytest.approx(1.0) and post.entropy == pytest.approx(0.0)


def test_a_trbj2_junction_leaves_both_d_genes_live():
    """The negative control for the test above: with a J2 the posterior must be uncertain."""
    post = posterior_d("CASSLAPGATSYEQYF", "TRBV5-1*01", "TRBJ2-7*01", "human")
    assert post is not None
    assert set(post.by_gene) == {"TRBD1", "TRBD2"}
    assert all(p > 0.0 for p in post.by_gene.values())
    assert post.entropy > 0.0


def test_a_junction_the_model_cannot_explain_returns_none():
    """``insVD + |D| + insDJ`` is bounded. A 180 nt middle is outside the model's support."""
    assert posterior_d("CASS" + "G" * 60 + "YEQYF", "TRBV5-1*01", "TRBJ2-7*01", "human") is None


def test_overlapping_v_and_j_templates_return_none():
    """When the germlines explain more than the junction holds, there is no middle to place."""
    assert posterior_d("CASSLF", "TRBV5-1*01", "TRBJ1-4*01", "human") is None


def test_an_empty_middle_is_prior_only_and_says_so():
    """V and J templates abut: nothing of the D survives, so nothing identifies it.

    The posterior is then the prior alone. For a TRBJ1 that is still a certainty -- genomic
    order leaves only TRBD1 -- but for IGH, with 35 D genes and no such constraint, the
    posterior stays diffuse and ``confident`` is False. ``n_middle_nt`` and ``support_aa``
    are the honest signals that no sequence evidence was used.
    """
    trb = posterior_d("CASSF", "TRBV5-1*01", "TRBJ1-4*01", "human")
    assert trb is not None
    assert trb.n_middle_nt == 0 and trb.support_aa == 0
    assert trb.d_call == "TRBD1" and trb.posterior == pytest.approx(1.0)
    assert trb.entropy == 0.0 and not repr(trb.entropy).startswith("-")

    igh = posterior_d("CARW", "IGHV1-18*01", "IGHJ4*02", "human")
    assert igh is not None
    assert igh.n_middle_nt == 0 and igh.support_aa == 0
    assert not igh.confident, "35 D genes, no constraint: an empty middle proves nothing"
    assert igh.entropy > 1.0


# --- `--d-prior`: scoring against a fitted table without adopting it ---------------------------
# `arda scenarios` fits exactly the table `load_d_prior` reads, and `docs/scenarios.rst` calls its
# output a drop-in for the shipped file -- which it was by FORMAT only. `load_d_prior` was
# `@lru_cache`d on the organism alone and read one fixed path, so the only way to use a fitted
# table was to overwrite a file inside the installed database: a change that silently alters every
# later run on the machine. These pin the separation between *using* an estimate and *adopting* one.

def _toy_trb_prior(tmp_path):
    """A flat, hand-made TRB prior -- enough support to place a real junction's middle.

    Deliberately uninformative: both D genes at 0.5 and every length equally likely, so anything
    these tests observe comes from the table being READ, never from it being well fitted.
    """
    rows = ["locus\tkind\tkey\tvalue"]
    for i in range(16):
        rows.append(f"TRB\tinsVD\t{i}\t{1 / 16:.6f}")
        rows.append(f"TRB\tinsDJ\t{i}\t{1 / 16:.6f}")
    for allele in ("TRBD1*01", "TRBD2*01"):
        for n in range(1, 13):
            rows.append(f"TRB\tdlen\t{allele}:{n}\t{1 / 12:.6f}")
        rows.append(f"TRB\td_marginal\t{allele}\t0.5")
    rows.append("TRB\tbeta\tbeta\t1.25")
    path = tmp_path / "fitted.tsv"
    path.write_text("\n".join(rows) + "\n")
    return path


def test_a_supplied_table_gives_a_posterior_to_an_organism_that_ships_none(tmp_path):
    """11 of the 13 shipped (organism, D-locus) pairs have no prior, and rhesus TRB is one.

    That is the whole point of the parameter: `arda scenarios` can fit a table from a real
    cohort for a locus OLGA never modelled, and until now there was nowhere to put it.
    """
    args = ("CASSLGMSEPRWETQYF", "TRBV11-1*01", "TRBJ2-5*01", "rhesus_monkey")
    assert posterior_d(*args) is None, "no shipped rhesus model: still nothing by default"

    post = posterior_d(*args, prior_path=_toy_trb_prior(tmp_path))
    assert post is not None
    assert post.locus == "TRB"
    assert set(post.by_gene) == {"TRBD1", "TRBD2"}
    assert sum(post.by_gene.values()) == pytest.approx(1.0)


def test_the_shipped_database_is_not_touched_by_reading_one(tmp_path):
    """Using an estimate is not adopting it: the default answer must be unchanged afterwards."""
    before = posterior_d(*HUMAN_TRB, "human")
    posterior_d("CASSLGMSEPRWETQYF", "TRBV11-1*01", "TRBJ2-5*01", "rhesus_monkey",
                prior_path=_toy_trb_prior(tmp_path))
    after = posterior_d(*HUMAN_TRB, "human")
    assert after is not None and before is not None
    assert after.by_gene == before.by_gene
    assert sorted(load_d_prior("rhesus_monkey")) == []


def test_a_path_that_does_not_exist_raises_rather_than_scoring_on_nothing(tmp_path):
    """Never: the shipped table is ALLOWED to be missing; a path the user typed is not.

    `load_d_prior` returns `{}` for an organism with no model and `posterior_d` then returns
    `None` -- correct, and indistinguishable from "your table was not found", which is why a
    caller-supplied path has to raise instead.
    """
    with pytest.raises(FileNotFoundError, match="D prior table not found"):
        load_d_prior("human", tmp_path / "never_written.tsv")


def test_a_fitted_table_can_move_the_call_off_the_shipped_answer(tmp_path):
    """The table is the model. Drop a D gene from it and the posterior must follow."""
    rows = [ln for ln in _toy_trb_prior(tmp_path).read_text().splitlines()
            if "TRBD1*01" not in ln]
    only_d2 = tmp_path / "only_d2.tsv"
    only_d2.write_text("\n".join(rows) + "\n")

    post = posterior_d("CASSLAPGATSYEQYF", "TRBV5-1*01", "TRBJ2-7*01", "human",
                       prior_path=only_d2)
    assert post is not None
    assert set(post.by_gene) == {"TRBD2"}, "the shipped human table leaves both genes live"
