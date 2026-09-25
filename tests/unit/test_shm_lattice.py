"""The SHM-aware lattice: `scenarios.lattice(..., shm=)` prices a mutated V tail.

Without a model the templated V length is bounded by an EXACT common prefix, so one substitution
in the V tail forces the rest of it to be re-read as N-region. These pin both readings and the
invariant that ties them together.
"""

from __future__ import annotations

import pytest

from arda.hmm import best_scenario, log_likelihood, model_for, posterior_d
from arda.scenarios import germlines_for, lattice
from arda.shmmodel import ShmModel

V_CALL, J_CALL = "IGHV3-30*18", "IGHJ4*02"


@pytest.fixture(scope="module")
def pair():
    g = germlines_for(V_CALL, J_CALL, "human")
    if g is None or len(g.v_nt) < 8:
        pytest.skip("the shipped human reference has no usable IGHV3-30*18 / IGHJ4*02 pair")
    return g


@pytest.fixture(scope="module")
def model():
    return model_for("human")


@pytest.fixture(scope="module")
def shm():
    """A model hot enough to price a mismatch at all; the numbers are the measured order."""
    return ShmModel(k=5, scale=0.03, context={}, region={})


def junctions(g):
    """``(clean, mutated)`` -- the same junction with one substitution in the templated V tail."""
    clean = g.v_nt + "GGGACTAC" + g.j_nt[-12:]
    i = 2
    sub = "A" if g.v_nt[i] != "A" else "C"
    return clean, clean[:i] + sub + clean[i + 1:]


def test_without_a_model_one_substitution_costs_the_whole_v_tail(pair, model):
    """The behaviour S2 exists to fix, pinned so the fix is visible as a change."""
    _clean, mutated = junctions(pair)
    sc = best_scenario(mutated, V_CALL, J_CALL, model, "human")
    assert sc.del_v >= len(pair.v_nt) - 3


def test_a_model_recovers_the_templated_tail(pair, model, shm):
    _clean, mutated = junctions(pair)
    exact = best_scenario(mutated, V_CALL, J_CALL, model, "human")
    scored = best_scenario(mutated, V_CALL, J_CALL, model, "human", shm)
    assert scored.del_v < exact.del_v
    assert scored.del_v == 0
    assert log_likelihood(mutated, V_CALL, J_CALL, model, "human", shm) > \
        log_likelihood(mutated, V_CALL, J_CALL, model, "human")


def test_a_clean_junction_keeps_its_reading(pair, model, shm):
    """The model must not rewrite a junction that never needed it."""
    clean, _mutated = junctions(pair)
    exact = best_scenario(clean, V_CALL, J_CALL, model, "human")
    scored = best_scenario(clean, V_CALL, J_CALL, model, "human", shm)
    assert (scored.del_v, scored.ins_vd, scored.d_call) == \
        (exact.del_v, exact.ins_vd, exact.d_call)


def test_a_zero_rate_model_is_exactly_the_exact_bound(pair, model):
    """Never: the exact-match bound IS the mu = 0 case of the emission, not a separate rule.

    With every rate 0 a mismatch costs `mu / 3 == 0`, so no templated length past the common
    prefix survives -- which is what `_common_prefix` computes. If these two ever disagree, the
    emission has stopped generalising the thing it replaced.
    """
    zero = ShmModel(k=5, scale=0.0, context={}, region={})
    for junction in junctions(pair):
        a = lattice(junction, V_CALL, J_CALL, model, "human")
        b = lattice(junction, V_CALL, J_CALL, model, "human", zero)
        assert [t[0] for t in a[2]] == pytest.approx([t[0] for t in b[2]])
        assert log_likelihood(junction, V_CALL, J_CALL, model, "human") == pytest.approx(
            log_likelihood(junction, V_CALL, J_CALL, model, "human", zero))


def test_an_unknown_allele_falls_back_rather_than_raising(pair, model, shm):
    """An allele with no full germline (or one that fails the suffix check) scores as today."""
    from arda import scenarios

    _clean, mutated = junctions(pair)
    original = scenarios._v_full_germlines
    scenarios._v_full_germlines = lambda organism: {}
    try:
        fallback = best_scenario(mutated, V_CALL, J_CALL, model, "human", shm)
    finally:
        scenarios._v_full_germlines = original
    exact = best_scenario(mutated, V_CALL, J_CALL, model, "human")
    assert (fallback.del_v, fallback.ins_vd) == (exact.del_v, exact.ins_vd)


def test_the_d_posterior_is_still_a_distribution(pair, model, shm):
    _clean, mutated = junctions(pair)
    post = posterior_d(mutated, V_CALL, J_CALL, model, "human", shm)
    assert post.probabilities
    assert sum(post.probabilities.values()) == pytest.approx(1.0)
    assert post.best in post.probabilities


def test_the_germline_map_keys_every_member_of_a_tie_group():
    """Never: a scaffold's `v_call` may be a comma-joined GROUP.

    `IGHV3-23*01` lives inside `IGHV3-23*01,IGHV3-23D*01`. Keying only the group string means a
    lookup by a single allele silently finds no germline -- which is how the first version of the
    SHM lattice scored exactly like the exact one and looked like a null result.
    """
    from arda.shmmodel import germline_map

    mapping = germline_map("human", "IGH")
    grouped = [k for k in mapping if "," in k]
    assert grouped, "the human IGH reference has tie-group v_calls"
    for key in grouped[:5]:
        for member in key.split(","):
            assert member in mapping
            assert mapping[member] == mapping[key]


def test_estimate_accepts_a_model_and_moves_the_fit(pair, shm):
    """The E-step's half of the same parameter: with a model, a mutated V tail stops being
    counted as trimming and insertion."""
    from arda.scenarios import estimate

    clean, mutated = junctions(pair)
    records = [(mutated, V_CALL, J_CALL, 1.0)] * 3 + [(clean, V_CALL, J_CALL, 1.0)]
    plain = estimate(records, organism="human", iterations=1)
    scored = estimate(records, organism="human", iterations=1, shm=shm)
    assert plain.records == scored.records == len(records)

    def mass(stats, kind):
        # `counts` is keyed (locus, kind, key) -> expected count.
        return {key: w for (_locus, k, key), w in stats.counts.items() if k == kind}

    assert mass(plain, "delV") != mass(scored, "delV")


def test_cli_scenarios_refuses_a_missing_shm_table(tmp_path):
    from typer.testing import CliRunner

    from arda.cli import app

    clones = tmp_path / "clones.tsv"
    clones.write_text("junction\tv_call\tj_call\nTGTGCGAAAGA\t%s\t%s\n" % (V_CALL, J_CALL))
    result = CliRunner().invoke(app, ["scenarios", "-i", str(clones), "-o",
                                      str(tmp_path / "out.tsv"), "--shm-model",
                                      str(tmp_path / "absent.tsv")])
    assert result.exit_code != 0
