"""`--error-rate auto` — measuring the substitution rate from the library's own error cloud."""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from arda.cli import _error_rate, app
from arda.rnaseq.correct import estimate_error_rate

PARENT = "TGTGCGAGAGGGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
V, J = "IGHV3-23*01", "IGHJ4*02"


def cloud(n_children: int, child_count: int, parent_count: int):
    """A parent plus ``n_children`` one-substitution neighbours, each at ``child_count``."""
    kids = [PARENT[:i] + ("A" if PARENT[i] != "A" else "C") + PARENT[i + 1:]
            for i in range(n_children)]
    junctions = [PARENT] + kids
    counts = [parent_count] + [child_count] * n_children
    return junctions, counts, [V] * len(junctions), [J] * len(junctions)


def test_it_recovers_a_planted_rate():
    """The estimator must invert the model's own definition: `p_sub = error_rate * L`."""
    junctions, counts, v, j = cloud(10, 50, 100_000)
    rate, report = estimate_error_rate(junctions, counts, v, j)
    assert rate == pytest.approx((500 / 100_000) / len(PARENT))
    assert report["parents"] == 1
    assert report["parent_reads"] == 100_000
    assert report["child_reads"] == 500


def test_a_library_too_shallow_to_measure_returns_none():
    """Never: no number is better than an invented one.

    A library where nothing reaches the depth floor is a normal outcome -- three bulk RNA-seq arms
    in benchmark round 35 hit it -- and the caller has to see that rather than a default wearing a
    measurement's clothes.
    """
    junctions, counts, v, j = cloud(5, 2, 40)
    rate, report = estimate_error_rate(junctions, counts, v, j, min_parent=100)
    assert rate is None
    assert report["parents"] == 0


def test_a_cloud_of_zero_is_unmeasurable_not_zero():
    """An estimate of exactly 0 is outside the model's domain (`0 < p_err < 1`), not inside it."""
    rate, _report = estimate_error_rate([PARENT], [5_000], [V], [J])
    assert rate is None


def test_children_must_be_less_abundant_than_the_parent():
    """Equal-or-greater neighbours have no defined direction, so they are not cloud."""
    junctions, counts, v, j = cloud(2, 50, 100_000)
    counts[1] = 100_000                      # a co-equal neighbour, not an error child
    rate, report = estimate_error_rate(junctions, counts, v, j)
    assert report["child_reads"] == 50       # only the genuinely smaller one counts


def test_require_vj_separates_clonotypes_that_merely_look_alike():
    junctions, counts, v, j = cloud(3, 100, 50_000)
    v[1] = "IGHV1-69*01"                     # same junction distance, different rearrangement
    strict, _ = estimate_error_rate(junctions, counts, v, j, require_vj=True)
    loose, _ = estimate_error_rate(junctions, counts, v, j, require_vj=False)
    assert loose > strict


def test_the_estimate_scales_with_the_cloud():
    a, _ = estimate_error_rate(*cloud(10, 20, 100_000))
    b, _ = estimate_error_rate(*cloud(10, 40, 100_000))
    assert b == pytest.approx(2 * a)


def test_a_deep_clonotype_with_no_cloud_is_evidence_of_a_LOW_rate():
    """It counts in the denominator, and that is deliberate.

    Every clonotype at or above `min_parent` acts as a parent, including one whose cloud is empty.
    Skipping those would average only over clonotypes that HAVE errors and bias the rate upward --
    the quiet ones are exactly the evidence that the library is clean.
    """
    junctions, counts, v, j = cloud(10, 50, 100_000)
    noisy, _ = estimate_error_rate(junctions, counts, v, j)
    quiet_junction = "TGTGCGAGAGGGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAA"
    with_quiet, report = estimate_error_rate(
        junctions + [quiet_junction], counts + [100_000], v + [V], j + [J])
    assert report["parents"] == 2
    assert with_quiet < noisy


def test_cli_parses_auto_and_refuses_nonsense():
    assert _error_rate("auto") is None
    assert _error_rate("AUTO") is None
    assert _error_rate("1e-5") == pytest.approx(1e-5)
    assert _error_rate("0.001") == pytest.approx(0.001)
    with pytest.raises(Exception, match="auto"):
        _error_rate("probably-fine")


def test_cli_correct_accepts_auto(tmp_path):
    """End to end: the flag reaches the pipeline and the run reports what it resolved to."""
    airr = tmp_path / "in.airr.tsv"
    cols = ["sequence_id", "locus", "v_call", "j_call", "junction", "junction_aa", "productive"]
    rows = []
    kids = [PARENT[:i] + ("A" if PARENT[i] != "A" else "C") + PARENT[i + 1:] for i in range(6)]
    for n, seq in [(400, PARENT)] + [(3, k) for k in kids]:
        for r in range(n):
            rows.append([f"r{len(rows)}_{r}", "IGH", V, J, seq, "CARXW", "T"])
    airr.write_text("\t".join(cols) + "\n"
                    + "\n".join("\t".join(r) for r in rows) + "\n")
    out = tmp_path / "clones.tsv"
    result = CliRunner().invoke(app, ["correct", "-i", str(airr), "-o", str(out),
                                      "--error-rate", "auto"])
    assert result.exit_code == 0, result.output
    assert out.exists()


def test_cli_correct_refuses_a_bad_error_rate(tmp_path):
    airr = tmp_path / "in.airr.tsv"
    airr.write_text("sequence_id\tlocus\tv_call\tj_call\tjunction\tjunction_aa\tproductive\n")
    result = CliRunner().invoke(app, ["correct", "-i", str(airr), "-o",
                                      str(tmp_path / "out.tsv"), "--error-rate", "nonsense"])
    assert result.exit_code != 0
