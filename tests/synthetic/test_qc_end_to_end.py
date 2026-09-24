"""A cohort's QC table, built from real runs rather than from hand-written rows.

The batch layer reads only what a run wrote, so the one thing that can silently break it is a
run writing something different from what the batch layer expects. That is asserted here against
`tests/data/rnaseq_real` split into samples — the same fixture every other end-to-end claim in
this repo is measured on.
"""

import json

import polars as pl
import pytest

from arda import qc
from arda.cluster import split_pairs
from arda.qcreport import render
from arda.rnaseq import pipeline
from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

READS_1 = "tests/data/rnaseq_real/reads_1.fq.gz"
READS_2 = "tests/data/rnaseq_real/reads_2.fq.gz"


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    """Three samples from contiguous cuts of the fixture, run through `arda rnaseq`."""
    chunks = split_pairs(READS_1, tmp_path_factory.mktemp("chunks"), shards=3, r2=READS_2)
    out = tmp_path_factory.mktemp("cohort")
    for i, pair in enumerate(chunks):
        pipeline.run([pair], out, f"S{i}", threads=2)
    sheet = out / "sheet.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\tproject\tbatch\n" + "".join(
        f"S{i}\t{r1}\t{r2}\tTRIAL\tRUN1\n" for i, (r1, r2) in enumerate(chunks)))
    return out, sheet


def test_every_run_writes_both_qc_artifacts_without_being_asked(cohort):
    out, _ = cohort
    for i in range(3):
        assert (out / f"S{i}.stats.tsv").exists()
        assert (out / f"S{i}.stats.json").exists()


def _read_stats(path):
    """The QC TSV, with the empty `key` restored.

    Never: an unkeyed scope writes an EMPTY key cell, and polars reads an empty cell as null. So
    `sample`'s key round-trips as None unless it is filled, and a join on `key` then silently
    matches nothing for the scope that carries most of the metrics.
    """
    return (pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
            .with_columns(pl.col(["scope", "key", "metric", "value"]).fill_null("")))


def test_the_two_per_sample_artifacts_carry_exactly_the_same_facts(cohort):
    """Written from one `rows` list, so a divergence means someone added a second pass."""
    out, _ = cohort
    tsv = _read_stats(out / "S0.stats.tsv")
    doc = json.loads((out / "S0.stats.json").read_text())
    flat = {(scope, key, metric)
            for scope, keys in doc["stats"].items()
            for key, metrics in keys.items() for metric in metrics}
    assert flat == set(zip(tsv["scope"], tsv["key"], tsv["metric"]))
    assert doc["rows"] == tsv.height


def test_the_batch_table_reproduces_each_samples_own_numbers(cohort):
    """The roll-up must not transform anything — it concatenates and labels, nothing else."""
    out, sheet = cohort
    df = qc.collect_batch(qc.find_stats(out), sheet=sheet)
    assert sorted(df["sample"].unique()) == ["S0", "S1", "S2"]
    for sid in ("S0", "S1", "S2"):
        own = _read_stats(out / f"{sid}.stats.tsv")
        mine = df.filter(pl.col("sample") == sid).select(own.columns)
        assert mine.to_dicts() == own.to_dicts()


def test_the_cohort_carries_the_distributions_and_the_sheets_labels(cohort):
    out, sheet = cohort
    df = qc.collect_batch(qc.find_stats(out), sheet=sheet)
    assert {"junction_aa_len", "read_len", "clone_size"} <= set(df["scope"].unique())
    assert set(df["project"]) == {"TRIAL"} and set(df["batch"]) == {"RUN1"}


def test_a_three_sample_group_is_too_small_for_a_z_and_says_so(cohort):
    """Three samples cannot yield a scale estimate; a median without a z is the honest answer."""
    out, sheet = cohort
    scored = qc.outliers(qc.collect_batch(qc.find_stats(out), sheet=sheet))
    assert scored["median"].null_count() == 0
    assert scored["z"].null_count() == scored.height


def test_the_rendered_page_holds_every_sample_and_still_fetches_nothing(cohort, tmp_path):
    out, sheet = cohort
    df = qc.collect_batch(qc.find_stats(out), sheet=sheet)
    *_, json_path = qc.write_batch(df, tmp_path / "cohort")
    page = render(json_path, tmp_path / "cohort.qc.html").read_text()
    for sid in ("S0", "S1", "S2"):
        assert f'"{sid}"' in page
    assert "<script src=" not in page and "<link" not in page
