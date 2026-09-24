"""``arda qc batch`` — many samples' QC tables as one cohort view.

What is pinned here is the arithmetic and the refusals, because both are silent when wrong:

* a metric one sample lacks must stay **null**, never 0 — the whole point of `stats` omitting a
  row rather than blanking it;
* a group too small for a scale estimate gets **no z**, not an invented one;
* a non-numeric value (a version, a file list) has no median and must not acquire one.
"""

from __future__ import annotations

import json

import polars as pl
import pytest

from arda import qc


def _stats(path, rows):
    pl.DataFrame(rows, schema=["scope", "key", "metric", "value"], orient="row").write_csv(
        path, separator="\t", quote_style="never")


@pytest.fixture
def batch_dir(tmp_path):
    """Six samples. S5 is a tenfold outlier on `reads`; S3 alone carries `pairing_rate`."""
    reads = {"S0": 100, "S1": 104, "S2": 96, "S3": 102, "S4": 98, "S5": 1000}
    for sid, n in reads.items():
        rows = [("sample", "", "reads", str(n)),
                ("sample", "", "clonotypes", "10"),
                ("run", "", "arda_version", "2.23.0"),
                ("run", "map", "mapped_fraction", "0.31"),
                ("chain", "IGH", "reads", str(n // 2)),
                ("junction_aa_len", "IGH:15", "reads", "3")]
        if sid == "S3":
            rows.append(("sample", "", "pairing_rate", "0.9"))
        _stats(tmp_path / f"{sid}.stats.tsv", rows)
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfastq_1\tproject\tbatch\n"
                     + "".join(f"{s}\t{s}_1.fq\tTRIAL9\tRUN1\n" for s in reads))
    for s in reads:
        (tmp_path / f"{s}_1.fq").write_text("")
    return tmp_path, sheet


def test_every_stats_table_in_the_directory_is_found_and_labelled(batch_dir):
    d, sheet = batch_dir
    df = qc.collect_batch(qc.find_stats(d), sheet=sheet)
    assert sorted(df["sample"].unique()) == ["S0", "S1", "S2", "S3", "S4", "S5"]
    assert set(df["project"]) == {"TRIAL9"} and set(df["batch"]) == {"RUN1"}


def test_a_sample_the_sheet_never_mentions_keeps_empty_labels_rather_than_vanishing(batch_dir):
    d, sheet = batch_dir
    _stats(d / "ORPHAN.stats.tsv", [("sample", "", "reads", "50")])
    df = qc.collect_batch(qc.find_stats(d), sheet=sheet)
    assert "ORPHAN" in df["sample"].to_list()
    assert df.filter(pl.col("sample") == "ORPHAN")["project"].to_list() == [""]


def test_the_wide_table_is_one_row_per_sample_and_a_missing_metric_stays_null(batch_dir):
    d, sheet = batch_dir
    wide = qc.pivot(qc.collect_batch(qc.find_stats(d), sheet=sheet))
    assert wide.height == 6
    assert wide.filter(pl.col("sample") == "S5")["reads"].to_list() == ["1000"]
    # Only S3 has it; the other five are null, and NOT "0" — a sample without single-cell
    # pairing has no pairing rate, which is a different statement from a rate of zero.
    assert wide.filter(pl.col("sample") == "S0")["pairing_rate"].to_list() == [None]


def test_the_run_scope_keeps_its_stage_in_the_column_name(batch_dir):
    d, sheet = batch_dir
    wide = qc.pivot(qc.collect_batch(qc.find_stats(d), sheet=sheet), scope="run")
    assert "map.mapped_fraction" in wide.columns
    assert "arda_version" in wide.columns


def test_the_robust_z_finds_the_tenfold_sample_and_flags_only_it(batch_dir):
    d, sheet = batch_dir
    scored = qc.outliers(qc.collect_batch(qc.find_stats(d), sheet=sheet))
    reads = scored.filter((pl.col("scope") == "sample") & (pl.col("metric") == "reads"))
    flagged = reads.filter(pl.col("outlier") == 1)["sample"].to_list()
    assert flagged == ["S5"]
    # median of 96,98,100,102,104,1000 is 101; MAD is 3 -> z = .6745*(1000-101)/3
    row = reads.filter(pl.col("sample") == "S5").row(0, named=True)
    assert row["median"] == pytest.approx(101.0)
    assert row["mad"] == pytest.approx(3.0)
    assert row["z"] == pytest.approx(0.6745 * (1000 - 101) / 3)


def test_a_group_too_small_for_a_scale_gets_a_median_but_no_z(tmp_path):
    for i, n in enumerate([10, 20, 30, 400]):
        _stats(tmp_path / f"T{i}.stats.tsv", [("sample", "", "reads", str(n))])
    scored = qc.outliers(qc.collect_batch(qc.find_stats(tmp_path)))
    assert scored["median"].to_list() == [25.0] * 4      # (20+30)/2
    assert scored["z"].to_list() == [None] * 4
    assert scored["outlier"].to_list() == [None] * 4


def test_a_metric_every_sample_agrees_on_gets_no_z_either(batch_dir):
    """MAD 0 would divide by zero, or flag whichever sample differs by a rounding error."""
    d, sheet = batch_dir
    scored = qc.outliers(qc.collect_batch(qc.find_stats(d), sheet=sheet))
    same = scored.filter(pl.col("metric") == "clonotypes")
    assert same["mad"].to_list() == [0.0] * 6
    assert same["z"].to_list() == [None] * 6


def test_a_version_string_never_acquires_a_median(batch_dir):
    d, sheet = batch_dir
    scored = qc.outliers(qc.collect_batch(qc.find_stats(d), sheet=sheet))
    assert "arda_version" not in scored["metric"].to_list()


def test_write_batch_lays_down_the_three_artifacts(batch_dir, tmp_path):
    d, sheet = batch_dir
    df = qc.collect_batch(qc.find_stats(d), sheet=sheet)
    long_p, wide_p, json_p = qc.write_batch(df, tmp_path / "cohort")
    assert long_p.name == "cohort.qc.tsv" and wide_p.name == "cohort.qc.wide.tsv"
    doc = json.loads(json_p.read_text())
    assert doc["samples"] == ["S0", "S1", "S2", "S3", "S4", "S5"]
    assert doc["labels"]["S0"] == {"project": "TRIAL9", "batch": "RUN1"}
    assert doc["stats"]["S5"]["sample"][""]["reads"] == 1000          # typed, not "1000"
    flagged = [o["sample"] for o in doc["outliers"]
               if o["metric"] == "reads" and o["scope"] == "sample"]
    assert flagged == ["S5"]


def test_an_empty_directory_is_an_empty_frame_not_a_crash(tmp_path):
    df = qc.collect_batch([])
    assert df.height == 0 and df.columns == qc.BATCH_COLUMNS
    assert qc.outliers(df).height == 0
    assert qc.pivot(df).height == 0
