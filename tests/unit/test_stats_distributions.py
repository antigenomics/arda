"""The QC distributions, and the contract that makes a QC table joinable across samples.

Two things are pinned here, and each exists because its absence is silent:

* **shape, not summary.** ``junction_aa_len`` / ``read_len`` / ``clone_size`` / ``isotype`` are the
  four scopes that can show a bimodal junction length, a read-length cliff, an over-amplified
  library or a failed isotype primer. min/max/mean cannot, and every one of those reads as an
  ordinary mean.
* **one fact, one address.** A sample delivered as several lane FASTQs goes through
  ``_merge_map_reports``, which renames ``wall_seconds`` -> ``wall_seconds_max`` and
  ``peak_rss_mb`` -> ``peak_rss_mb_max``. Left alone, a cross-sample join on ``run/map/wall_seconds``
  drops exactly the lane-split samples and reports nothing.
"""

from __future__ import annotations

import json

import polars as pl
import pytest

from arda.stats import _flatten_report, collect, write_stats_json

# One IGH pair and one TRB read. `mmseqs2_qlen` is written by mmseqs as a FLOAT string (`407.0`),
# which is exactly how the read-length distribution first came out empty.
_AIRR_COLS = ["sequence_id", "locus", "v_call", "j_call", "junction", "junction_aa",
              "productive", "mmseqs2_qlen", "rev_comp", "junction_completed_nt"]
_AIRR_ROWS = [
    ["r1", "IGH", "IGHV3-9*02", "IGHJ4*02", "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF",
     "T", "407.0", "F", "0"],
    ["r2", "IGH", "IGHV3-9*02", "IGHJ4*02", "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF",
     "T", "402.0", "T", "3"],
    ["r3", "TRB", "TRBV19*01", "TRBJ2-7*01", "TGTGCCAGCAGCTTAGACGGGACATTC", "CASSLDGTF",
     "T", "75.0", "F", "0"],
    # A truncated junction: it has a LENGTH, and counting it would put a second, fake mode on the
    # left of every junction-length histogram. Its READ length still counts -- it is a read.
    ["r4", "TRB", "TRBV19*01", "TRBJ2-7*01", "TGTGCCAGCAGC", "CASSL", "F", "78.0", "F", "0"],
]
_CLONE_COLS = ["junction", "junction_aa", "v_call", "j_call", "c_call", "locus",
               "duplicate_count", "consensus_count", "chimera_parents"]
_CLONE_ROWS = [
    ["TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF", "IGHV3-9*02", "IGHJ4*02", "IGHG", "IGH",
     "17", "9", ""],
    ["TGTGCCAGCAGCTTAGACGGGACAGGGTGG", "CASSLDGTGW", "IGHV3-9*02", "IGHJ4*02", "IGHM", "IGH",
     "3", "3", ""],
    ["TGTGCCAGCAGCTTAGACGGGACATTC", "CASSLDGTF", "TRBV19*01", "TRBJ2-7*01", "", "TRB",
     "1", "1", ""],
]


@pytest.fixture
def run_dir(tmp_path):
    # Never: `quote_style="never"`, as every arda writer does — the default turns a blank
    # `chimera_parents` into the two-character value `""` and reads every clonotype as chimeric.
    pl.DataFrame(_AIRR_ROWS, schema=_AIRR_COLS, orient="row").write_csv(
        tmp_path / "s.airr.tsv", separator="\t", quote_style="never")
    pl.DataFrame(_CLONE_ROWS, schema=_CLONE_COLS, orient="row").write_csv(
        tmp_path / "s.clones.tsv", separator="\t", quote_style="never")
    return tmp_path


def _rows(run_dir):
    return collect(airr=run_dir / "s.airr.tsv", clones=run_dir / "s.clones.tsv")


def _scope(rows, scope, metric=None):
    return {key: value for s, key, m, value in rows
            if s == scope and (metric is None or m == metric)}


def test_junction_length_counts_only_reads_that_span_both_anchors(run_dir):
    """r4's junction is truncated, so it has no length to contribute."""
    rows = _rows(run_dir)
    assert _scope(rows, "junction_aa_len", "reads") == {"IGH:10": "2", "TRB:9": "1"}
    assert _scope(rows, "junction_aa_len", "clonotypes") == {"IGH:10": "2", "TRB:9": "1"}


def test_read_length_buckets_survive_a_float_valued_qlen_column(run_dir):
    """mmseqs writes `407.0`; a direct Int64 cast of that string is null, not 407."""
    # ...and the truncated read IS here: it was aligned, it just has no junction.
    assert _scope(_rows(run_dir), "read_len", "reads") == {"IGH:400": "2", "TRB:70": "2"}


def test_clone_size_buckets_are_powers_of_two_keyed_by_lower_bound(run_dir):
    rows = _rows(run_dir)
    assert _scope(rows, "clone_size", "clonotypes") == {"IGH:16": "1", "IGH:2": "1",
                                                        "TRB:1": "1"}
    # ...and the same buckets weighted by reads, which is the axis that shows over-amplification.
    assert _scope(rows, "clone_size", "reads") == {"IGH:16": "17", "IGH:2": "3", "TRB:1": "1"}


def test_isotype_is_reported_per_locus_and_only_where_there_is_evidence(run_dir):
    rows = _rows(run_dir)
    assert _scope(rows, "isotype", "clonotypes") == {"IGH:IGHG": "1", "IGH:IGHM": "1"}
    assert _scope(rows, "isotype", "reads") == {"IGH:IGHG": "17", "IGH:IGHM": "3"}


def test_strand_balance_and_completed_junctions_reach_the_sample_scope(run_dir):
    rows = _rows(run_dir)
    sample = {m: v for s, k, m, v in rows if s == "sample"}
    assert sample["reads_rev_comp"] == "1"
    assert sample["reads_junction_completed"] == "1"


def test_no_metric_is_ever_emitted_blank(run_dir):
    """A blank value casts to 0 downstream, which is the opposite of "no input"."""
    assert [r for r in _rows(run_dir) if r[3] == ""] == []


# ── the join contract ─────────────────────────────────────────────────────────────────────────

_SEGMENT = {"implied": 10, "rescued": 2, "fast_fraction": 0.1736,
            "reasons": {"v_only": 3, "no_combination": 7}}
_SINGLE = {"map": {"total_reads": 100, "wall_seconds": 12.5, "peak_rss_mb": 300,
                   "rss_gain_mb": 20, "segment_search": _SEGMENT}}
_MERGED = {"map": {"total_reads": 100, "shards": 4, "wall_seconds_max": 12.5,
                   "wall_seconds_sum": 40.0, "peak_rss_mb_max": 300,
                   "segment_search": _SEGMENT}}


def _metrics(report):
    out: list[tuple] = []
    _flatten_report(report, out)
    return {m for _, _, m, _ in out}


def test_a_lane_split_sample_reports_wall_time_under_the_same_name_as_a_single_file_one():
    """Never: the one defect that would corrupt every cohort table rather than failing."""
    shared = {"total_reads", "wall_seconds", "peak_rss_mb"}
    assert shared <= _metrics(_SINGLE)
    assert shared <= _metrics(_MERGED)
    # The renamed forms do not survive alongside the canonical one...
    assert "wall_seconds_max" not in _metrics(_MERGED)
    assert "peak_rss_mb_max" not in _metrics(_MERGED)
    # ...but `wall_seconds_sum` is a different fact and keeps its own name.
    assert "wall_seconds_sum" in _metrics(_MERGED)


def test_an_alias_never_overwrites_a_name_the_report_already_carries():
    both = {"map": {"wall_seconds": 1.0, "wall_seconds_max": 9.0}}
    out: list[tuple] = []
    _flatten_report(both, out)
    assert ("run", "map", "wall_seconds", "1") in out
    assert ("run", "map", "wall_seconds_max", "9") in out


def test_the_reasons_a_read_missed_the_fast_path_reach_the_table():
    """`reasons` is a dict inside a dict, and one level of flattening dropped it silently."""
    assert {"segment_search.reasons.v_only",
            "segment_search.reasons.no_combination"} <= _metrics(_SINGLE)


# ── the JSON twin ─────────────────────────────────────────────────────────────────────────────

def test_the_json_twin_is_typed_and_carries_exactly_the_tsv_rows(run_dir, tmp_path):
    rows = _rows(run_dir)
    write_stats_json(rows, tmp_path / "s.stats.json", arda_version="9.9.9")
    doc = json.loads((tmp_path / "s.stats.json").read_text())

    assert doc["arda_version"] == "9.9.9"
    assert doc["rows"] == len(rows)
    flat = {(scope, key, metric)
            for scope, keys in doc["stats"].items()
            for key, metrics in keys.items() for metric in metrics}
    assert flat == {(s, k, m) for s, k, m, _ in rows}
    # Typed, so a consumer can tell a count from a version string.
    assert doc["stats"]["sample"][""]["reads"] == 4
    assert isinstance(doc["stats"]["chain"]["IGH"]["junction_aa_mean"], (int, float))
