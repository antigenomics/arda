"""``arda stats`` — the run QC table.

The contract this file pins is the SHAPE, because a QC table is consumed by joining across
samples and every one of these is a silent break for that consumer:

* four columns, ``scope`` / ``key`` / ``metric`` / ``value``, one value per cell -- never a
  ``134/62`` hybrid;
* integers stay exact integers, not ``1.0e3``;
* a metric with no input is OMITTED, never emitted as 0 (a sample with no ``junction_quality``
  column must not read as "quality 0");
* the same fact has one address -- ``per_locus`` folds into the metric name whether it came from
  a bare ``map --report`` JSON or from a merged ``.arda.json``.
"""

from __future__ import annotations

import json

import polars as pl
import pytest

from arda.stats import STATS_COLUMNS, collect, write_stats

# Two IGH reads and one TRB, one of the IGH reads carrying a stop codon and a truncated junction.
_AIRR_COLS = ["sequence_id", "locus", "v_call", "j_call", "junction", "junction_aa",
              "productive", "stop_codon", "vj_in_frame", "v_identity", "v_mutations",
              "j_mutations", "junction_quality", "v_mutation_quality"]
_AIRR_ROWS = [
    ["r1", "IGH", "IGHV3-9*02", "IGHJ4*02", "TGTGCCAGCAGCTTAGACGGGACAGGGTTC",
     "CASSLDGTGF", "T", "F", "T", "0.98", "G100A,G200C", "", "IIIIIIIIII5IIIIIIIIIIIIIIIIIII",
     "40,20"],
    ["r2", "IGH", "IGHV3-9*02", "IGHJ4*02", "TGTGCCAGCAGCTTAGAC",
     "CASSLD", "F", "T", "T", "0.90", "G100A", "", "IIIIIIIIIIIIIIIIII", "38"],
    ["r3", "TRB", "TRBV19*01", "TRBJ2-7*01", "TGTGCCAGCAGCTTAGACGGGACAGGGTTC",
     "CASSLDGTGF", "T", "F", "T", "1.0", "", "", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII", ""],
]
_CLONE_COLS = ["junction", "junction_aa", "v_call", "j_call", "c_call", "locus",
               "duplicate_count", "consensus_count", "chimera_parents"]
_CLONE_ROWS = [
    ["TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF", "IGHV3-9*02", "IGHJ4*02", "IGHG", "IGH",
     "7", "5", ""],
    ["TGTGCCAGCAGCTTAGACGGGACAGGGTTA", "CASSLDGTGL", "IGHV3-9*02", "IGHJ4*02", "IGHG", "IGH",
     "2", "2", "IGHV3-9*02_IGHJ4*02_CASSLDGTGF"],
    ["TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF", "TRBV19*01", "TRBJ2-7*01", "", "TRB",
     "3", "3", ""],
]


@pytest.fixture
def run_dir(tmp_path):
    # ⛔ `quote_style="never"`, as every arda writer does. The default quotes an empty cell as
    # `""`, and `stats` reads with `quote_char=None` (so does `correct`) -- so a defaulted fixture
    # turns every blank `chimera_parents` into a 2-character value and every read chimeric.
    pl.DataFrame(_AIRR_ROWS, schema=_AIRR_COLS, orient="row").write_csv(
        tmp_path / "s.airr.tsv", separator="\t", quote_style="never")
    pl.DataFrame(_CLONE_ROWS, schema=_CLONE_COLS, orient="row").write_csv(
        tmp_path / "s.clones.tsv", separator="\t", quote_style="never")
    (tmp_path / "s.arda.json").write_text(json.dumps({
        "arda_version": "9.9.9", "wall_seconds": 12.5,
        "map": {"total_reads": 1000, "mapped_reads": 3, "paired": True, "input_bytes": 4096,
                "read_length_mean": 100.0, "threads": 8, "peak_rss_mb": 301.5,
                "per_locus": {"IGH": 2, "TRB": 1}},
        "correct": {"clonotypes_out": 3, "reads_assigned": 12},
        "assemble": None,
    }))
    return tmp_path


def _index(rows) -> dict:
    return {(scope, key, metric): value for scope, key, metric, value in rows}


def test_the_table_is_four_flat_columns(run_dir, tmp_path):
    rows = collect(airr=run_dir / "s.airr.tsv", clones=run_dir / "s.clones.tsv",
                   report=run_dir / "s.arda.json")
    write_stats(rows, tmp_path / "out.tsv")
    df = pl.read_csv(tmp_path / "out.tsv", separator="\t", infer_schema_length=0)
    assert df.columns == STATS_COLUMNS
    assert all("/" not in v for v in df["value"] if v is not None)
    assert df.height == len(rows)


def test_counts_are_exact_integers(run_dir):
    got = _index(collect(airr=run_dir / "s.airr.tsv", clones=run_dir / "s.clones.tsv"))
    assert got[("sample", "", "reads")] == "3"
    assert got[("chain", "IGH", "reads")] == "2"
    assert got[("chain", "TRB", "reads")] == "1"
    assert got[("sample", "", "clonotypes")] == "3"
    assert got[("sample", "", "clonotype_reads")] == "12"      # 7 + 2 + 3, not 12.0


def test_functionality_flags_are_counted_separately_from_truncation(run_dir):
    """A stop codon, an out-of-frame junction and a junction that does not reach [FW]118 are three
    different defects; ``_COMPLETE`` folds them together and a QC table must not."""
    got = _index(collect(airr=run_dir / "s.airr.tsv"))
    assert got[("chain", "IGH", "reads_truncated_junction")] == "1"    # r2, no trailing F/W
    assert got[("chain", "IGH", "reads_stop_codon")] == "0"            # junction_aa carries no *
    assert got[("chain", "IGH", "reads_nonfunctional")] == "1"         # r2, productive == F
    assert got[("chain", "IGH", "reads_out_of_frame")] == "0"


def test_junction_length_ignores_truncated_junctions(run_dir):
    """A prefix of a junction is not a short junction. Min/max over incomplete rows would report a
    length no rearrangement has."""
    got = _index(collect(airr=run_dir / "s.airr.tsv"))
    assert got[("chain", "IGH", "junction_nt_min")] == "30"            # not r2's 18
    assert got[("chain", "IGH", "junction_nt_max")] == "30"


def test_chimeras_are_counted_from_the_flag_column(run_dir):
    got = _index(collect(clones=run_dir / "s.clones.tsv"))
    assert got[("sample", "", "clonotypes_chimeric")] == "1"
    assert got[("sample", "", "chimeric_reads")] == "2"
    assert got[("chain", "TRB", "clonotypes_chimeric")] == "0"


def test_shm_rate_is_one_minus_mean_v_identity(run_dir):
    got = _index(collect(airr=run_dir / "s.airr.tsv"))
    assert float(got[("chain", "TRB", "shm_rate")]) == pytest.approx(0.0)
    assert float(got[("chain", "IGH", "shm_rate")]) == pytest.approx(1 - (0.98 + 0.90) / 2)


def test_mutation_quality_is_read_as_integers_not_phred33(run_dir):
    """``v_mutation_quality`` is comma-joined INTEGERS while ``junction_quality`` is Phred+33
    characters. Reading one as the other gives plausible numbers off by 33."""
    got = _index(collect(airr=run_dir / "s.airr.tsv", allele_min_frac=0.4, allele_min_reads=2))
    # G100A is on both IGHV3-9*02 reads at Phred 40 and 38 -> a candidate; G200C is on one.
    assert got[("allele_candidate", "IGHV3-9*02:G100A", "reads")] == "2"
    assert got[("allele_candidate", "IGHV3-9*02:G100A", "frequency")] == "1"
    assert float(got[("allele_candidate", "IGHV3-9*02:G100A", "mean_quality")]) == pytest.approx(39)
    assert got[("sample", "", "shm_variants")] == "1"                  # G200C, seen once
    assert float(got[("sample", "", "shm_variant_mean_quality")]) == pytest.approx(20)


def test_junction_quality_is_read_as_phred33(run_dir):
    got = _index(collect(airr=run_dir / "s.airr.tsv"))
    # r1 carries one '5' (Phred 20) among 'I' (Phred 40) over 30 bases.
    assert float(got[("chain", "IGH", "junction_quality_min_mean")]) == pytest.approx(30.0)


def test_a_missing_column_omits_its_metrics_rather_than_reporting_zero(tmp_path):
    """A sample run without ``--junction-quality`` must not read as "mean quality 0"."""
    pl.DataFrame([["r1", "TRB", "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF"]],
                 schema=["sequence_id", "locus", "junction", "junction_aa"],
                 orient="row").write_csv(tmp_path / "bare.airr.tsv", separator="\t",
                                         quote_style="never")
    metrics = {m for _, _, m, _ in collect(airr=tmp_path / "bare.airr.tsv")}
    assert "reads" in metrics
    assert "junction_quality_mean" not in metrics
    assert "shm_rate" not in metrics


def test_the_run_report_is_flattened_under_one_address(run_dir):
    got = _index(collect(report=run_dir / "s.arda.json"))
    assert got[("run", "map", "total_reads")] == "1000"
    assert got[("run", "map", "per_locus.IGH")] == "2"
    assert got[("run", "map", "peak_rss_mb")] == "301.5"
    assert got[("run", "correct", "reads_assigned")] == "12"
    assert got[("run", "", "wall_seconds")] == "12.5"


def test_a_bare_stage_report_addresses_per_locus_the_same_way(tmp_path):
    """A ``map --report`` JSON has no stage wrapper. ``per_locus`` must still be a metric name, or
    the same number lives in the key column for one file and the metric column for the other."""
    (tmp_path / "m.json").write_text(json.dumps(
        {"total_reads": 10, "per_locus": {"IGH": 4}}))
    got = _index(collect(report=tmp_path / "m.json"))
    assert got[("run", "", "per_locus.IGH")] == "4"
    assert got[("run", "", "total_reads")] == "10"


def test_gene_coverage_is_measured_against_the_shipped_germline_set(run_dir):
    got = _index(collect(airr=run_dir / "s.airr.tsv", clones=run_dir / "s.clones.tsv"))
    total = int(got[("sample", "", "v_genes_reference")])
    assert total > 100                                   # human carries hundreds of V genes
    # Two V genes were seen on reads; only IGHV3-9 on more than one.
    assert got[("v_gene", "IGHV3-9", "reads")] == "2"
    assert float(got[("sample", "", "v_gene_coverage_reads")]) == pytest.approx(2 / total)
    assert float(got[("sample", "", "v_gene_coverage_reads_multi")]) == pytest.approx(1 / total)


def test_library_shape_comes_from_the_inputs_not_the_airr(run_dir, tmp_path):
    """The AIRR holds the reads that MAPPED, so its row count and mate suffixes describe the
    receptor subset. Size and pairedness have to be recorded from the FASTQs."""
    r1, r2 = tmp_path / "a_1.fq", tmp_path / "a_2.fq"
    r1.write_text("@x\nACGT\n+\nIIII\n")
    r2.write_text("@x\nACGT\n+\nIIII\n")
    got = _index(collect(airr=run_dir / "s.airr.tsv", r1=r1, r2=r2))
    assert got[("sample", "", "paired")] == "1"
    assert got[("sample", "", "input_files")] == "2"
    assert got[("sample", "", "input_bytes")] == str(r1.stat().st_size + r2.stat().st_size)


def test_it_runs_on_any_one_input_alone(run_dir):
    for kw in ({"airr": run_dir / "s.airr.tsv"}, {"clones": run_dir / "s.clones.tsv"},
               {"report": run_dir / "s.arda.json"}):
        assert collect(**kw)
