"""The sharded RNA-seq path must be byte-identical to a single-node run.

The real proof needs mmseqs and a reference. These tests instead pin the *pieces* the
guarantee rests on, so CI catches a regression with no DB installed -- this repo has a
documented history of reference regressions slipping through precisely because the only tests
covering them were gated on `requires_human_db` and silently skipped.

The pieces:
  1. Stage 2-3 are byte-stable when their input is split contiguously and re-merged in order
     (with a shuffled-order negative control, so the test cannot pass vacuously).
  2. The generated SLURM script runs Stage 2-3 exactly once, in the reduce step, and never
     inside the array body.
  3. Per-shard Stage-1 reports merge into sums, and wall/RSS are never collapsed into a single
     misleading number.
"""

from __future__ import annotations

import json
import random

import polars as pl
import pytest

from arda.cluster import merge, render_rnaseq_submit_script
from arda.rnaseq import pipeline


def _airr_rows(n_clones: int = 8, per_clone: int = 4) -> list[dict]:
    """A Stage-1-shaped AIRR: complete, in-frame, canonical junctions."""
    rows = []
    for c in range(n_clones):
        junc = "TGTGCCAGCAGCTTAGACGGGACAGG" + ("TTC", "GTC", "CTC")[c % 3] + "T" * (c % 2)
        junc = junc[: (len(junc) // 3) * 3]
        for k in range(per_clone):
            rows.append({
                "sequence_id": f"read{c:02d}_{k}", "junction": junc,
                "junction_aa": "CASSLDGTF", "v_call": f"TRBV20-1*0{1 + c % 2}",
                "j_call": "TRBJ2-1*01", "locus": "TRB",
            })
    return rows


def _write(path, rows):
    pl.DataFrame(rows).write_csv(path, separator="\t")
    return path


def _partition_contiguously(rows, k):
    size = -(-len(rows) // k)
    return [rows[i:i + size] for i in range(0, len(rows), size)]


def test_stage23_identical_whether_input_arrived_whole_or_in_contiguous_shards(tmp_path):
    """The core equivalence property, without mmseqs."""
    rows = _airr_rows()
    whole = _write(tmp_path / "whole.airr.tsv", rows)
    pipeline.finish(whole, tmp_path / "single", "S", assemble=False)

    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    parts = _partition_contiguously(rows, 4)
    files = [_write(shard_dir / f"shard_{i:05d}.airr.tsv", p) for i, p in enumerate(parts)]
    merged = merge(sorted(files), tmp_path / "merged.airr.tsv")
    assert merged.read_bytes() == whole.read_bytes(), "ordered merge must rebuild the input exactly"

    pipeline.finish(merged, tmp_path / "sharded", "S", assemble=False)
    assert ((tmp_path / "sharded" / "S.clones.tsv").read_bytes()
            == (tmp_path / "single" / "S.clones.tsv").read_bytes())


def test_merging_shards_out_of_order_does_change_the_merged_airr(tmp_path):
    """Negative control: proves the previous test is not passing vacuously.

    If order never mattered, ordered merge would not be load-bearing and this suite would be
    asserting nothing. Shard names are zero-padded so `sorted()` is numeric; reverse that and
    the bytes must move.
    """
    rows = _airr_rows()
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()
    files = [_write(shard_dir / f"shard_{i:05d}.airr.tsv", p)
             for i, p in enumerate(_partition_contiguously(rows, 4))]
    in_order = merge(sorted(files), tmp_path / "a.tsv").read_bytes()
    reversed_ = merge(sorted(files, reverse=True), tmp_path / "b.tsv").read_bytes()
    assert in_order != reversed_


def test_reduce_reads_shards_in_numeric_not_lexicographic_order(tmp_path):
    """`shard_10` must not sort before `shard_2`; with >=10 shards this is the real trap."""
    rows = _airr_rows(n_clones=12, per_clone=1)
    shard_dir = tmp_path / "map"
    shard_dir.mkdir()
    for i, part in enumerate(_partition_contiguously(rows, 12)):
        _write(shard_dir / f"shard_{i:05d}.airr.tsv", part)

    pipeline.reduce(shard_dir, tmp_path / "out", "S", assemble=False)
    merged = pl.read_csv(tmp_path / "out" / "S.airr.tsv", separator="\t", infer_schema_length=0)
    assert merged["sequence_id"].to_list() == [r["sequence_id"] for r in rows]


def test_reduce_refuses_an_empty_shard_dir(tmp_path):
    with pytest.raises(FileNotFoundError, match="shard"):
        pipeline.reduce(tmp_path, tmp_path / "out", "S")


def test_finish_report_records_provenance(tmp_path):
    """A cross-mode divergence must be diagnosable: which arda, which mmseqs, which reference."""
    rows = _airr_rows()
    report = pipeline.finish(_write(tmp_path / "in.airr.tsv", rows), tmp_path, "S", assemble=False)
    for key in ("arda_version", "mmseqs_version", "reference", "wall_seconds", "correct"):
        assert key in report, key
    on_disk = json.loads((tmp_path / "S.arda.json").read_text())
    assert on_disk["arda_version"] == report["arda_version"]


def test_merge_map_reports_sums_counts_and_never_fakes_a_single_wall_time(tmp_path):
    shards = []
    for i, (total, mapped, wall, rss) in enumerate([(100, 10, 5.0, 300.0), (100, 20, 9.0, 280.0)]):
        p = tmp_path / f"shard_{i:05d}.map.json"
        p.write_text(json.dumps({
            "organism": "human", "total_reads": total, "mapped_reads": mapped,
            "per_locus": {"IGH": mapped}, "constant_only_fragments": 1, "isotype_from_mate": 2,
            "min_score": 75.0, "threads": 8, "wall_seconds": wall, "peak_rss_mb": rss}))
        shards.append(p)

    m = pipeline._merge_map_reports(shards)
    assert m["shards"] == 2
    assert m["total_reads"] == 200 and m["mapped_reads"] == 30
    assert m["per_locus"] == {"IGH": 30}
    assert m["constant_only_fragments"] == 2 and m["isotype_from_mate"] == 4
    assert m["wall_seconds_max"] == 9.0 and m["wall_seconds_sum"] == 14.0
    assert m["peak_rss_mb_max"] == 300.0
    # Summing 40 array tasks' wall time and calling it "wall_seconds" would be a lie.
    assert "wall_seconds" not in m and "peak_rss_mb" not in m


def test_submit_script_runs_stage23_once_and_never_in_the_array(tmp_path):
    s = render_rnaseq_submit_script("/d/r1.fq.gz", "SAMP", tmp_path, shards=8,
                                    r2="/d/r2.fq.gz", out_dir="/o", partition="medium")
    assert s.count("arda rnaseq reduce") == 1
    # The whole point: these must not appear as their own array steps.
    assert "arda rnaseq correct" not in s
    assert "arda rnaseq assemble" not in s
    assert "arda rnaseq map" in s and "--array=0-7" in s
    assert "--dependency=afterok:$SPLIT_JID" in s
    assert "--dependency=afterok:$ARRAY_JID" in s
    assert 'printf "%05d"' in s          # numeric shard names
    assert '[ -s "$f" ] || exit 0' in s  # an empty shard must not break the afterok chain
    assert "_R2.fastq" in s


def test_submit_script_omits_r2_when_single_end(tmp_path):
    s = render_rnaseq_submit_script("/d/r1.fq.gz", "SAMP", tmp_path, shards=3)
    assert "_R2.fastq" not in s and "--r2" not in s


def test_submit_script_threads_flags_through(tmp_path):
    s = render_rnaseq_submit_script("/d/r1.fq.gz", "S", tmp_path, shards=2, organism="mouse",
                                    kmer=13, min_score=0.0, reconstruct=True,
                                    assemble=False, map_d=False)
    assert "--organism mouse" in s and "--kmer 13" in s and "--min-score 0.0" in s
    assert "--reconstruct" in s and "--no-map-d" in s and "--no-assemble" in s
