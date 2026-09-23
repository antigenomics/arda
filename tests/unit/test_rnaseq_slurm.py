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

import polars as pl
import pytest

from pathlib import Path

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


def _shard_report(total, mapped, wall, rss, *, length=100, nbytes=1000):
    return {"organism": "human", "total_reads": total, "mapped_reads": mapped,
            "per_locus": {"IGH": mapped}, "constant_only_fragments": 1, "isotype_from_mate": 2,
            "min_score": 75.0, "threads": 8, "wall_seconds": wall, "peak_rss_mb": rss,
            "paired": True, "input": "r1.fq", "input_bytes": nbytes,
            "read_length_min": length, "read_length_max": length,
            "read_length_mean": float(length),
            "prefilter_stats": {"seen": total, "passed": mapped}, "segment_search": {}}


def test_merge_map_reports_sums_counts_and_never_fakes_a_single_wall_time():
    shards = [_shard_report(100, 10, 5.0, 300.0), _shard_report(100, 20, 9.0, 280.0)]

    m = pipeline._merge_map_reports(shards)
    assert m["shards"] == 2 and m["read_groups"] == 2
    assert m["total_reads"] == 200 and m["mapped_reads"] == 30
    assert m["per_locus"] == {"IGH": 30}
    assert m["constant_only_fragments"] == 2 and m["isotype_from_mate"] == 4
    assert m["wall_seconds_max"] == 9.0 and m["wall_seconds_sum"] == 14.0
    assert m["peak_rss_mb_max"] == 300.0
    # Summing 40 array tasks' wall time and calling it "wall_seconds" would be a lie.
    assert "wall_seconds" not in m and "peak_rss_mb" not in m


def test_merge_map_reports_keeps_the_library_shape_and_the_prefilter_accounting():
    """Never: these survive the merge or a sharded run silently describes a different library.

    `paired`, `input_bytes` and the read lengths are recorded by Stage 1 because nothing
    downstream can recover them -- the AIRR holds only the reads that mapped, so its row count and
    its sequence lengths describe the receptor subset. And `passed / seen` is the only number that
    says whether the prefilter earned its keep; dropping it reported nothing at all.
    """
    shards = [_shard_report(100, 10, 5.0, 300.0, length=100, nbytes=1000),
              _shard_report(300, 20, 9.0, 280.0, length=150, nbytes=4000)]
    m = pipeline._merge_map_reports(shards)
    assert m["paired"] is True
    assert m["input"] == ["r1.fq", "r1.fq"]
    assert m["input_bytes"] == 5000
    assert m["read_length_min"] == 100 and m["read_length_max"] == 150
    # Weighted by reads, not by shard: 100 reads at 100 nt and 300 at 150 nt is 137.5, not 125.
    assert m["read_length_mean"] == 137.5
    assert m["prefilter_stats"] == {"seen": 400, "passed": 30}
    assert m["segment_search"] == {}


def test_segment_search_is_merged_not_summed_because_it_carries_a_ratio():
    """Never: `fast_fraction` is a RATIO and `reasons` is nested; `segment_search` is not a counter dict.

    Summed over four read groups, a library whose true fast fraction is 0.1717 would report
    0.6868 -- and that is the number the regime rule is read off. It would say "these reads span V
    into J, use the amplicon preset" about a library where they do not. (Summing it also just
    crashes on the nested `reasons`, which is how this was found.)
    """
    def seg(implied, rescued, v_only):
        n = implied + rescued
        return {"implied": implied, "rescued": rescued, "no_segment_hit": 7,
                "v_only_on_segment": 12, "reasons": {"v_only": v_only, "j_only": 3},
                "fast_fraction": round(implied / n, 4)}

    shards = [dict(_shard_report(100, 10, 5.0, 300.0), segment_search=seg(17, 83, 40)),
              dict(_shard_report(100, 20, 9.0, 280.0), segment_search=seg(25, 75, 60))]
    s = pipeline._merge_map_reports(shards)["segment_search"]
    assert s["implied"] == 42 and s["rescued"] == 158
    assert s["no_segment_hit"] == 14 and s["v_only_on_segment"] == 24
    assert s["reasons"] == {"v_only": 100, "j_only": 6}
    # Recomputed from the totals: 42 / (42 + 158), NOT 0.17 + 0.25.
    assert s["fast_fraction"] == 0.21


def test_segment_search_stays_empty_when_the_two_pass_never_ran():
    shards = [_shard_report(100, 10, 5.0, 300.0), _shard_report(100, 20, 9.0, 280.0)]
    assert pipeline._merge_map_reports(shards)["segment_search"] == {}


def test_a_shard_that_mapped_nothing_does_not_drag_the_read_length_to_zero():
    """An empty shard has a wall time but no reads; its zeroed length fields are not a measurement."""
    shards = [_shard_report(100, 10, 5.0, 300.0, length=100),
              _shard_report(0, 0, 0.5, 60.0, length=0, nbytes=0)]
    m = pipeline._merge_map_reports(shards)
    assert m["read_length_min"] == 100 and m["read_length_mean"] == 100.0


def test_submit_script_runs_stage23_once_and_never_in_the_array(tmp_path):
    s = render_rnaseq_submit_script("/d/r1.fq.gz", "SAMP", tmp_path, shards=8,
                                    r2="/d/r2.fq.gz", out_dir="/o", partition="medium")
    assert s.count('"$ARDA" cluster reduce') == 1
    # The whole point: these must not appear as their own array steps.
    assert "arda rnaseq correct" not in s
    assert "arda rnaseq assemble" not in s
    assert '"$ARDA" map' in s and "--array=0-7" in s
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


# ---------------------------------------------------------- the multi-sample / read-group path

def _samples(tmp_path, spec):
    """Build real FASTQs and resolve them, so `plan` sees what the CLI would hand it."""
    from arda.samples import load

    r1, r2, ids = [], [], []
    for sid, n in spec:
        for i in range(n):
            a = tmp_path / f"{sid}_{i}_1.fq"
            b = tmp_path / f"{sid}_{i}_2.fq"
            a.write_text("@r\nACGT\n+\nIIII\n")
            b.write_text("@r\nACGT\n+\nIIII\n")
            r1.append(a)
            r2.append(b)
            ids.append(sid)
    return load(r1=r1, r2=r2, ids=ids)


def test_plan_writes_one_row_per_read_group_and_one_per_sample(tmp_path):
    from arda.cluster import READGROUP_COLUMNS, SAMPLE_COLUMNS, plan

    samples = _samples(tmp_path, [("A", 3), ("B", 1)])
    rg_path, s_path = plan(samples, tmp_path / "work")

    rg = [r.split("\t") for r in rg_path.read_text().splitlines()]
    assert tuple(rg[0]) == READGROUP_COLUMNS
    assert len(rg) == 1 + 4                       # 3 read groups for A, 1 for B
    assert [r[0] for r in rg[1:]] == ["A", "A", "A", "B"]
    assert [r[1] for r in rg[1:]] == ["0", "1", "2", "0"]

    sm = [r.split("\t") for r in s_path.read_text().splitlines()]
    assert tuple(sm[0]) == SAMPLE_COLUMNS
    assert [r[0] for r in sm[1:]] == ["A", "B"]
    # Per-sample shard dirs: two samples in one work dir must not collide on part names.
    assert len({r[1] for r in sm[1:]}) == 2
    for r in sm[1:]:
        assert Path(r[1]).is_dir()


def test_plan_part_names_sort_numerically_within_a_sample(tmp_path):
    """`reduce` merges parts in NAME order and that order IS read order."""
    from arda.cluster import plan

    samples = _samples(tmp_path, [("A", 12)])
    rg_path, _ = plan(samples, tmp_path / "work")
    parts = [Path(r.split("\t")[4]).name for r in rg_path.read_text().splitlines()[1:]]
    assert parts == sorted(parts), "shard_10 must not sort before shard_2"
    assert parts[0] == "shard_00000.airr.tsv" and parts[-1] == "shard_00011.airr.tsv"


def test_regime_flags_wires_stage1_to_the_denoising_preset():
    """Never: `--ec-mode rnaseq` reads a column Stage 1 only writes when asked for it.

    The mode commands make that link themselves because both stages happen in one call. A path
    that runs `map` and `reduce` as separate jobs cannot, so the flags must carry it -- otherwise
    the gate reads a column that is not there and silently does nothing.
    """
    from arda.cluster import regime_flags

    m, r = regime_flags("rnaseq", threads=4)
    assert "--prefilter" in m and "--two-pass" not in m
    assert "--junction-quality" in m
    assert "--ec-mode rnaseq" in r

    m, r = regime_flags("amplicon", threads=4)
    assert "--two-pass --fast-segments --v-only-on-segment" in m and "--prefilter" not in m
    assert "--ec-mode amplicon" in r

    # `fast` is the one preset that reads no quality column, so it must not demand one.
    m, _ = regime_flags("rnaseq", ec_mode="fast")
    assert "--junction-quality" not in m

    with pytest.raises(ValueError, match="regime must be one of"):
        regime_flags("bulk")


def test_samples_submit_arrays_over_read_groups_then_reduces_per_sample(tmp_path):
    from arda.cluster import plan, render_samples_submit_script

    samples = _samples(tmp_path, [("A", 3), ("B", 1)])
    rg_path, s_path = plan(samples, tmp_path / "work")
    s = render_samples_submit_script(rg_path, s_path, n_read_groups=4, n_samples=2,
                                     out_dir="/o", partition="medium", regime="rnaseq")

    # The map array is sized by READ GROUPS across every sample, the reduce array by samples.
    assert "--array=0-3" in s and "--job-name=arda-map" in s
    assert "--array=0-1" in s and "--job-name=arda-reduce" in s
    # One reduce per sample, gated on the WHOLE map array -- not `aftercorr`, whose
    # index-for-index pairing is wrong when one sample consumes many map tasks.
    assert "--dependency=afterok:$MAP_JID" in s
    assert "aftercorr" not in s
    # No split step: a sample delivered as several files is already sharded.
    assert '"$ARDA" cluster split ' not in s
    # Stages 2-3 appear exactly once, and never inside the map array.
    assert s.count('"$ARDA" cluster reduce') == 1
    assert "--junction-quality" in s and "--ec-mode rnaseq" in s


def test_samples_submit_reads_its_inputs_from_the_manifest_not_from_filenames(tmp_path):
    from arda.cluster import plan, render_samples_submit_script

    samples = _samples(tmp_path, [("A", 2)])
    rg_path, s_path = plan(samples, tmp_path / "work")
    s = render_samples_submit_script(rg_path, s_path, n_read_groups=2, n_samples=1)
    assert str(rg_path) in s and str(s_path) in s
    # Row i is line i+2: line 1 is the header and sed is 1-based.
    assert "+ 2 ))p" in s
    # Never: fields come out with `cut -f`, never with `IFS=$'\t' read`. The body is passed to
    # sbatch inside `--wrap '...'`, and a single-quoted string cannot contain `$'\t'` -- bash
    # closes the quote at the `$` and it degrades to `IFS=$t`. Under the script's own `set -u`
    # that is an unbound variable and the task dies; without it, default IFS collapses a
    # single-end row's blank `r2` and every later field shifts left.
    assert "IFS=" not in s
    assert "cut -f3" in s and "cut -f4" in s
    # A blank r2 must drop the flag, not pass an empty path.
    assert '${r2:+--r2 "$r2"}' in s


def test_merge_map_reports_carries_a_throughput_number_like_a_single_node_run():
    """A cohort mixing one-file and many-file samples needs ONE throughput column, not two.

    Stage 1 reports `reads_per_second` for a single pair. The merge reported none at all, so in a
    cohort where some samples arrive as one pair and some as four lanes the QC table had a
    half-empty column -- and the missing half was exactly the samples big enough to be split.
    The denominator is the SUMMED per-shard wall, the same quantity the single-node number
    divides by; `wall_seconds_max` would be the rate of a fan-out this function cannot know about.
    """
    shards = [_shard_report(100, 10, 5.0, 300.0), _shard_report(200, 20, 5.0, 280.0)]
    m = pipeline._merge_map_reports(shards)
    assert m["wall_seconds_sum"] == 10.0
    assert m["reads_per_second"] == 30.0           # 300 reads / 10.0 s
