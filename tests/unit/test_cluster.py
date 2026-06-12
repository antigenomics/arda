"""Unit tests for SLURM sharding (split/merge/script) — no cluster required."""

from pathlib import Path

from arda.cluster import split, merge, render_submit_script
from arda.annotate.io import read_sequences


def _write_fasta(path: Path, n: int):
    path.write_text("".join(f">s{i}\nACGT{i:04d}AAA\n" for i in range(n)))


def test_split_partitions_all_records_disjointly(tmp_path):
    src = tmp_path / "in.fasta"
    _write_fasta(src, 53)                      # not a multiple of shard count
    paths = split(src, tmp_path / "shards", shards=8)
    assert len(paths) == 8
    seen = []
    for p in paths:
        ids = [sid for sid, _ in read_sequences(p)]
        seen.extend(ids)
    # Every record appears exactly once across all shards.
    assert sorted(seen) == sorted(f"s{i}" for i in range(53))
    assert len(seen) == 53
    # Round-robin balance: shard sizes differ by at most 1.
    sizes = [sum(1 for _ in read_sequences(p)) for p in paths]
    assert max(sizes) - min(sizes) <= 1


def test_merge_keeps_single_header(tmp_path):
    d = tmp_path / "out"
    d.mkdir()
    (d / "out_0.tsv").write_text("sequence_id\tv_call\nq0\tIGHV1\n")
    (d / "out_1.tsv").write_text("sequence_id\tv_call\nq1\tIGHV2\nq2\tIGHV3\n")
    combined = merge(d, tmp_path / "all.tsv")
    lines = combined.read_text().splitlines()
    assert lines[0] == "sequence_id\tv_call"          # header once
    assert lines.count("sequence_id\tv_call") == 1
    assert lines[1:] == ["q0\tIGHV1", "q1\tIGHV2", "q2\tIGHV3"]


def test_render_submit_script_chains_steps(tmp_path):
    s = render_submit_script(
        "big.fastq", "out.airr.tsv", tmp_path / "work",
        shards=50, organism="human", seqtype="nt", threads=8, partition="cpu")
    assert "arda split" in s
    assert "--array=0-49" in s
    assert "arda annotate" in s and "shard_${SLURM_ARRAY_TASK_ID}.fasta" in s
    assert "--dependency=afterok:$ARRAY_JID" in s
    assert "arda merge" in s
    assert "--partition=cpu" in s and "--cpus-per-task=8" in s
