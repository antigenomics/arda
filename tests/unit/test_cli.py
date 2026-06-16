"""Unit tests for the Typer CLI surface — no external tools or cluster required.

These exercise the CLI layer end-to-end via ``CliRunner`` for the commands that
delegate to already-tested pure helpers (``info``, ``split``/``merge``, ``slurm``).
"""

from pathlib import Path

from typer.testing import CliRunner

from arda import __version__
from arda.cli import app
from arda.annotate.io import read_sequences

runner = CliRunner()


def _write_fasta(path: Path, n: int):
    path.write_text("".join(f">s{i}\nACGT{i:04d}AAA\n" for i in range(n)))


def test_info_reports_version_and_paths():
    result = runner.invoke(app, ["info"])
    assert result.exit_code == 0
    assert __version__ in result.stdout
    assert "project_root" in result.stdout
    assert "database_dir" in result.stdout


def test_split_then_merge_roundtrip(tmp_path):
    src = tmp_path / "in.fasta"
    _write_fasta(src, 10)
    shards_dir = tmp_path / "shards"

    res = runner.invoke(app, ["split", str(src), str(shards_dir), "--shards", "3"])
    assert res.exit_code == 0
    assert "wrote 3 shards" in res.stdout
    shard_files = sorted(shards_dir.glob("*.fasta"))
    assert len(shard_files) == 3
    # Every record lands in exactly one shard.
    seen = [sid for p in shard_files for sid, _ in read_sequences(p)]
    assert sorted(seen) == sorted(f"s{i}" for i in range(10))

    # merge fake per-shard AIRR TSVs back into one with a single header.
    out_dir = tmp_path / "airr"
    out_dir.mkdir()
    (out_dir / "out_0.tsv").write_text("sequence_id\tv_call\nq0\tIGHV1\n")
    (out_dir / "out_1.tsv").write_text("sequence_id\tv_call\nq1\tIGHV2\n")
    combined = tmp_path / "all.tsv"
    res = runner.invoke(app, ["merge", str(out_dir), str(combined)])
    assert res.exit_code == 0
    lines = combined.read_text().splitlines()
    assert lines.count("sequence_id\tv_call") == 1
    assert lines[1:] == ["q0\tIGHV1", "q1\tIGHV2"]


def test_slurm_writes_executable_submit_script(tmp_path):
    work = tmp_path / "work"
    res = runner.invoke(
        app,
        ["slurm", "-i", "big.fastq", "-o", "out.airr.tsv",
         "--work-dir", str(work), "--shards", "4"],
    )
    assert res.exit_code == 0
    submit = work / "submit.sh"
    assert submit.exists()
    assert submit.stat().st_mode & 0o111  # executable
    body = submit.read_text()
    assert "arda split" in body and "--array=0-3" in body and "arda merge" in body
