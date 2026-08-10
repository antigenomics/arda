"""Unit tests for the Typer CLI surface — no external tools or cluster required.

These exercise the CLI layer end-to-end via ``CliRunner`` for the commands that delegate to
already-tested pure helpers (``info``, ``cluster split-fasta``/``merge``/``submit-fasta``), plus
the 2.16.0 mode surface: three named modes, each carrying its own speed preset.
"""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from arda import __version__
from arda.cli import app, _MODE_SPEED
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

    res = runner.invoke(app, ["cluster", "split-fasta", str(src), str(shards_dir), "--shards", "3"])
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
    res = runner.invoke(app, ["cluster", "merge", str(out_dir), str(combined)])
    assert res.exit_code == 0
    lines = combined.read_text().splitlines()
    assert lines.count("sequence_id\tv_call") == 1
    assert lines[1:] == ["q0\tIGHV1", "q1\tIGHV2"]


def test_slurm_writes_executable_submit_script(tmp_path):
    work = tmp_path / "work"
    res = runner.invoke(
        app,
        ["cluster", "submit-fasta", "-i", "big.fastq", "-o", "out.airr.tsv",
         "--work-dir", str(work), "--shards", "4"],
    )
    assert res.exit_code == 0
    submit = work / "submit.sh"
    assert submit.exists()
    assert submit.stat().st_mode & 0o111  # executable
    body = submit.read_text()
    # The script must name commands that EXIST — the whole point of moving them into one group
    # is that a generated script cannot invoke a command the CLI no longer has.
    assert "arda cluster split-fasta" in body and "--array=0-3" in body
    assert "arda cluster merge" in body


# ── the 2.16.0 mode surface ───────────────────────────────────────────────────────────────────

def test_rnaseq_run_is_gone(monkeypatch, tmp_path):
    """The hard break. `arda rnaseq run` was the amplicon entry point too; it must not resolve.

    `rnaseq` is now a COMMAND, not a group, so `run` arrives as a stray positional argument
    rather than a subcommand — and a command with no arguments must reject it instead of
    quietly ignoring it and running the pipeline anyway.
    """
    monkeypatch.setattr("arda.rnaseq.pipeline.run", lambda **kw: None)
    res = runner.invoke(app, ["rnaseq", "run", "--r1", "r1.fq", "-p", "S", "-d", str(tmp_path)])
    assert res.exit_code != 0


def test_singlecell_is_reserved_and_refuses():
    res = runner.invoke(app, ["singlecell"])
    assert res.exit_code == 2
    assert "not implemented" in (res.stdout + str(res.stderr))


@pytest.mark.parametrize("mode,expected", [
    ("amplicon", {"two_pass": True, "fast_segments": True, "segment_only_v": True,
                  "prefilter": False}),
    ("rnaseq", {"two_pass": False, "fast_segments": False, "segment_only_v": False,
                "prefilter": True}),
])
def test_mode_presets_are_the_measured_configurations(mode, expected):
    """⛔ The two configurations do NOT compose, and `--two-pass` alone is a loss on both regimes.

    Pinning the table is the point: if `rnaseq` ever gains `two_pass` without `fast_segments`, it
    silently ships the dominated config (0.762x on bulk) under a name that promises the opposite.
    """
    assert _MODE_SPEED[mode] == expected


def test_exact_clears_every_speedup(monkeypatch, tmp_path):
    seen = {}

    def fake_run(**kw):
        seen.update(kw)

    monkeypatch.setattr("arda.rnaseq.pipeline.run", fake_run)
    res = runner.invoke(app, ["amplicon", "--r1", "r1.fq", "-p", "S", "-d", str(tmp_path),
                              "--exact"])
    assert res.exit_code == 0, res.output
    assert not any(seen[k] for k in _MODE_SPEED["amplicon"])


def test_mode_passes_its_preset_through(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr("arda.rnaseq.pipeline.run", lambda **kw: seen.update(kw))
    res = runner.invoke(app, ["rnaseq", "--r1", "r1.fq", "-p", "S", "-d", str(tmp_path)])
    assert res.exit_code == 0, res.output
    assert seen["prefilter"] is True and seen["fast_segments"] is False
    # The mode's own denoising default, not the historical `fast`.
    assert seen["ec_mode"] == "rnaseq"
    assert seen["shm"] == "framework"


def test_indel_rescue_without_fast_segments_raises(monkeypatch, tmp_path):
    """⛔ A flag that is accepted and silently does nothing is the failure this project keeps
    hitting. `--indel-rescue` needs the fast segment pass, so `--exact` must reject it.

    Asserted on BEHAVIOUR — non-zero exit, and the pipeline never started — not on the message.
    typer renders errors through rich, in a box wrapped to the terminal width, so a substring
    assertion over that passes on an 80-column laptop and fails on a CI runner that wraps
    `--indel-rescue` across the line break. It did exactly that.
    """
    called = []
    monkeypatch.setattr("arda.rnaseq.pipeline.run", lambda **kw: called.append(kw))
    res = runner.invoke(app, ["amplicon", "--r1", "r1.fq", "-p", "S", "-d", str(tmp_path),
                              "--indel-rescue", "--exact"])
    assert res.exit_code != 0
    assert not called, "the pipeline ran despite an unsatisfiable flag combination"


def test_the_two_version_literals_agree():
    """⛔ `arda.__version__` and `pyproject.toml`'s `version` are TWO literals with no link.

    `publish.yml` asserts pyproject == the release tag, and nothing asserted this one — so a
    release could ship with `arda --version` reporting the PREVIOUS release, and every
    `arda.json` provenance block would record the wrong version. Caught exactly that way: 2.16.0's
    pyproject was bumped and `src/arda/__init__.py` still said 2.15.0, which `setup.sh` printed
    without anything failing.
    """
    import re
    from pathlib import Path

    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    m = re.search(r'^version = "([^"]+)"', pyproject.read_text(), re.M)
    assert m, "no version in pyproject.toml"
    assert m.group(1) == __version__, (
        f"pyproject.toml says {m.group(1)}, arda.__version__ says {__version__}")
