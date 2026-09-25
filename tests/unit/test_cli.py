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
    assert '"$ARDA" cluster split-fasta' in body and "--array=0-3" in body
    assert '"$ARDA" cluster merge' in body


# ── the 2.16.0 mode surface ───────────────────────────────────────────────────────────────────

def test_rnaseq_run_is_gone(monkeypatch, tmp_path):
    """The hard break. `arda rnaseq run` was the amplicon entry point too; it must not resolve.

    `rnaseq` is now a COMMAND, not a group, so `run` arrives as a stray positional argument
    rather than a subcommand — and a command with no arguments must reject it instead of
    quietly ignoring it and running the pipeline anyway.
    """
    monkeypatch.setattr("arda.rnaseq.pipeline.run", lambda pairs, **kw: None)
    res = runner.invoke(app, ["rnaseq", "run", "--r1", "r1.fq", "-p", "S", "-d", str(tmp_path)])
    assert res.exit_code != 0



@pytest.fixture
def r1(tmp_path):
    """A real FASTQ: the mode commands check their inputs exist before the pipeline starts."""
    p = tmp_path / "r1.fq"
    p.write_text("@r\nACGT\n+\nIIII\n")
    return str(p)

def test_singlecell_is_reserved_and_points_at_the_command_that_exists():
    # The MODE name stays reserved -- its only sensible preset is what `--exact` already gives.
    # The single-cell WORK is `arda cells`, so refusing without naming it is the useless half.
    res = runner.invoke(app, ["singlecell"])
    assert res.exit_code == 2
    out = res.stdout + str(res.stderr)
    assert "reserved" in out
    assert "arda cells" in out


@pytest.mark.parametrize("mode,expected", [
    ("amplicon", {"two_pass": True, "fast_segments": True, "segment_only_v": True,
                  "prefilter": False}),
    ("rnaseq", {"two_pass": False, "fast_segments": False, "segment_only_v": False,
                "prefilter": True}),
])
def test_mode_presets_are_the_measured_configurations(mode, expected):
    """Never: The two configurations do NOT compose, and `--two-pass` alone is a loss on both regimes.

    Pinning the table is the point: if `rnaseq` ever gains `two_pass` without `fast_segments`, it
    silently ships the dominated config (0.762x on bulk) under a name that promises the opposite.
    """
    assert _MODE_SPEED[mode] == expected


def test_exact_clears_every_speedup(monkeypatch, tmp_path, r1):
    seen = {}

    def fake_run(pairs, **kw):
        seen.update(kw, pairs=pairs)

    monkeypatch.setattr("arda.rnaseq.pipeline.run", fake_run)
    res = runner.invoke(app, ["amplicon", "--r1", r1, "-p", "S", "-d", str(tmp_path),
                              "--exact"])
    assert res.exit_code == 0, res.output
    assert not any(seen[k] for k in _MODE_SPEED["amplicon"])
    assert seen["pairs"] == [(Path(r1), None)]


def test_mode_passes_its_preset_through(monkeypatch, tmp_path, r1):
    seen = {}
    monkeypatch.setattr("arda.rnaseq.pipeline.run",
                        lambda pairs, **kw: seen.update(kw, pairs=pairs))
    res = runner.invoke(app, ["rnaseq", "--r1", r1, "-p", "S", "-d", str(tmp_path)])
    assert res.exit_code == 0, res.output
    assert seen["prefilter"] is True and seen["fast_segments"] is False
    # The mode's own denoising default, not the historical `fast`.
    assert seen["ec_mode"] == "rnaseq"
    assert seen["shm"] == "framework"


def test_indel_rescue_without_fast_segments_raises(monkeypatch, tmp_path, r1):
    """Never: A flag that is accepted and silently does nothing is the failure this project keeps
    hitting. `--indel-rescue` needs the fast segment pass, so `--exact` must reject it.

    Asserted on BEHAVIOUR — non-zero exit, and the pipeline never started — not on the message.
    typer renders errors through rich, in a box wrapped to the terminal width, so a substring
    assertion over that passes on an 80-column laptop and fails on a CI runner that wraps
    `--indel-rescue` across the line break. It did exactly that.
    """
    called = []
    monkeypatch.setattr("arda.rnaseq.pipeline.run",
                        lambda pairs, **kw: called.append(kw))
    res = runner.invoke(app, ["amplicon", "--r1", r1, "-p", "S", "-d", str(tmp_path),
                              "--indel-rescue", "--exact"])
    assert res.exit_code != 0
    assert not called, "the pipeline ran despite an unsatisfiable flag combination"


def test_the_two_version_literals_agree():
    """Never: `arda.__version__` and `pyproject.toml`'s `version` are TWO literals with no link.

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


# ------------------------------------------------- read-group scheduling surface
#
# `cluster plan` and `cluster submit-samples` are how a foreign scheduler drives arda at READ
# GROUP granularity. `cluster.plan` and the script renderer are unit-tested; the CLI wiring
# around them -- the regime check, sheet loading, and the recipe it prints -- was not, and that
# wiring is where a typo reaches a stakeholder's job array.


@pytest.fixture
def two_sample_sheet(tmp_path):
    (tmp_path / "a_1.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "a_2.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "b_1.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "b_2.fq").write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\n"
                     "A\ta_1.fq\ta_2.fq\n"
                     "A\tb_1.fq\tb_2.fq\n"
                     "B\tb_1.fq\tb_2.fq\n")
    return sheet


def test_cluster_plan_writes_both_manifests_and_prints_the_flags_filled_in(two_sample_sheet,
                                                                          tmp_path):
    """Never: the printed recipe carries REAL flags, not `<flags>`.

    `arda rnaseq` wires Stage 1 to Stage 2 itself -- `--ec-mode rnaseq` reads a column Stage 1
    writes only under `--junction-quality`. A scheduler running the two as separate jobs cannot
    know that, so `plan` prints both halves resolved from the regime.
    """
    work = tmp_path / "work"
    res = runner.invoke(app, ["cluster", "plan", "--samples", str(two_sample_sheet),
                              "--work-dir", str(work), "-d", str(tmp_path / "res"),
                              "--regime", "rnaseq", "--threads", "16"])
    assert res.exit_code == 0, res.output
    rg = (work / "readgroups.tsv").read_text().splitlines()
    assert rg[0].split("\t") == ["sample", "read_group", "r1", "r2", "part_airr", "part_report"]
    assert len(rg) == 4                                   # header + 3 read groups
    assert [ln.split("\t")[0] for ln in rg[1:]] == ["A", "A", "B"]
    assert [ln.split("\t")[1] for ln in rg[1:]] == ["0", "1", "0"]   # per-sample, from zero
    assert (work / "samples.tsv").read_text().splitlines()[1].startswith("A\t")
    assert "2 sample(s), 3 read group(s)" in res.output
    assert "--junction-quality" in res.output and "--ec-mode rnaseq" in res.output
    assert "<flags>" not in res.output


def test_cluster_plan_refuses_a_regime_it_has_no_preset_for(two_sample_sheet, tmp_path):
    res = runner.invoke(app, ["cluster", "plan", "--samples", str(two_sample_sheet),
                              "--work-dir", str(tmp_path / "w"), "--regime", "bulk"])
    assert res.exit_code != 0
    assert "rnaseq" in res.output and "amplicon" in res.output


def test_cluster_submit_samples_renders_two_arrays_without_submitting(two_sample_sheet, tmp_path):
    """No `--submit`, no sbatch: the script is written and the sizes come from the sheet."""
    work = tmp_path / "work"
    res = runner.invoke(app, ["cluster", "submit-samples", "--samples", str(two_sample_sheet),
                              "--work-dir", str(work), "-d", str(tmp_path / "res"),
                              "--regime", "amplicon", "--threads", "8",
                              "--partition", "medium"])
    assert res.exit_code == 0, res.output
    script = (work / "submit.sh").read_text()
    assert "--array=0-2" in script                        # 3 read groups
    assert "--array=0-1" in script                        # 2 samples
    assert "--dependency=afterok:$MAP_JID" in script      # not aftercorr: one sample, many tasks
    assert "--partition=medium" in script
    assert '"$ARDA" map' in script and '"$ARDA" cluster reduce' in script
    assert "IFS=" not in script and "cut -f3" in script   # see render_samples_submit_script


def test_markup_d_prior_implies_the_posterior_rather_than_being_accepted_and_ignored(tmp_path):
    """Never: a parameter accepted while doing nothing is this project's recurring failure.

    `--d-prior` on its own used to be unrepresentable -- the prior columns were gated on
    `--d-posterior` alone -- so a user who passed a table they had just fitted would get a file
    with no D columns in it and exit 0.
    """
    records = tmp_path / "in.tsv"
    records.write_text("cdr3\tv\tj\tspecies\n"
                       "CASSLAPGATNEKLFF\tTRBV5-1*01\tTRBJ1-4*01\thuman\n")
    prior = tmp_path / "prior.tsv"
    rows = ["locus\tkind\tkey\tvalue"]
    for i in range(16):
        rows += [f"TRB\tinsVD\t{i}\t{1 / 16:.6f}", f"TRB\tinsDJ\t{i}\t{1 / 16:.6f}"]
    for allele in ("TRBD1*01", "TRBD2*01"):
        rows += [f"TRB\tdlen\t{allele}:{n}\t{1 / 12:.6f}" for n in range(1, 13)]
        rows.append(f"TRB\td_marginal\t{allele}\t0.5")
    rows.append("TRB\tbeta\tbeta\t1.25")
    prior.write_text("\n".join(rows) + "\n")

    out = tmp_path / "out.tsv"
    result = runner.invoke(app, ["markup", "-i", str(records), "-o", str(out),
                                 "--d-prior", str(prior)])
    assert result.exit_code == 0, result.output
    assert "d_call" in out.read_text().splitlines()[0]


def test_markup_d_prior_refuses_a_path_that_is_not_there(tmp_path):
    records = tmp_path / "in.tsv"
    records.write_text("cdr3\tv\tj\tspecies\n"
                       "CASSLAPGATNEKLFF\tTRBV5-1*01\tTRBJ1-4*01\thuman\n")
    result = runner.invoke(app, ["markup", "-i", str(records), "-o", str(tmp_path / "o.tsv"),
                                 "--d-prior", str(tmp_path / "absent.tsv")])
    assert result.exit_code != 0
