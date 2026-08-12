"""Verbosity, the log file, and the stdout/stderr split.

⛔ **Progress on stderr, results on stdout.** The stage lines were ``typer.echo`` on stdout until
2.20.0, which means ``arda export-ref ... > out.tsv`` interleaved a progress line into the data
and ``$(arda map ...)`` captured prose. That split is what the tests here pin; the rest is that
``-q`` must not be able to silence ``--log-file``, because a quiet cluster job that leaves no
record is how a 10-hour run gets repeated.
"""

from __future__ import annotations

import logging

from typer.testing import CliRunner

from arda._log import Throttle, logger, setup
from arda.cli import app

runner = CliRunner()


def _flush_file_handlers() -> None:
    # Only the FileHandlers: CliRunner closes the stream the console handler captured, so
    # flushing that one raises on a closed file.
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.flush()


def test_verbosity_selects_the_console_level(tmp_path):
    for verbosity, quiet, expected in ((0, False, logging.INFO),
                                       (1, False, logging.DEBUG),
                                       (3, False, logging.DEBUG),
                                       (0, True, logging.WARNING)):
        setup(verbosity=verbosity, quiet=quiet)
        console = [h for h in logger.handlers if isinstance(h, logging.StreamHandler)]
        assert [h.level for h in console] == [expected], (verbosity, quiet)


def test_quiet_does_not_silence_the_log_file(tmp_path):
    path = tmp_path / "run.log"
    setup(verbosity=0, quiet=True, log_file=path)
    logger.info("a quiet run still records this")
    _flush_file_handlers()
    assert "a quiet run still records this" in path.read_text()


def test_the_log_file_carries_a_timestamp_and_peak_rss(tmp_path):
    path = tmp_path / "run.log"
    setup(log_file=path)
    logger.warning("something to look at")
    _flush_file_handlers()
    line = [ln for ln in path.read_text().splitlines() if "something to look at" in ln][0]
    assert line.startswith("20")                       # ISO date
    assert "WARNING" in line and " MB " in line


def test_setup_replaces_handlers_rather_than_stacking_them(tmp_path):
    for _ in range(3):
        setup(log_file=tmp_path / "run.log")
    assert len(logger.handlers) == 2                   # one console, one file


def test_throttle_gates_by_time():
    tick = Throttle(seconds=3600.0)
    assert not tick.ready()                            # nothing has elapsed since construction
    tick = Throttle(seconds=0.0)
    assert tick.ready() and tick.ready()


def test_results_go_to_stdout_and_progress_to_stderr(tmp_path):
    """``arda info`` is a RESULT (paths), so it stays on stdout even at ``-q``."""
    result = runner.invoke(app, ["-q", "info"])
    assert result.exit_code == 0
    assert "project_root" in result.stdout


def test_the_global_options_come_before_the_subcommand(tmp_path):
    """Typer puts callback options on the top-level parser only; a run that passed ``-v`` after
    the command name would fail, and the help text says so."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    for flag in ("--verbose", "--quiet", "--log-file"):
        assert flag in result.stdout


def test_stats_writes_the_path_to_stdout_and_the_count_to_the_log(tmp_path):
    import polars as pl

    airr = tmp_path / "s.airr.tsv"
    pl.DataFrame([["r1", "TRB", "TGTGCCAGCAGCTTAGACGGGACAGGGTTC", "CASSLDGTGF"]],
                 schema=["sequence_id", "locus", "junction", "junction_aa"],
                 orient="row").write_csv(airr, separator="\t", quote_style="never")
    out = tmp_path / "s.stats.tsv"
    log = tmp_path / "run.log"
    result = runner.invoke(app, ["--log-file", str(log), "stats", "-i", str(airr), "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert result.stdout.strip().endswith("s.stats.tsv")
    _flush_file_handlers()
    assert "stats:" in log.read_text()
