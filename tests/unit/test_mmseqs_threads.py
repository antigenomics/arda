"""Every threaded MMseqs phase must have an explicit bounded thread count."""

import subprocess

import pytest

from arda import mmseqs


@pytest.mark.parametrize("command", ["createdb", "filterdb", "createindex"])
def test_preparation_and_reduction_do_not_default_to_every_host_core(command, monkeypatch):
    calls = []
    monkeypatch.setattr(mmseqs, "mmseqs_binary", lambda: "mmseqs")
    monkeypatch.setattr(mmseqs.subprocess, "run", lambda cmd, **kw:
                        calls.append((cmd, kw)) or subprocess.CompletedProcess(cmd, 0))
    mmseqs.run([command, "in", "out"])
    cmd, kw = calls[0]
    assert cmd[-2:] == ["--threads", "1"]
    assert kw["env"]["OMP_THREAD_LIMIT"] == "1"


def test_search_keeps_budget_and_disables_nested_openmp(monkeypatch):
    calls = []
    monkeypatch.setattr(mmseqs, "mmseqs_binary", lambda: "mmseqs")
    monkeypatch.setattr(mmseqs.subprocess, "run", lambda cmd, **kw:
                        calls.append((cmd, kw)) or subprocess.CompletedProcess(cmd, 0))
    monkeypatch.setenv("OMP_NUM_THREADS", "1000")
    mmseqs.run(["search", "q", "t", "r", "tmp", "--threads", "2"])
    cmd, kw = calls[0]
    assert cmd.count("--threads") == 1
    assert kw["env"]["OMP_NUM_THREADS"] == "2"
    assert kw["env"]["OMP_MAX_ACTIVE_LEVELS"] == "1"
    assert mmseqs.os.environ["OMP_NUM_THREADS"] == "1000"
