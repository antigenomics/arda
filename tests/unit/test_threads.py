"""Automatic mapping must respect the CPUs granted to this process/job."""

import pytest

from arda import _threads


def test_all_uses_smallest_process_and_scheduler_limit(monkeypatch):
    monkeypatch.setattr(_threads, "_cgroup_cpus", lambda: None)
    monkeypatch.setattr(_threads.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(_threads.os, "sched_getaffinity", lambda _: set(range(12)), raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    assert _threads.resolve_threads(0) == 4
    assert _threads.resolve_threads(2) == 2
    for invalid in (-1, 5):
        with pytest.raises(ValueError, match="threads"):
            _threads.resolve_threads(invalid)


def test_without_affinity_or_scheduler_uses_host_count(monkeypatch):
    monkeypatch.setattr(_threads, "_cgroup_cpus", lambda: None)
    monkeypatch.delattr(_threads.os, "sched_getaffinity", raising=False)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setattr(_threads.os, "cpu_count", lambda: 8)
    assert _threads.resolve_threads(0) == 8


def test_invalid_scheduler_allocation_fails(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "invalid")
    with pytest.raises(ValueError, match="SLURM_CPUS_PER_TASK"):
        _threads.resolve_threads(0)


def test_package_import_does_not_initialize_a_host_sized_table_pool():
    import subprocess
    import sys

    result = subprocess.run([sys.executable, "-c",
                             "import sys, arda; assert 'polars' not in sys.modules"],
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("version", [1, 2])
def test_quota_only_container_limits_threads(tmp_path, monkeypatch, version):
    monkeypatch.setattr(_threads, "Path", lambda p: tmp_path / p.lstrip("/"))
    monkeypatch.setattr(_threads.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(_threads.os, "sched_getaffinity", lambda _: set(range(64)), raising=False)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    root = tmp_path / "sys/fs/cgroup"
    if version == 2:
        root.mkdir(parents=True)
        (root / "cpu.max").write_text("200000 100000\n")
    else:
        root = root / "cpu"
        root.mkdir(parents=True)
        (root / "cpu.cfs_quota_us").write_text("200000\n")
        (root / "cpu.cfs_period_us").write_text("100000\n")
    assert _threads.resolve_threads(0) == 2
    with pytest.raises(ValueError, match="threads"):
        _threads.resolve_threads(3)


def test_nested_v1_quota_does_not_depend_on_controller_order(tmp_path, monkeypatch):
    monkeypatch.setattr(_threads, "Path", lambda p: tmp_path / p.lstrip("/"))
    proc = tmp_path / "proc/self"
    proc.mkdir(parents=True)
    (proc / "cgroup").write_text("2:cpuacct,cpu:/workload\n")
    mount = tmp_path / "sys/fs/cgroup/cpu,cpuacct"
    child = mount / "workload"
    child.mkdir(parents=True)
    (mount / "cpu.cfs_quota_us").write_text("-1")
    (mount / "cpu.cfs_period_us").write_text("100000")
    (child / "cpu.cfs_quota_us").write_text("200000")
    (child / "cpu.cfs_period_us").write_text("100000")
    assert _threads._cgroup_cpus() == 2


def test_streaming_keeps_input_synchronous_when_table_pool_fills_budget(tmp_path, monkeypatch):
    import threading
    from arda.annotate import mapper

    caller = threading.get_ident()
    monkeypatch.setattr(mapper, "_prep", lambda *args: (None, None, 2, 7, 2))
    monkeypatch.setattr(mapper.pl, "thread_pool_size", lambda: 2)

    def reads(_):
        assert threading.get_ident() == caller
        yield "q", "ACGT"

    def annotate(records, *args, **kwargs):
        assert kwargs["threads"] == 2
        return []

    monkeypatch.setattr(mapper.seqio, "read_sequences", reads)
    monkeypatch.setattr(mapper, "_annotate_chunk", annotate)
    mapper.annotate_file("unused", tmp_path / "out.tsv", threads=2)


def test_library_rejects_an_existing_oversized_table_pool():
    import os
    import subprocess
    import sys

    program = '''
import polars as pl
assert pl.thread_pool_size() == 4
import arda
try:
    arda.annotate_sequences([], threads=1)
except ValueError as exc:
    assert "Polars is already initialized" in str(exc)
else:
    raise AssertionError("Oversized existing table pool was accepted")
'''
    result = subprocess.run([sys.executable, "-c", program],
                            env=dict(os.environ, POLARS_MAX_THREADS="4"),
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
