"""Resolve one search-thread budget within the process and SLURM allocation."""

import os
from pathlib import Path


def _cgroup_cpus() -> int | None:
    """Read v2/v1 CPU quotas, including parent limits visible in the cgroup mount."""
    root = Path("/sys/fs/cgroup")
    bases = [root, root / "cpu", root / "cpu,cpuacct", root / "cpuacct,cpu"]
    paths = list(bases)
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            _, controllers, relative = line.split(":", 2)
            if not controllers or "cpu" in controllers.split(","):
                # Controller order in /proc need not match the mount directory.
                for base in bases if controllers else [root]:
                    path = base / relative.lstrip("/")
                    if ".." not in path.parts:
                        paths += [path, *[p for p in path.parents if p == base or base in p.parents]]
    except OSError:
        pass
    limits = []
    for path in set(paths):
        try:
            quota, period = (path / "cpu.max").read_text().split()
        except OSError:
            try:
                quota = (path / "cpu.cfs_quota_us").read_text().strip()
                period = (path / "cpu.cfs_period_us").read_text().strip()
            except OSError:
                continue
        if quota != "max" and int(quota) > 0 and int(period) > 0:
            limits.append(max(1, int(quota) // int(period)))
    return min(limits) if limits else None


def streaming_threads(threads: int, table_threads: int) -> tuple[int, bool]:
    """Reserve a reader CPU only when both table and search phases can spare it."""
    read_ahead = threads > max(1, table_threads)
    return threads - int(read_ahead), read_ahead


def import_polars():
    """Initialize ARDA's table pool serially without changing caller environment."""
    previous = os.environ.get("POLARS_MAX_THREADS")
    os.environ.setdefault("POLARS_MAX_THREADS", "1")
    try:
        import polars as pl

        pl.thread_pool_size()  # force initialization before restoring the environment
        return pl
    finally:
        if previous is None:
            os.environ.pop("POLARS_MAX_THREADS", None)
        else:
            os.environ["POLARS_MAX_THREADS"] = previous


def resolve_threads(threads: int = 0) -> int:
    """Zero selects all allocated CPUs; reject requests outside that allocation."""
    available = os.cpu_count() or 1
    if hasattr(os, "sched_getaffinity"):
        available = min(available, len(os.sched_getaffinity(0)))
    if quota := _cgroup_cpus():
        available = min(available, quota)
    allocated = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated is not None:
        try:
            limit = int(allocated)
        except ValueError:
            raise ValueError("SLURM_CPUS_PER_TASK must be a positive integer") from None
        if limit < 1:
            raise ValueError("SLURM_CPUS_PER_TASK must be a positive integer")
        available = min(available, limit)
    if not isinstance(threads, int) or isinstance(threads, bool) or not 0 <= threads <= available:
        raise ValueError(f"threads must be 0 (all allocated CPUs) or 1..{available}, got {threads!r}")
    return threads or available
