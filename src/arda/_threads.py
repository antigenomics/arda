"""Resolve one search-thread budget within the process and SLURM allocation."""

import os
import subprocess
import sys
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


def _windows_cpu_counts():
    """Return process-affinity logical/physical counts using Windows topology.

    Multiple processor groups cannot be represented by GetProcessAffinityMask;
    reject them instead of silently counting only the primary group.
    The current job rate is honored; Windows does not expose parent-job rates here.
    """
    import ctypes as C
    from ctypes import wintypes as W
    kernel = C.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = W.HANDLE
    process = kernel.GetCurrentProcess()
    groups = W.USHORT(64)
    group_ids = (W.USHORT * 64)()
    kernel.GetProcessGroupAffinity.argtypes = [W.HANDLE, C.POINTER(W.USHORT), C.POINTER(W.USHORT)]
    if not kernel.GetProcessGroupAffinity(process, C.byref(groups), group_ids):
        raise RuntimeError("Windows CPU topology or affinity query failed")
    if groups.value != 1:
        raise RuntimeError("CPU budgeting requires a single Windows processor group")
    mask, system = C.c_size_t(), C.c_size_t()
    kernel.GetProcessAffinityMask.argtypes = [W.HANDLE, C.POINTER(C.c_size_t), C.POINTER(C.c_size_t)]
    if not kernel.GetProcessAffinityMask(process, C.byref(mask), C.byref(system)) or not mask.value:
        raise RuntimeError("Windows CPU topology or affinity query failed")
    size = W.DWORD()
    query = kernel.GetLogicalProcessorInformationEx
    query.argtypes = [W.DWORD, C.c_void_p, C.POINTER(W.DWORD)]
    query(0, None, C.byref(size))  # RelationProcessorCore
    if not size.value:
        raise RuntimeError("Windows CPU topology or affinity query failed")
    buffer = C.create_string_buffer(size.value)
    if not query(0, buffer, C.byref(size)):
        raise RuntimeError("Windows CPU topology or affinity query failed")
    data = buffer.raw
    offset, physical = 0, 0
    width = C.sizeof(C.c_size_t)
    while offset < len(data):
        length = int.from_bytes(data[offset + 4:offset + 8], "little")
        if length < 32 or offset + length > len(data):
            raise RuntimeError("Windows CPU topology or affinity query failed")
        count = int.from_bytes(data[offset + 30:offset + 32], "little")
        for index in range(count):
            start = offset + 32 + index * (width + 8)
            if start + width + 8 > offset + length:
                raise RuntimeError("Windows CPU topology or affinity query failed")
            bits = int.from_bytes(data[start:start + width], "little")
            group = int.from_bytes(data[start + width:start + width + 2], "little")
            if group == group_ids[0] and bits & mask.value:
                physical += 1
                break
        offset += length
    logical = mask.value.bit_count()
    # Job CPU rate is a fraction of the host capacity, not the affinity mask.
    class Rate(C.Structure):
        _fields_ = [("flags", W.DWORD), ("rate", W.DWORD)]
    rate = Rate()
    query_job = kernel.QueryInformationJobObject
    query_job.argtypes = [W.HANDLE, C.c_int, C.c_void_p, W.DWORD, C.c_void_p]
    in_job = W.BOOL()
    kernel.IsProcessInJob.argtypes = [W.HANDLE, W.HANDLE, C.POINTER(W.BOOL)]
    if not kernel.IsProcessInJob(process, None, C.byref(in_job)):
        raise RuntimeError("Cannot query Windows job membership")
    if in_job.value and not query_job(None, 15, C.byref(rate), C.sizeof(rate), None):
        raise RuntimeError("Cannot query Windows CPU job limit")
    if in_job.value and rate.flags & 1:
        maximum = (rate.rate >> 16) if rate.flags & 16 else rate.rate
        if rate.flags & (4 | 16) and maximum:
            kernel.GetActiveProcessorCount.argtypes = [W.WORD]
            kernel.GetActiveProcessorCount.restype = W.DWORD
            host = kernel.GetActiveProcessorCount(0xffff) or os.cpu_count() or logical
            logical = min(logical, max(1, host * maximum // 10000))
    if not physical:
        raise RuntimeError("Windows affinity contains no physical cores")
    return logical, min(physical, logical)


def _physical_cpus() -> int:
    """Count distinct cores in the affinity mask; never guess the SMT ratio."""
    if sys.platform == "win32":
        return _windows_cpu_counts()[1]
    if sys.platform == "darwin":
        try:
            return max(1, int(subprocess.check_output(
                ["sysctl", "-n", "hw.physicalcpu"], text=True, timeout=5)))
        except (OSError, ValueError, subprocess.SubprocessError):
            return os.cpu_count() or 1
    allowed = (os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity")
               else range(os.cpu_count() or 1))
    cores = set()
    try:
        for cpu in allowed:
            root = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
            package = int((root / "physical_package_id").read_text())
            core = int((root / "core_id").read_text())
            if package < 0 or core < 0:
                return len(allowed) or 1
            cores.add((package, core))
    except (OSError, ValueError):
        return len(allowed) or 1
    return len(cores) or 1


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
    """Zero selects physical cores within allocated CPUs; reject requests outside that allocation."""
    windows = _windows_cpu_counts() if sys.platform == "win32" else None
    available = windows[0] if windows else os.cpu_count() or 1
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
        raise ValueError(f"threads must be 0 (automatic physical-core budget) or 1..{available}, got {threads!r}")
    return threads or min(available, windows[1] if windows else _physical_cpus())
