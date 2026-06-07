"""Benchmarks: runtime, peak memory, and multiprocessing/thread scalability.

Guarded by ``RUN_BENCHMARK=1`` so it never runs in the normal suite. Builds a
scaled query set from the human reference scaffolds and annotates it at several
thread counts, reporting a structured summary table.

    env RUN_BENCHMARK=1 ARDA_MMSEQS=$(which mmseqs) \\
        python -m pytest tests/benchmark -s
"""

from __future__ import annotations

import os
import time
import tracemalloc
import random

import pytest

from arda.paths import vdj_dir
from arda.annotate.mapper import annotate_records
from tests.conftest import requires_mmseqs, requires_human_db

pytestmark = [
    pytest.mark.skipif(not os.getenv("RUN_BENCHMARK"), reason="set RUN_BENCHMARK=1"),
    requires_mmseqs,
    requires_human_db,
]


def _query_set(n: int, loci: tuple[str, ...] | None = None) -> list[tuple[str, str]]:
    from pathlib import Path
    from arda.refbuild.imgt import read_fasta

    fa = read_fasta(Path(vdj_dir("human") / "alleles.fasta"))
    if loci:
        fa = [(sid, seq) for sid, seq in fa if sid.split("_")[0] in loci]
    rng = random.Random(0)
    out = []
    for k in range(n):
        sid, seq = fa[rng.randrange(len(fa))]
        out.append((f"q{k}", seq))
    return out


@pytest.mark.parametrize("n", [1000, 5000])
def test_throughput_and_memory(n):
    queries = _query_set(n)
    tracemalloc.start()
    t0 = time.perf_counter()
    recs = annotate_records(queries, "human", "nt", threads=os.cpu_count() or 1)
    dt = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert len(recs) == n
    print(f"\n[bench] n={n:>6}  time={dt:6.2f}s  "
          f"rate={n/dt:8.1f} seq/s  peak_py={peak/1e6:7.1f} MB")


def test_thread_scaling():
    queries = _query_set(5000)
    print("\n[bench] thread scaling (n=5000):")
    base = None
    for threads in (1, 2, 4, max(1, os.cpu_count() or 1)):
        t0 = time.perf_counter()
        annotate_records(queries, "human", "nt", threads=threads)
        dt = time.perf_counter() - t0
        base = base or dt
        print(f"  threads={threads:>2}  time={dt:6.2f}s  speedup={base/dt:4.2f}x")


def _timed(queries, *, map_d, threads, repeats=3):
    """Best-of-`repeats` wall time for annotating `queries` (warms caches first)."""
    annotate_records(queries, "human", "nt", threads=threads, map_d=map_d)  # warm-up
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        annotate_records(queries, "human", "nt", threads=threads, map_d=map_d)
        best = min(best, time.perf_counter() - t0)
    return best


@pytest.mark.parametrize("n", [2000, 10000])
def test_d_mapping_overhead(n):
    """Extra wall time from D mapping, measured on a VDJ-only (IGH/TRB/TRD) query
    set — the worst case, since D mapping only runs for VDJ-locus hits. Reports
    absolute and per-sequence overhead so the cost of the option is on record."""
    threads = os.cpu_count() or 1
    queries = _query_set(n, loci=("IGH", "TRB", "TRD"))
    off = _timed(queries, map_d=False, threads=threads)
    on = _timed(queries, map_d=True, threads=threads)
    overhead = on - off
    pct = 100.0 * overhead / off if off else 0.0
    print(f"\n[bench] D-mapping overhead (VDJ-only, n={n:>6}, threads={threads}): "
          f"off={off:6.3f}s  on={on:6.3f}s  "
          f"extra={overhead*1e3:7.1f} ms ({pct:+5.1f}%)  "
          f"per_seq={overhead/n*1e6:6.2f} us")
    # Sanity: D mapping is a short per-hit local alignment, not a new mmseqs pass.
    assert overhead < off  # never more than doubles end-to-end time
