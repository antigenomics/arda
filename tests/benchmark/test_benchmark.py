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


def _query_set(n: int) -> list[tuple[str, str]]:
    from pathlib import Path
    from arda.refbuild.imgt import read_fasta

    fa = read_fasta(Path(vdj_dir("human") / "alleles.fasta"))
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
