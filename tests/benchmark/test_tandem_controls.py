"""Issue 157 end-to-end annotation benchmark; optional raw FASTQ stays external.

RUN_BENCHMARK=1 ARDA_BENCH_FASTQ=/path/to/reads.fq.gz python -m pytest \
    tests/benchmark/test_tandem_controls.py -s
Reports wall time and process peak RSS; output comparisons exclude D/NP fields.
"""
import gzip
import os
import resource
import statistics
import sys
import time
from pathlib import Path

import pytest

from arda.annotate.mapper import annotate_records
from arda.refbuild.imgt import read_fasta

pytestmark = pytest.mark.skipif(not os.getenv('RUN_BENCHMARK'), reason='set RUN_BENCHMARK=1')


def test_tandem_modes():
    source = os.getenv('ARDA_BENCH_FASTQ')
    if source:
        opener = gzip.open if source.endswith('.gz') else open
        records = []
        with opener(source, 'rt') as fh:
            for i in range(5000):
                name, seq, plus, quality = [fh.readline().rstrip() for _ in range(4)]
                if not name:
                    break
                assert name.startswith('@') and plus.startswith('+') and len(seq) == len(quality)
                records.append((f'q{i}', seq))
    else:
        seeds = list(read_fasta(Path(__file__).resolve().parents[2] / 'examples/dd.fasta'))
        records = [(f'q{i}', seeds[i % len(seeds)][1]) for i in range(1000)]
    assert records
    baseline = None
    print('\nmode\tqueries\tmedian_seconds\tpeak_rss_mb')
    for mode, options in [('tandem', {}), ('single_d', {'map_dd': False}), ('no_d', {'map_d': False})]:
        timings = []
        for _ in range(3):
            start = time.perf_counter()
            rows = annotate_records(records, threads=0, **options)
            timings.append(time.perf_counter() - start)
            stable = [{k: v for k, v in row.items() if not k.startswith(('d_', 'd2_', 'np'))}
                      for row in rows]
            if baseline is None:
                baseline = stable
            assert stable == baseline
            if mode != 'tandem':
                assert all(not r['d2_call'] for r in rows)
            if mode == 'no_d':
                assert all(not r['d_call'] for r in rows)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1e6 if sys.platform == 'darwin' else 1e3)
        print(f'{mode}\t{len(records)}\t{statistics.median(timings):.6f}\t{rss:.3f}')
