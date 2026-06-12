# arda — achieved performance & accuracy

Measured on Apple M3 (16 threads), MMseqs2 18.x, IgBLAST 1.22.0, against the
committed reference DB and test fixtures. All numbers are reproducible with the
scripts/tests noted; this file is the on-record snapshot.

## Accuracy vs IgBLAST (gold standard)

Offline, on ~7.3k real GenBank mRNA across all five organisms
(`tests/assets/realworld/`, compared on IgBLAST-productive records).

| organism | region concordance (productive) |
|---|---|
| human | 98.6% |
| mouse | 99.6% |
| rat | 98.0% |
| rabbit | 99.7% |
| rhesus_monkey | 99.4% |

- **V gene** assignment agrees ~**100%**.
- **junction_aa / cdr3_aa** match IgBLAST ~**99–100%** on productive-canonical
  records, and satisfy the AIRR invariant `cdr3_aa == junction_aa[1:-1]` exactly
  for *every* emitted junction (including out-of-frame ones, rendered with a `_`).
- **D gene** (where both tools call a D): TRB/TRD ~**97%** gene agreement; IGH
  ~**46–69%** — IGH D is inherently ambiguous (≈50 paralogous germlines + SHM),
  consistent with inter-tool reports.
- Remaining region diffs are one-residue FR1/CDR3 boundary conventions, not errors.

Reproduce: `pytest tests/realworld -s` (offline); raw vs-IgBLAST run on synthetic
IGH: `scripts/bench_vs_igblast.py` (98.9% region concordance at n=10k).

## Speed vs IgBLAST

Synthetic human IGH, 16 threads (`scripts/bench_vs_igblast.py`):

| sequences | arda | arda rate | speedup vs IgBLAST |
|---:|---:|---:|---:|
| 10,000 | 5.5 s | ~1.8k/s | **4.4×** |
| 50,000 | 16 s | ~3.0k/s | **7.3×** |
| 100,000 | 30 s | ~3.3k/s | **7.9×** |

IgBLAST runs at ~0.44k seq/s on the same data. arda's rate rises with input size
as fixed costs (DB load) amortize; for small batches that overhead dominates.

## Bulk RNA-seq (prefilter)

mmseqs k-mer prefilter rejects non-receptor reads before alignment, so realistic
bulk RNA-seq (~1–5% receptor) is several-fold faster than amplicon. 150 nt reads,
16 threads (`scripts/bench_prefilter.py`):

| receptor content | throughput |
|---:|---:|
| 100% (amplicon) | ~5.7k reads/s |
| 10% | ~19k reads/s |
| 1% (blood RNA-seq) | ~25k reads/s |

**32-core / 30M-read estimate:** ~10–20 min for a ~1%-receptor library — same order
of magnitude as a STAR genome-alignment pass (STAR is faster per read but indexes
the whole genome). Memory is flat (bounded-chunk streaming); multi-node via
`arda slurm`.

## D-mapping overhead

D mapping is a short per-hit gapless local alignment (C++ `d_local_align`), not a
second mmseqs pass. On a VDJ-only query set (worst case), 16 threads
(`tests/benchmark/test_d_mapping_overhead`): **~±1%** end-to-end, **~7 µs/seq**.
On by default; `--no-map-d` to disable.

## How to reproduce

```bash
# accuracy (offline, committed fixtures)
ARDA_MMSEQS=$(which mmseqs) pytest tests/realworld -s

# guarded benchmarks (timing / memory / scaling / D-overhead)
RUN_BENCHMARK=1 ARDA_MMSEQS=$(which mmseqs) pytest tests/benchmark -s

# headline speed + concordance vs IgBLAST
ARDA_MMSEQS=$(which mmseqs) python scripts/bench_vs_igblast.py
ARDA_MMSEQS=$(which mmseqs) python scripts/bench_prefilter.py
```
