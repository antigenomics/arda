"""Does mostly-non-AIRR input run faster? (bulk RNA-seq scenario)

mmseqs prefilters by k-mer matching: reads with no k-mer hit to any reference are
rejected before alignment. So a file that is mostly random (like blood bulk
RNA-seq, where ~1-5% of reads are actually from rearranged receptors) should be
faster per read than an all-receptor file. Measured here at fixed N.
"""
from __future__ import annotations
import os, random, time
from arda.paths import vdj_dir
from arda.refbuild.imgt import read_fasta
from arda.annotate.mapper import annotate_records

N = 100_000
RLEN = 150  # typical RNA-seq read length


def real_pool():
    return [s for i, s in read_fasta(vdj_dir("human") / "alleles.fasta") if i.startswith("IGH_")]


def make(frac_real: float, seed: int):
    rng = random.Random(seed)
    pool = real_pool()
    out = []
    for k in range(N):
        if rng.random() < frac_real:
            s = rng.choice(pool)
            # emulate an RNA-seq read: random RLEN-window of the transcript
            if len(s) > RLEN:
                p = rng.randrange(len(s) - RLEN)
                s = s[p:p + RLEN]
        else:
            s = "".join(rng.choice("ACGT") for _ in range(RLEN))
        out.append((f"q{k}", s))
    return out


def timeit(records, threads):
    t0 = time.perf_counter()
    recs = annotate_records(records, "human", "nt", threads=threads)
    dt = time.perf_counter() - t0
    hits = sum(1 for r in recs if r["v_call"])
    return dt, hits


def main():
    threads = os.cpu_count() or 1
    print(f"\n=== prefilter scenario (N={N}, read={RLEN}nt, threads={threads}) ===\n")
    print(f"{'%receptor':>10} {'time(s)':>8} {'seq/s':>9} {'annotated':>10}")
    for frac in (1.0, 0.10, 0.01):
        recs = make(frac, seed=int(frac * 100))
        dt, hits = timeit(recs, threads)
        print(f"{frac*100:9.0f}% {dt:8.2f} {N/dt:9.0f} {hits:10d}")


if __name__ == "__main__":
    main()
