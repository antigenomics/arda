#!/usr/bin/env python3
"""Benchmark arda vs IgBLAST on synthetic IGH queries: speedup + concordance.

Generates N synthetic full-length IGH rearrangements (random human IGH scaffold +
substitutions + occasional 3-nt indel), then:
  * times arda annotation at several scales,
  * times igblastn on a capped subset,
  * reports CDR3/region concordance on the overlap.

    env ARDA_MMSEQS=$(which mmseqs) python scripts/bench_vs_igblast.py
"""

from __future__ import annotations

import os
import random
import time
from pathlib import Path

import polars as pl

from arda.paths import vdj_dir, data_dir
from arda.refbuild.imgt import read_fasta
from arda.refbuild.translate import translate
from arda import igblast
from arda.annotate.mapper import annotate_records

SIZES = [10_000, 50_000, 100_000]
IGBLAST_CAP = 10_000          # cap igblast timing/concordance set
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3")


def synth_queries(n: int, seed: int = 0) -> list[tuple[str, str]]:
    fa = [(i, s) for i, s in read_fasta(vdj_dir("human") / "alleles.fasta")
          if i.startswith("IGH_")]
    rng = random.Random(seed)
    out = []
    for k in range(n):
        _, s = fa[rng.randrange(len(fa))]
        b = list(s)
        for _ in range(rng.randint(0, 6)):           # substitutions
            p = rng.randrange(len(b))
            b[p] = rng.choice("ACGT")
        s = "".join(b)
        if rng.random() < 0.15:                       # occasional in-frame 3-nt indel
            p = rng.randrange(30, len(s) - 30)
            s = s[:p] + "AAA" + s[p:] if rng.random() < 0.5 else s[:p] + s[p + 3:]
        out.append((f"q{k}", s))
    return out


def run_arda(queries, threads):
    t0 = time.perf_counter()
    recs = annotate_records(queries, "human", "nt", threads=threads)
    return time.perf_counter() - t0, recs


def run_igblast(queries, threads):
    work = data_dir() / "bench"
    work.mkdir(parents=True, exist_ok=True)
    fa = work / "bench.fasta"
    fa.write_text("".join(f">{i}\n{s}\n" for i, s in queries))
    db = data_dir() / "blastdb" / "Homo_sapiens"
    aux = Path(igblast.bin_dir()) / "optional_file" / "human_gl.aux"
    out = work / "bench.igblast.airr.tsv"
    t0 = time.perf_counter()
    igblast.igblastn_airr(fa, out, organism="human",
                          germline_db_v=db / "IGHV", germline_db_j=db / "IGHJ",
                          germline_db_d=db / "IGHD",
                          auxiliary_data=aux if aux.exists() else None,
                          ig_seqtype="Ig", num_threads=threads)
    dt = time.perf_counter() - t0
    df = pl.read_csv(out, separator="\t", infer_schema_length=0)
    return dt, {r["sequence_id"]: r for r in df.iter_rows(named=True)}


def concordance(arda_recs, igb):
    arda = {r["sequence_id"]: r for r in arda_recs}
    comp = agree = 0
    for sid, ig in igb.items():
        ar = arda.get(sid)
        if not ar:
            continue
        for r in REGIONS:
            i = (ig.get(f"{r}_aa") or "").strip()
            a = (ar.get(f"{r}_aa") or "").strip()
            if not i or not a:
                continue
            comp += 1
            agree += (i in a or a in i) if r in ("fwr1", "cdr3") else (i == a)
    return agree, comp


def main():
    threads = os.cpu_count() or 1
    print(f"\n=== arda vs IgBLAST benchmark (threads={threads}) ===\n")
    print(f"{'N':>8} {'arda(s)':>9} {'arda seq/s':>11} {'igblast(s)':>11} "
          f"{'igb seq/s':>10} {'speedup':>8}")

    big = synth_queries(max(SIZES))
    igb_dt, igb = run_igblast(big[:IGBLAST_CAP], threads)
    igb_rate = IGBLAST_CAP / igb_dt

    for n in SIZES:
        q = big[:n]
        a_dt, recs = run_arda(q, threads)
        a_rate = n / a_dt
        speedup = a_rate / igb_rate
        igb_col = f"{igb_dt:11.2f}" if n == IGBLAST_CAP else f"{'~'+format(n/igb_rate,'.1f'):>11}"
        print(f"{n:>8} {a_dt:9.2f} {a_rate:11.0f} {igb_col} {igb_rate:10.0f} {speedup:7.1f}x")

    agree, comp = concordance(recs if len(big) == IGBLAST_CAP else run_arda(big[:IGBLAST_CAP], threads)[1], igb)
    print(f"\nregion concordance (n={IGBLAST_CAP}): {agree}/{comp} = {agree/comp:.1%}")


if __name__ == "__main__":
    main()
