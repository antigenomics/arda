#!/usr/bin/env python
# Fit arda's RNA-seq wall-time cost model from completed runs, or measure it directly.
#
# The model, fitted on 10 full-depth cluster runs (754.7M reads, 41,461 s, 16 threads):
#
#     wall_map ~= total_reads / A  +  mapped_reads / B        A ~= 44,470/s   B ~= 681/s
#
# Two terms, and telling them apart is the whole point:
#
#   * the SCAN term (reads/A) is what mmseqs pays to reject a non-receptor read;
#   * the ALIGN term (mapped/B) is what it pays to align one that hits.
#
# A mapped read costs ~65x a non-mapped one. So the split between the two terms is set by the
# library's receptor fraction, not by anything arda chooses -- and that decides which
# optimisation can possibly help:
#
#   receptor fraction   scan share   align share   what a k-mer prefilter can reach
#   0.0003% (HepG2)        ~100%          ~0%      almost all of it
#   0.18%  (GM12878)         89%          11%      most of it
#   2.7%   (bulk, rich)      35%          65%      a third of it
#   48%    (amplicon)         9%          91%      almost none of it
#
# A prefilter only ever attacks the scan term. On a high-receptor library the align term
# dominates and no amount of prefiltering touches it.
#
# Usage:
#   # fit from runs that already exist (any dir of <sample>/<sample>.arda.json + arda.time)
#   python scripts/bench_cost_model.py --fit results/round4/results
#
#   # measure here: one input, several depths
#   ARDA_MMSEQS=$(which mmseqs) python scripts/bench_cost_model.py \
#       --r1 R1.fq.gz --r2 R2.fq.gz --depths 50000,200000,800000 --threads 8
#
# 2026-08-05
"""Fit or measure arda's RNA-seq map cost model."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def _fit(rows: list[tuple[int, int, float]]) -> tuple[float, float, list[float]]:
    """Least-squares fit of ``wall = reads/A + mapped/B``; returns (A, B, residual %)."""
    try:
        import numpy as np
    except ImportError:  # keep the script usable without numpy
        raise SystemExit("--fit needs numpy (pip install numpy)")
    M = np.array([[r, m] for r, m, _ in rows], dtype=float)
    y = np.array([w for _, _, w in rows], dtype=float)
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    pred = M @ coef
    resid = [100 * (p - o) / o for p, o in zip(pred, y)]
    return 1 / coef[0], 1 / coef[1], resid


def fit_from_runs(root: Path) -> None:
    rows, names = [], []
    for jf in sorted(root.glob("*/*.arda.json")):
        rep = json.loads(jf.read_text())
        m = rep.get("map") or {}
        if not m.get("total_reads") or not m.get("wall_seconds"):
            continue
        rows.append((m["total_reads"], m.get("mapped_reads", 0), m["wall_seconds"]))
        names.append(jf.parent.name)
    if len(rows) < 3:
        raise SystemExit(f"need >=3 runs with a map report under {root}, found {len(rows)}")

    A, B, resid = _fit(rows)
    print(f"wall_map = reads/{A:,.0f} + mapped_reads/{B:,.0f}      "
          f"(a mapped read costs {A / B:.0f}x a non-mapped one)\n")
    print(f"{'sample':<14}{'reads':>14}{'mapped':>12}{'recep%':>8}"
          f"{'wall_s':>10}{'resid%':>8}{'scan%':>7}{'align%':>7}")
    for (r, mp, w), n, res in zip(rows, names, resid):
        scan, align = r / A, mp / B
        print(f"{n:<14}{r:>14,}{mp:>12,}{100 * mp / r:>7.3f}%{w:>10.1f}{res:>+8.1f}"
              f"{100 * scan / (scan + align):>7.0f}{100 * align / (scan + align):>7.0f}")


def measure(r1: Path, r2: Path | None, depths: list[int], threads: int, chunk: int) -> None:
    rows = []
    print(f"{'depth':>10}{'reads':>12}{'mapped':>10}{'wall_s':>10}{'reads/s':>10}{'rss_mb':>9}")
    with tempfile.TemporaryDirectory(prefix="arda_cost_") as td:
        for n in depths:
            out, rep = Path(td) / f"{n}.tsv", Path(td) / f"{n}.json"
            cmd = ["arda", "rnaseq", "map", "--r1", str(r1), "-o", str(out),
                   "--report", str(rep), "--threads", str(threads),
                   "--chunk-size", str(chunk), "-n", str(n)]
            if r2:
                cmd += ["--r2", str(r2)]
            t0 = time.perf_counter()
            proc = subprocess.run(cmd, capture_output=True, text=True)
            wall = time.perf_counter() - t0
            if proc.returncode != 0:
                print(f"{n:>10}  FAILED: {proc.stderr.strip().splitlines()[-1:]}")
                continue
            m = json.loads(rep.read_text())
            rows.append((m["total_reads"], m["mapped_reads"], m["wall_seconds"]))
            print(f"{n:>10}{m['total_reads']:>12,}{m['mapped_reads']:>10,}"
                  f"{m['wall_seconds']:>10.1f}{m['reads_per_second']:>10,.0f}"
                  f"{m['peak_rss_mb']:>9.0f}   (external {wall:.1f}s)")
    if len(rows) >= 3:
        A, B, resid = _fit(rows)
        print(f"\nwall_map = reads/{A:,.0f} + mapped_reads/{B:,.0f}   "
              f"(max |residual| {max(abs(x) for x in resid):.1f}%)")
    else:
        print("\n(need >=3 successful depths to fit)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fit", type=Path, help="Directory of completed runs to fit from.")
    ap.add_argument("--r1", type=Path)
    ap.add_argument("--r2", type=Path)
    ap.add_argument("--depths", default="50000,200000,800000")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--chunk-size", type=int, default=200_000)
    args = ap.parse_args()

    if args.fit:
        fit_from_runs(args.fit)
    elif args.r1:
        measure(args.r1, args.r2, [int(x) for x in args.depths.split(",")],
                args.threads, args.chunk_size)
    else:
        ap.error("pass --fit <dir> or --r1 <fastq>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
