#!/usr/bin/env python
"""Audit `arda.cdr3fix` over VDJdb's curation corpus -- issue #141.

2026-09-30. Reads every VDJdb chunk TSV, reduces it to distinct `(species, cdr3, v, j)` keys,
marks each up, and writes one parquet of decisions plus a summary to stdout.

The point is a row-level A/B: run it on two checkouts and diff the parquets. Comparing counts is
not enough -- a count-equal, decision-different leg is invisible in a summary table (CLAUDE.md,
"compare a call digest, not a count"), so the parquet carries the key and every decision column
and the summary prints a digest over all of them.

    python scripts/audit_cdr3fix.py --chunks ~/vcs/code/vdjdb-db --out /tmp/audit_before.parquet

Nothing is cached: the key set is rebuilt from the chunk TSVs on every run.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import sys
import time
from pathlib import Path

import polars as pl

from arda.cdr3fix import markup_batch

#: Every directory of curation chunks a VDJdb checkout carries. Each row holds an alpha and a beta
#: chain in separate column families, so one row yields up to two keys.
CHUNK_DIRS = ("chunks", "chunks_unformatted", "chunks_with_unconventional_aa", "chunks_negative")

#: The decision columns an A/B compares. `errors` is deliberately included: a change that leaves
#: every verdict alone but alters what it REPORTS is still a change the reader has to see.
DECISION = ["cdr3_repaired", "v_call", "j_call", "v_end", "j_start",
            "v_end_nt", "j_start_nt", "v_fix", "j_fix", "good", "n_errors", "errors"]


def read_keys(root: Path) -> pl.DataFrame:
    """Distinct `(species, cdr3, v, j)` over every chunk TSV under `root`."""
    frames = []
    for d in CHUNK_DIRS:
        for p in sorted(glob.glob(str(root / d / "*.txt"))):
            try:
                df = pl.read_csv(p, separator="\t", infer_schema_length=0, quote_char=None,
                                 truncate_ragged_lines=True, ignore_errors=True)
            except Exception as exc:                      # a chunk mid-curation may not parse
                print(f"  skip {Path(p).name}: {exc}", file=sys.stderr)
                continue
            for chain in ("alpha", "beta"):
                need = [f"cdr3.{chain}", f"v.{chain}", f"j.{chain}", "species"]
                if all(c in df.columns for c in need):
                    frames.append(df.select(
                        pl.col(f"cdr3.{chain}").alias("cdr3"),
                        pl.col(f"v.{chain}").alias("v"),
                        pl.col(f"j.{chain}").alias("j"),
                        pl.col("species")))
    if not frames:
        raise SystemExit(f"no chunk TSVs under {root}")
    return (pl.concat(frames, how="vertical_relaxed")
            .with_columns(pl.all().str.strip_chars().fill_null(""))
            .filter((pl.col("cdr3") != "") & (pl.col("v") != "") & (pl.col("j") != ""))
            .unique(subset=["species", "cdr3", "v", "j"])
            .sort(["species", "cdr3", "v", "j"]))


def digest(df: pl.DataFrame) -> str:
    """One hash over the key and every decision column, in sorted key order."""
    cols = ["species", "cdr3", "v", "j"] + [c for c in DECISION if c in df.columns]
    h = hashlib.sha256()
    for row in df.select(cols).sort(["species", "cdr3", "v", "j"]).iter_rows():
        h.update(("\t".join("" if x is None else str(x) for x in row) + "\n").encode())
    return h.hexdigest()[:16]


def summarise(df: pl.DataFrame) -> None:
    n = df.height
    print(f"\nkeys                : {n}")
    print(f"good                : {df['good'].sum()} ({df['good'].mean():.4%})")
    print(f"fix_needed          : {df['fix_needed'].sum()}")
    print(f"carries an error    : {(df['n_errors'] > 0).sum()}")
    # The defect-4 metric is `good` alongside an error that was NOT applied. `good` alongside an
    # APPLIED error is the normal, correct case: the record was repaired and the edit is its record.
    unrepaired = pl.col("errors").str.contains("reported, not repaired")
    print(f"good + unrepaired   : {df.filter(pl.col('good') & unrepaired).height}"
          f"    <- issue #141 defect 4, must be 0")
    print(f"good + a repair     : {(df['good'] & (df['n_errors'] > 0)).sum()}")
    print(f"decision digest     : {digest(df)}")

    for side in ("v", "j"):
        if f"{side}_flags" in df.columns:
            print(f"\n{side}_flags:")
            counts = (df.select(pl.col(f"{side}_flags").str.split(",").alias("f"))
                      .explode("f").group_by("f").len().sort("len", descending=True))
            for flag, c in counts.iter_rows():
                print(f"  {flag or '-':<14} {c:>8}")

    print("\nlegacy fix types (v_fix x j_fix), top 15:")
    top = (df.group_by(["v_fix", "j_fix"]).len().sort("len", descending=True).head(15))
    for vf, jf, c in top.iter_rows():
        print(f"  {vf:<20} {jf:<20} {c:>8}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunks", type=Path, default=Path.home() / "vcs/code/vdjdb-db",
                    help="a vdjdb-db checkout (reads its chunk*/ directories)")
    ap.add_argument("--out", type=Path, required=True, help="parquet of per-key decisions")
    ap.add_argument("--limit", type=int, default=0, help="first N keys only (smoke runs)")
    args = ap.parse_args()

    t0 = time.perf_counter()
    keys = read_keys(args.chunks)
    if args.limit:
        keys = keys.head(args.limit)
    t_read = time.perf_counter() - t0
    print(f"corpus              : {keys.height} keys from {args.chunks} in {t_read:.1f}s")

    t0 = time.perf_counter()
    out = markup_batch(keys, cdr3="cdr3", v="v", j="j", species="species")
    t_mark = time.perf_counter() - t0
    # `markup_batch` already carries `cdr3` and `species`; keep the SUBMITTED calls beside the
    # resolved ones, because "the allele moved" is itself one of the verdicts being audited.
    out = out.with_columns(v=keys["v"], j=keys["j"])
    print(f"markup              : {t_mark:.1f}s "
          f"({keys.height / max(t_mark, 1e-9):,.0f} keys/s)")
    out.write_parquet(args.out)
    print(f"wrote               : {args.out}")
    summarise(out)


if __name__ == "__main__":
    main()
