#!/usr/bin/env python
"""A/B `arda.cdr3fix` against a shipped VDJdb release -- issues #141, #711.

2026-09-30. The `cdr3fix` column of any VDJdb release is the retired `Cdr3Fixer`'s own answer:
`cdr3_old` (what the curator submitted), `cdr3` (what it repaired that to), `vId`/`jId` (the
segments it used), `vEnd`/`jStart` (the boundaries it found). Feeding arda the SAME four inputs
makes the comparison input-for-input with no join and no coverage loss -- which is the point, since
joining a release against a checkout's chunks loses 2.55 % of keys to curation drift and those keys
are not a random sample.

    python scripts/compare_vdjdb_release.py \
        --release ~/vcs/code/vdjdb-db/database/vdjdb.txt --out /tmp/ab.parquet

`--release` also takes the release zip. Nothing is cached: the key set and every arda decision are
recomputed from the release table on each run (CLAUDE.md §0d).
"""
from __future__ import annotations

import argparse
import io
import json
import time
import zipfile
from pathlib import Path

import polars as pl

from arda.cdr3fix import markup_batch

#: Legacy fix types that mean "the junction was rewritten".
_REPAIRED = ("FixAdd", "FixTrim", "FixReplace", "Realign")


def read_release(path: Path) -> pl.DataFrame:
    """Distinct `(species, cdr3_old, vId, jId)` keys with the release's own verdict."""
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            name = next(n for n in zf.namelist() if n.endswith("vdjdb.txt"))
            raw = zf.read(name)
    else:
        raw = path.read_bytes()
    df = pl.read_csv(io.BytesIO(raw), separator="\t", infer_schema_length=0, quote_char=None,
                     columns=["cdr3", "v.segm", "j.segm", "species", "cdr3fix"])
    fix = [json.loads(x) for x in df["cdr3fix"]]
    out = pl.DataFrame({
        "species": df["species"],
        "cdr3_old": [f.get("cdr3_old", "") for f in fix],
        "v": [f.get("vId", "") for f in fix],
        "j": [f.get("jId", "") for f in fix],
        "cdr3_rel": [f.get("cdr3", "") for f in fix],
        "v_fix_rel": [f.get("vFixType", "") for f in fix],
        "j_fix_rel": [f.get("jFixType", "") for f in fix],
        "v_end_rel": [int(f.get("vEnd", -1)) for f in fix],
        "j_start_rel": [int(f.get("jStart", -1)) for f in fix],
        "good_rel": [bool(f.get("good", False)) for f in fix],
    })
    return (out.filter((pl.col("cdr3_old") != "") & (pl.col("v") != "") & (pl.col("j") != ""))
            .unique(subset=["species", "cdr3_old", "v", "j"])
            .sort(["species", "cdr3_old", "v", "j"]))


def compare(rel: pl.DataFrame) -> pl.DataFrame:
    got = markup_batch(rel.rename({"cdr3_old": "cdr3"}), cdr3="cdr3", v="v", j="j",
                       species="species")
    return pl.concat([rel, got.select(pl.exclude("species", "cdr3"))], how="horizontal")


def report(df: pl.DataFrame) -> None:
    n = df.height
    rel_fixed = pl.col("v_fix_rel").is_in(_REPAIRED) | pl.col("j_fix_rel").is_in(_REPAIRED)
    agree = pl.col("cdr3_repaired") == pl.col("cdr3_rel")
    arda_fixed = pl.col("cdr3_repaired") != pl.col("cdr3_old")

    print(f"\nrelease keys (distinct species, cdr3_old, vId, jId) : {n:,}")
    print(f"repaired junction agrees with the release           : "
          f"{df.filter(agree).height:,} ({df.select(agree.mean()).item():.4%})")
    for side, col in (("vEnd", "v_end"), ("jStart", "j_start")):
        sub = df.filter(pl.col(f"{col}_rel") >= 0)
        ok = sub.filter(pl.col(col) == pl.col(f"{col}_rel")).height
        print(f"{side} agrees                                       : "
              f"{ok:,}/{sub.height:,} ({ok / sub.height:.4%})")
    rf = df.filter(rel_fixed)
    print(f"the release's own repairs reproduced                : "
          f"{rf.filter(agree).height:,} of {rf.height:,}")
    over = df.filter(~rel_fixed & arda_fixed)
    print(f"repairs the release did NOT make                    : {over.height:,}")
    print(f"  of those, into an already-canonical junction      : "
          f"{over.filter(pl.col('cdr3_old').str.starts_with('C')).height:,}")
    print(f"good beside an unrepaired disagreement              : "
          f"{df.filter(pl.col('good') & pl.col('errors').str.contains('not repaired')).height:,}")
    print(f"arda good / release good                            : "
          f"{df['good'].sum():,} / {df['good_rel'].sum():,}")

    print("\nover-fixes by flag pair, top 12:")
    top = (over.group_by(["v_flags", "j_flags"]).len().sort("len", descending=True).head(12))
    for vf, jf, c in top.iter_rows():
        print(f"  V[{vf or '-'}] J[{jf or '-'}]  {c:>7}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", type=Path, required=True,
                    help="a shipped vdjdb.txt (or the release zip)")
    ap.add_argument("--out", type=Path, required=True, help="parquet of the row-level A/B")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    rel = read_release(args.release)
    if args.limit:
        rel = rel.head(args.limit)
    print(f"release             : {args.release} -> {rel.height:,} distinct keys")
    t0 = time.perf_counter()
    df = compare(rel)
    dt = time.perf_counter() - t0
    print(f"markup              : {dt:.1f}s ({rel.height / max(dt, 1e-9):,.0f} keys/s)")
    df.write_parquet(args.out)
    print(f"wrote               : {args.out}")
    report(df)


if __name__ == "__main__":
    main()
