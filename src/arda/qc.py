"""Batch QC: many samples' :mod:`arda.stats` tables as one cohort view.

A run already reduces gigabytes of AIRR to a few hundred short rows (``<prefix>.stats.tsv``), and
:mod:`arda.stats` says in its own docstring that the long format exists so a metric can be
``join``ed across samples. Nothing performed that join. This does.

**It reads only the QC tables.** Never an AIRR, never a clonotype table. So a 1,000-sample cohort
is ~20 MB and one ``pl.concat``, it runs on a laptop against results copied off a cluster, and it
cannot be the slow step. That is the same split as Stage 1 (per read, shards) against Stages 2-3
(global, do not) -- the expensive reduction happens once, per sample, where the data already is.

**Flags, never filters, and no shipped threshold.** ``stats`` refuses to decide whether a sample
is usable and so does this: there is no pass/fail column and no constant anywhere that says what
a good ``mapped_fraction`` is, because that number depends on the library, the organism and the
depth. What a batch *can* say is whether a sample looks like the batch it came in with, so each
metric gets its group's median, its MAD, and a robust z -- and the operator decides.

Repertoire biology is deliberately absent: no diversity, no clonality, no rarefaction, no overlap,
no cross-sample clonotype matching. Those are ``vdjtools``', which takes arda as a dependency; QC
here answers "did this run work, and is this sample like its batch".
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from ._log import logger
from .stats import STATS_COLUMNS, _typed

__all__ = ["collect_batch", "pivot", "outliers", "write_batch", "find_stats", "BATCH_COLUMNS",
           "MIN_GROUP", "Z_FLAG"]

#: The long batch table: a sample column and a label pair in front of the four `stats` ones.
BATCH_COLUMNS = ["sample", "project", "batch"] + STATS_COLUMNS

#: A group smaller than this gets a median but no z. MAD over four points is not a scale estimate,
#: and a z computed from one would flag whichever sample happened to be furthest from the middle.
MIN_GROUP = 5

#: Iglewicz-Hoaglin: |z| >= 3.5 on the MAD-based modified z is the standard outlier call. It is
#: reported, not applied -- nothing here drops or rewrites a sample.
Z_FLAG = 3.5

#: 0.6745 is the standard normal's MAD, so `0.6745 * (x - median) / MAD` puts the robust z on the
#: same scale as an ordinary one when the data happen to be normal.
_MAD_TO_SIGMA = 0.6745


def find_stats(out_dir: str | Path) -> list[Path]:
    """Every ``*.stats.tsv`` in ``out_dir``, sorted. The sample id is the file's stem."""
    return sorted(Path(out_dir).glob("*.stats.tsv"))


def _sample_id(path: Path) -> str:
    """``PT01.stats.tsv`` -> ``PT01``. The prefix a run was given IS the sample id."""
    return path.name[: -len(".stats.tsv")]


def _labels(sheet: str | Path | None) -> dict[str, tuple[str, str]]:
    if sheet is None:
        return {}
    from .samples import read_sheet

    return {s.id: (s.project, s.batch) for s in read_sheet(sheet)}


def collect_batch(stats: list[Path] | list[str], *,
                  sheet: str | Path | None = None) -> pl.DataFrame:
    """Concatenate per-sample QC tables into one long frame.

    Args:
        stats: the ``*.stats.tsv`` paths, one per sample.
        sheet: the sample sheet, read ONLY for its optional ``project``/``batch`` labels. A sample
            the sheet does not mention keeps empty labels rather than being dropped -- a cohort
            assembled from two sheets is normal, and silently losing half of it is not.

    Every value stays Utf8, exactly as ``stats`` wrote it. Casting happens in :func:`outliers`,
    per metric, because a QC table legitimately holds versions and file lists next to counts.
    """
    labels = _labels(sheet)
    frames = []
    for path in stats:
        path = Path(path)
        sid = _sample_id(path)
        df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
        if not df.height:
            logger.warning("%s: no QC rows", path.name)
            continue
        project, batch = labels.get(sid, ("", ""))
        frames.append(df.select(
            pl.lit(sid).alias("sample"), pl.lit(project).alias("project"),
            pl.lit(batch).alias("batch"),
            *[pl.col(c).cast(pl.Utf8).fill_null("") for c in STATS_COLUMNS]))
    if not frames:
        return pl.DataFrame(schema={c: pl.Utf8 for c in BATCH_COLUMNS})
    return pl.concat(frames, how="vertical")


def pivot(df: pl.DataFrame, scope: str = "sample") -> pl.DataFrame:
    """One row per sample, one column per metric, for a single scope.

    Never: the ``run`` scope is keyed by STAGE, so its metric names collide across stages
    (``map``/``correct``/``assemble`` all carry ``wall_seconds``). Pivoting it prefixes the stage,
    which is why the column is ``map.wall_seconds`` here and ``wall_seconds`` in the long table.
    """
    sub = df.filter(pl.col("scope") == scope)
    if not sub.height:
        return pl.DataFrame(schema={"sample": pl.Utf8, "project": pl.Utf8, "batch": pl.Utf8})
    name = (pl.when(pl.col("key") == "").then(pl.col("metric"))
            .otherwise(pl.col("key") + "." + pl.col("metric")))
    wide = (sub.with_columns(name.alias("_col"))
            .pivot(on="_col", index=["sample", "project", "batch"], values="value",
                   aggregate_function="first")
            .sort("sample"))
    return wide.select(["sample", "project", "batch"]
                       + sorted(c for c in wide.columns
                                if c not in ("sample", "project", "batch")))


def outliers(df: pl.DataFrame, by: tuple[str, ...] = ("project", "batch")) -> pl.DataFrame:
    """Per metric within each group: the group's median and MAD, and each sample's robust z.

    Rows whose value is not numeric are dropped -- a version string has no median. Groups smaller
    than :data:`MIN_GROUP`, and metrics whose MAD is 0 (every sample identical, which is the
    normal case for ``organism`` or ``threads``), get ``null`` rather than a z: an invented scale
    would flag the one sample that differs by a rounding error.

    Returns ``sample, *by, scope, key, metric, value, median, mad, z, outlier``.
    """
    numeric = (df.with_columns(pl.col("value").cast(pl.Float64, strict=False).alias("_v"))
               .filter(pl.col("_v").is_not_null() & pl.col("_v").is_finite()))
    if not numeric.height:
        return pl.DataFrame(schema={c: pl.Utf8 for c in BATCH_COLUMNS}
                            | {"median": pl.Float64, "mad": pl.Float64, "z": pl.Float64,
                               "outlier": pl.Int64})
    keys = list(by) + ["scope", "key", "metric"]
    stats = numeric.group_by(keys).agg(
        pl.col("_v").median().alias("median"), pl.len().alias("n"))
    joined = numeric.join(stats, on=keys, how="left")
    mad = (joined.with_columns((pl.col("_v") - pl.col("median")).abs().alias("_d"))
           .group_by(keys).agg(pl.col("_d").median().alias("mad")))
    joined = joined.join(mad, on=keys, how="left")
    z = pl.when((pl.col("n") >= MIN_GROUP) & (pl.col("mad") > 0)) \
          .then(_MAD_TO_SIGMA * (pl.col("_v") - pl.col("median")) / pl.col("mad")) \
          .otherwise(None)
    return (joined.with_columns(z.alias("z"))
            .with_columns(pl.when(pl.col("z").is_null()).then(None)
                          .otherwise((pl.col("z").abs() >= Z_FLAG).cast(pl.Int64))
                          .alias("outlier"))
            .select(["sample", *by, "scope", "key", "metric", "value",
                     "median", "mad", "z", "outlier"])
            .sort(["scope", "key", "metric", "sample"]))


def write_batch(df: pl.DataFrame, out_prefix: str | Path) -> list[Path]:
    """Write the three batch artifacts and return their paths.

    ``<prefix>.qc.tsv``
        the long table, every sample x scope x metric, with median/MAD/z where there is a scale
    ``<prefix>.qc.wide.tsv``
        one row per sample, ``sample`` scope only -- the table an operator actually reads
    ``<prefix>.qc.json``
        typed and nested, what ``arda qc report`` inlines
    """
    out_prefix = Path(out_prefix)
    scored = outliers(df)
    long_path = out_prefix.with_name(out_prefix.name + ".qc.tsv")
    wide_path = out_prefix.with_name(out_prefix.name + ".qc.wide.tsv")
    json_path = out_prefix.with_name(out_prefix.name + ".qc.json")

    scored.write_csv(long_path, separator="\t", quote_style="never", float_precision=6)
    pivot(df).write_csv(wide_path, separator="\t", quote_style="never")
    json_path.write_text(json.dumps(_document(df, scored), indent=2, sort_keys=True))
    return [long_path, wide_path, json_path]


def _document(df: pl.DataFrame, scored: pl.DataFrame) -> dict:
    """The batch as one typed JSON document: ``{samples, labels, stats, outliers}``.

    Nested by sample and scope, the way ``stats.write_stats_json`` nests one sample, so a reader
    that already parses a single run's ``.stats.json`` needs no second shape.
    """
    samples = sorted(df["sample"].unique().to_list())
    labels = {row["sample"]: {"project": row["project"], "batch": row["batch"]}
              for row in df.select(["sample", "project", "batch"]).unique()
              .iter_rows(named=True)}
    stats: dict[str, dict[str, dict[str, dict[str, object]]]] = {}
    for sample, scope, key, metric, value in df.select(
            ["sample", "scope", "key", "metric", "value"]).iter_rows():
        stats.setdefault(sample, {}).setdefault(scope, {}).setdefault(key, {})[metric] = \
            _typed(value)
    flagged = [
        {"sample": r["sample"], "scope": r["scope"], "key": r["key"], "metric": r["metric"],
         "value": r["value"], "median": r["median"], "z": r["z"]}
        for r in scored.filter(pl.col("outlier") == 1).iter_rows(named=True)
    ]
    return {"samples": samples, "labels": labels, "stats": stats, "outliers": flagged,
            "z_flag": Z_FLAG, "min_group": MIN_GROUP}
