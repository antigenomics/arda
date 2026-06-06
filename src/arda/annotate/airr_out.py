"""Assemble and write AIRR-formatted output."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .transfer import AIRR_COLUMNS

__all__ = ["write_airr"]


def write_airr(records: list[dict], path: str | Path) -> Path:
    """Write annotation records to an AIRR-style TSV with stable column order."""
    path = Path(path)
    if not records:
        path.write_text("\t".join(AIRR_COLUMNS) + "\n")
        return path
    df = pl.DataFrame(records).select(AIRR_COLUMNS)
    df.write_csv(path, separator="\t")
    return path
