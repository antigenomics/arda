"""Assemble and write AIRR-formatted output."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .transfer import AIRR_COLUMNS

__all__ = ["write_airr", "airr_header", "format_rows"]


def write_airr(records: list[dict], path: str | Path) -> Path:
    """Write annotation records to an AIRR-style TSV with stable column order."""
    path = Path(path)
    if not records:
        path.write_text("\t".join(AIRR_COLUMNS) + "\n")
        return path
    df = pl.DataFrame(records).select(AIRR_COLUMNS)
    df.write_csv(path, separator="\t")
    return path


def airr_header() -> str:
    """The AIRR TSV header line (no trailing newline)."""
    return "\t".join(AIRR_COLUMNS)


def format_rows(records: list[dict]) -> str:
    """Format records as TSV rows (column order, trailing newline per row).

    Used by the streaming writer in ``mapper.annotate_file`` to append chunks
    incrementally without holding the whole output in memory.
    """
    out = []
    for rec in records:
        out.append("\t".join(
            "" if rec.get(c) is None else str(rec.get(c, "")) for c in AIRR_COLUMNS))
    return "\n".join(out) + ("\n" if out else "")
