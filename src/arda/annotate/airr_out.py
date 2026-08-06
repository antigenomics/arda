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

    This is per-record per-COLUMN, so it is one of the few Python loops whose call count scales
    with the output: profiled on a 4.58 %-receptor library it was 4.6 M generator steps and 9.2 M
    dict lookups for 54,876 rows, ~6 % of the run. Two of those lookups per column were the same
    lookup done twice (once to test for None, once to fetch), and `str()` was called on values that
    are already strings -- `transfer_hit` writes strings for every field it fills.
    """
    cols = AIRR_COLUMNS
    out = []
    for rec in records:
        vals = []
        for c in cols:
            v = rec.get(c)
            if v is None:
                vals.append("")
            elif type(v) is str:          # exact type, not isinstance: no subclass dispatch
                vals.append(v)
            else:
                vals.append(str(v))
        out.append("\t".join(vals))
    return "\n".join(out) + ("\n" if out else "")
