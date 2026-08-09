"""Assemble and write AIRR-formatted output."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from .transfer import AIRR_COLUMNS

__all__ = ["write_airr", "airr_header", "format_rows", "read_airr"]

# Hoisted out of the call: a tuple of the SAME str objects every time, so the extension's per-column
# dict lookups hit cached string hashes instead of re-interning 52 names per record.
_COLUMNS = tuple(AIRR_COLUMNS)

try:                                              # optional: built by scikit-build-core
    from .._markup import format_rows as _format_rows_cpp
except ImportError:                               # pragma: no cover - source checkout without ext
    _format_rows_cpp = None


def write_airr(records: list[dict], path: str | Path) -> Path:
    """Write annotation records to an AIRR-style TSV with stable column order."""
    path = Path(path)
    if not records:
        path.write_text("\t".join(AIRR_COLUMNS) + "\n")
        return path
    df = pl.DataFrame(records).select(AIRR_COLUMNS)
    df.write_csv(path, separator="\t", quote_style="never")
    return path


def airr_header(extra_columns: tuple[str, ...] = ()) -> str:
    """The AIRR TSV header line (no trailing newline).

    ``extra_columns`` appends non-schema fields (today only ``junction_quality``, emitted by
    ``arda rnaseq map --junction-quality``). They go at the END, after every shipped column, so a
    consumer reading the shipped set by position is unaffected; with the default ``()`` the header
    is byte-identical to what it has always been.
    """
    return "\t".join(list(AIRR_COLUMNS) + list(extra_columns))


def format_rows(records: list[dict], extra_columns: tuple[str, ...] = ()) -> str:
    """Format records as TSV rows (column order, trailing newline per row).

    ``extra_columns`` must match what :func:`airr_header` was given. A record missing that key
    renders as an empty field in both implementations (pinned by ``tests/unit/test_airr_out.py``).

    Used by the streaming writer in ``mapper.annotate_file`` to append chunks
    incrementally without holding the whole output in memory.

    This is per-record per-COLUMN, so it is one of the few Python loops whose call count scales
    with the output: profiled on a 4.58 %-receptor library it was 4.6 M generator steps and 9.2 M
    dict lookups for 54,876 rows, ~6 % of the run. Two of those lookups per column were the same
    lookup done twice (once to test for None, once to fetch), and `str()` was called on values that
    are already strings -- `transfer_hit` writes strings for every field it fills.

    Re-profiled on the amplicon path once the search term collapsed, it is **12.7 % of the wall**
    (0.697 s for 54,178 records x 52 columns), which is why the loop below now lives in
    ``_markup.format_rows``. The Python version is kept as the reference implementation and the
    fallback when the extension is not built; the two are asserted byte-identical by
    ``tests/unit/test_airr_out.py``.
    """
    cols = _COLUMNS + tuple(extra_columns) if extra_columns else _COLUMNS
    if _format_rows_cpp is not None:
        return _format_rows_cpp(records, cols)
    return _format_rows_py(records, cols)


def _format_rows_py(records: list[dict], cols=_COLUMNS) -> str:
    """Reference implementation of :func:`format_rows` (see its docstring)."""
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

def read_airr(path):
    """Read an AIRR TSV, tolerating BOTH dialects arda has written.

    ⛔ ``quote_char=None`` is not a style choice. ``junction_quality`` is a Phred+33 string and
    **chr 34 is ``"``, i.e. Q1** -- a legitimate score any low-quality base produces. polars' reader
    treats it as a quote character, so ONE such base collapses the parse of the whole file
    (``CSV malformed: expected 1 rows, actual 155 rows``). Measured on a real Raji run: exactly one
    row contained a ``"`` and the entire table became unreadable.

    The complication is that arda has written the format two ways. The streaming writer
    (``_markup.format_rows``, which produces every ``map`` output) emits **raw** fields and a truly
    empty string for a missing value. polars' ``write_csv`` quotes -- it renders an empty string as
    the two characters ``""`` and doubles an embedded quote. So reading unquoted is right for the
    big files and would turn every empty field of an older polars-written one into a literal ``""``.
    Hence the normalisation below: an AIRR field is never legitimately the two-character string
    ``""``, so mapping it back to empty is unambiguous and makes both dialects read the same.
    """
    df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
    return df.with_columns(
        pl.when(pl.col(c) == '""').then(pl.lit("")).otherwise(pl.col(c)).alias(c)
        for c, t in df.schema.items() if t == pl.Utf8)
