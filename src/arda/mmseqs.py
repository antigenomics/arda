"""Thin wrapper around the ``mmseqs`` binary.

Inspired by pymmseqs (MIT) but deliberately dependency-free: we only need
binary discovery, a subprocess runner, and the ``createdb`` / ``search`` /
``convertalis`` (and ``easy-search``) pipeline used by the annotator.

Discovery order for the binary: ``$ARDA_MMSEQS`` → ``<project>/bin/mmseqs`` →
``mmseqs`` on ``PATH``. If none are found, a static binary is auto-fetched into
``<project>/bin/mmseqs`` (one-time, transparent) unless ``$ARDA_NO_AUTO_FETCH``
is set — so neither pip nor conda users need to install mmseqs manually.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from functools import lru_cache

from .paths import bin_dir

__all__ = [
    "MMseqsError",
    "mmseqs_binary",
    "run",
    "version",
    "version_key",
    "createdb",
    "search",
    "convertalis",
    "top_hit",
    "easy_search",
    "DEFAULT_FORMAT_OUTPUT",
]

# Fields needed to transfer reference markup onto a query. 1-based inclusive
# coords for qstart/qend/tstart/tend; cigar + qaln/taln drive the projection.
DEFAULT_FORMAT_OUTPUT = (
    "query,target,qstart,qend,tstart,tend,qlen,tlen,"
    "alnlen,mismatch,gapopen,cigar,qaln,taln,evalue,bits,pident"
)

# search-type values (see `mmseqs search --help`).
SEARCH_TYPE_AUTO = 0
SEARCH_TYPE_PROTEIN = 1   # aa query vs aa target
SEARCH_TYPE_TRANSLATED = 2  # nt query vs aa target (blastx-like)
SEARCH_TYPE_NUCLEOTIDE = 3  # nt query vs nt target


class MMseqsError(RuntimeError):
    """Raised when an ``mmseqs`` invocation exits non-zero."""


@lru_cache(maxsize=1)
def mmseqs_binary() -> str:
    """Locate the ``mmseqs`` executable, auto-fetching a static build if needed.

    Resolution: ``$ARDA_MMSEQS`` → ``<project>/bin/mmseqs`` → ``mmseqs`` on
    ``PATH``. If still not found, download a static binary into
    ``<project>/bin/mmseqs`` (one-time) unless ``$ARDA_NO_AUTO_FETCH`` is set.
    """
    env = os.environ.get("ARDA_MMSEQS")
    if env:
        return env
    local = bin_dir() / "mmseqs"
    if local.exists():
        return str(local)
    found = shutil.which("mmseqs")
    if found:
        return found
    if "ARDA_NO_AUTO_FETCH" not in os.environ:
        fetched = _auto_fetch()
        if fetched is not None:
            return fetched
    raise MMseqsError(
        "mmseqs binary not found. Install it (conda install -c bioconda mmseqs2), "
        "set $ARDA_MMSEQS, or allow auto-fetch (unset $ARDA_NO_AUTO_FETCH)."
    )


def _auto_fetch() -> str | None:
    """Download a static mmseqs binary into ``bin/``; return its path or None."""
    try:
        from ._mmseqs_fetch import fetch

        return str(fetch(bin_dir()))
    except Exception as exc:  # noqa: BLE001 - network/layout failures are non-fatal here
        import warnings
        warnings.warn(f"mmseqs auto-fetch failed: {exc}", RuntimeWarning, stacklevel=2)
        return None


def run(args: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Run ``mmseqs <args>`` capturing stdout/stderr."""
    cmd = [mmseqs_binary(), *map(str, args)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise MMseqsError(
            f"`{' '.join(cmd)}` failed (exit {proc.returncode}):\n{proc.stderr}"
        )
    return proc


@lru_cache(maxsize=1)
def version() -> str:
    """Return the mmseqs version string (e.g. ``18-8cc5c``)."""
    return run(["version"]).stdout.strip()


def version_key(v: str) -> str:
    """Canonical form of a version string, for comparing an index marker to the running mmseqs.

    ``mmseqs version`` prints the same release+commit with different punctuation across builds:
    the official static binary says ``18-8cc5c``, the bioconda build ``18.8cc5c``. They are the
    same mmseqs and produce byte-compatible indexes; only the separator differs. Comparing the raw
    strings with ``==`` therefore rejected every committed index the moment the toolchain moved from
    conda's mmseqs to the static one -- so the precompiled DBs shipped in ``database/`` were never
    used and every run rebuilt a private cache instead.

    Fold each run of separator characters to a single ``-`` and lowercase. This bridges the cosmetic
    difference while still distinguishing genuinely different versions (``17-b804f`` != ``18-8cc5c``),
    so an incompatible index is never accepted.
    """
    return re.sub(r"[\s._-]+", "-", v.strip().lower())


def createdb(fasta: str | Path, db: str | Path, *, dbtype: int | None = None) -> Path:
    """Create an mmseqs sequence DB from a FASTA file.

    ``dbtype``: ``None`` auto-detect, ``1`` amino-acid, ``2`` nucleotide.
    """
    args = ["createdb", str(fasta), str(db)]
    if dbtype is not None:
        args += ["--dbtype", str(dbtype)]
    run(args)
    return Path(db)


def search(
    query_db: str | Path,
    target_db: str | Path,
    result_db: str | Path,
    tmp_dir: str | Path,
    *,
    search_type: int = SEARCH_TYPE_AUTO,
    sensitivity: float = 5.7,
    evalue: float = 1e-3,
    max_seqs: int = 300,
    threads: int = 1,
    kmer: int | None = None,
    extra: list[str] | None = None,
) -> Path:
    """Run ``mmseqs search`` with backtrace enabled (``-a``).

    Args:
        kmer: MMseqs2 ``-k``. **This is the memory knob.** The nucleotide prefilter allocates a
            k-mer index table of 4**k entries, so the default k=15 costs 4**15 * 8 B ~ 8.6 GB
            regardless of database size, thread count or chunk size. k=13 costs ~0.7 GB. Leave
            ``None`` for MMseqs2's own default.
    """
    args = [
        "search", str(query_db), str(target_db), str(result_db), str(tmp_dir),
        "--search-type", str(search_type),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "--threads", str(threads),
        "-a",  # keep backtrace so convertalis can emit cigar/qaln/taln
    ]
    if kmer is not None:
        args += ["-k", str(kmer)]
    if extra:
        args += extra
    run(args)
    return Path(result_db)


def top_hit(result_db: str | Path, out_db: str | Path) -> Path:
    """Reduce an alignment DB to the single best-scoring hit per query.

    MMseqs2 already stores each query's results sorted by descending score, so taking the first line
    per entry *is* the best hit -- this is the idiom mmseqs' own `filterdb` usage message shows.
    Verified on 100 k reads against the human reference: 4,101 queries, identical target and identical
    bit score to a full sort-and-dedupe in polars, on every one.

    Why it matters: with `--max-seqs 300`, 4,101 hitting queries produced **804,341** alignment rows
    (194 MB of TSV, each row carrying `cigar`/`qaln`/`taln`). Parsing that dominated arda's peak RSS --
    877 MB, against 284 MB for the mmseqs subprocess itself. Reducing before `convertalis` writes
    1.0 MB instead, and costs 0.04 s.
    """
    run(["filterdb", str(result_db), str(out_db), "--extract-lines", "1"])
    return Path(out_db)


def convertalis(
    query_db: str | Path,
    target_db: str | Path,
    result_db: str | Path,
    out_tsv: str | Path,
    *,
    format_output: str = DEFAULT_FORMAT_OUTPUT,
    threads: int = 1,
    search_type: int | None = None,
) -> Path:
    """Convert an alignment result DB to a TSV with the requested columns.

    ``search_type`` must be passed for nucleotide results (3) so convertalis can
    interpret the alignment; otherwise mmseqs cannot tell nt from translated.
    """
    args = [
        "convertalis", str(query_db), str(target_db), str(result_db), str(out_tsv),
        "--format-output", format_output,
        "--threads", str(threads),
    ]
    if search_type is not None:
        args += ["--search-type", str(search_type)]
    run(args)
    return Path(out_tsv)


def easy_search(
    query_fasta: str | Path,
    target_fasta_or_db: str | Path,
    out_tsv: str | Path,
    tmp_dir: str | Path,
    *,
    search_type: int = SEARCH_TYPE_AUTO,
    sensitivity: float = 5.7,
    evalue: float = 1e-3,
    max_seqs: int = 300,
    threads: int = 1,
    format_output: str = DEFAULT_FORMAT_OUTPUT,
    strand: int | None = None,
    extra: list[str] | None = None,
) -> Path:
    """One-shot createdb+search+convertalis producing a TSV.

    ``strand`` (nucleotide search only): 1 forward, 2 both strands; ``None`` lets
    mmseqs default (forward).
    """
    args = [
        "easy-search", str(query_fasta), str(target_fasta_or_db),
        str(out_tsv), str(tmp_dir),
        "--search-type", str(search_type),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "--threads", str(threads),
        "--format-output", format_output,
    ]
    if strand is not None:
        args += ["--strand", str(strand)]
    if extra:
        args += extra
    run(args)
    return Path(out_tsv)
