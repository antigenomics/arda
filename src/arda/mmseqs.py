"""Thin wrapper around the ``mmseqs`` binary.

Inspired by pymmseqs (MIT) but deliberately dependency-free: we only need
binary discovery, a subprocess runner, and the ``createdb`` / ``search`` /
``convertalis`` (and ``easy-search``) pipeline used by the annotator.

Discovery order for the binary: ``$ARDA_MMSEQS`` → ``<project>/bin/mmseqs`` →
``mmseqs`` on ``PATH``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from functools import lru_cache

from .paths import bin_dir

__all__ = [
    "MMseqsError",
    "mmseqs_binary",
    "run",
    "createdb",
    "search",
    "convertalis",
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
    """Locate the ``mmseqs`` executable."""
    env = os.environ.get("ARDA_MMSEQS")
    if env:
        return env
    local = bin_dir() / "mmseqs"
    if local.exists():
        return str(local)
    found = shutil.which("mmseqs")
    if found:
        return found
    raise MMseqsError(
        "mmseqs binary not found. Install it (conda install -c bioconda mmseqs2) "
        "or set $ARDA_MMSEQS."
    )


def run(args: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Run ``mmseqs <args>`` capturing stdout/stderr."""
    cmd = [mmseqs_binary(), *map(str, args)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise MMseqsError(
            f"`{' '.join(cmd)}` failed (exit {proc.returncode}):\n{proc.stderr}"
        )
    return proc


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
    extra: list[str] | None = None,
) -> Path:
    """Run ``mmseqs search`` with backtrace enabled (``-a``)."""
    args = [
        "search", str(query_db), str(target_db), str(result_db), str(tmp_dir),
        "--search-type", str(search_type),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "--threads", str(threads),
        "-a",  # keep backtrace so convertalis can emit cigar/qaln/taln
    ]
    if extra:
        args += extra
    run(args)
    return Path(result_db)


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
