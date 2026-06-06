"""Wrapper around the downloaded IgBLAST release.

Used only at *build time* (Phase 1) to construct the curated reference DB; the
runtime annotator does not depend on IgBLAST.

The IgBLAST release is expected under ``<project>/bin`` (placed there by
``setup.sh``), laid out as::

    bin/
      igblastn  igblastp  makeblastdb  edit_imgt_file.pl
      internal_data/   optional_file/

``$IGDATA`` is pointed at ``bin/`` so IgBLAST finds ``internal_data`` and the
per-organism ``optional_file/<organism>_gl.aux`` auxiliary files.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from .paths import bin_dir

__all__ = [
    "IgBlastError",
    "igdata_env",
    "tool",
    "edit_imgt_file",
    "makeblastdb",
    "igblastn_airr",
    "SUPPORTED_ORGANISMS",
]

# Organisms shipping with IgBLAST internal_data that arda builds DBs for.
SUPPORTED_ORGANISMS = ("human", "mouse", "rat", "rabbit", "rhesus_monkey")


class IgBlastError(RuntimeError):
    """Raised when an IgBLAST tool invocation fails or is missing."""


def tool(name: str) -> Path:
    """Resolve an IgBLAST tool path under ``bin/``."""
    p = bin_dir() / name
    if not p.exists():
        raise IgBlastError(
            f"IgBLAST tool {name!r} not found at {p}. Run setup.sh to download "
            "the IgBLAST release."
        )
    return p


def has_internal_annotation(organism: str, group: str) -> bool:
    """Whether IgBLAST ships V-region annotation for this organism + group.

    IG uses the generic ``<org>_V`` database; TR needs ``<org>_TR_V``, which only
    human and mouse ship. Missing annotation means IgBLAST cannot assign FR/CDR
    regions for that group, so the locus must be skipped during the build.
    """
    stem = f"{organism}_TR_V" if group == "TR" else f"{organism}_V"
    return (bin_dir() / "internal_data" / organism / f"{stem}.nin").exists()


def igdata_env() -> dict[str, str]:
    """Environment with ``IGDATA`` pointing at the IgBLAST data root."""
    env = dict(os.environ)
    env["IGDATA"] = str(bin_dir())
    return env


def _run(cmd: list[str], *, stdout_path: Path | None = None) -> subprocess.CompletedProcess:
    out = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    try:
        proc = subprocess.run(
            list(map(str, cmd)),
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
            env=igdata_env(),
        )
    finally:
        if stdout_path:
            out.close()
    if proc.returncode != 0:
        raise IgBlastError(
            f"`{' '.join(map(str, cmd))}` failed (exit {proc.returncode}):\n{proc.stderr}"
        )
    return proc


def edit_imgt_file(imgt_fasta: str | Path, out_fasta: str | Path) -> Path:
    """Ungap an IMGT germline FASTA via ``edit_imgt_file.pl``."""
    _run(["perl", tool("edit_imgt_file.pl"), str(imgt_fasta)],
         stdout_path=Path(out_fasta))
    return Path(out_fasta)


def makeblastdb(in_fasta: str | Path, out_db: str | Path, *, dbtype: str = "nucl") -> Path:
    """Build a germline BLAST database from an ungapped FASTA."""
    _run([
        tool("makeblastdb"),
        "-in", str(in_fasta),
        "-parse_seqids",
        "-dbtype", dbtype,
        "-out", str(out_db),
    ])
    return Path(out_db)


def igblastn_airr(
    query_fasta: str | Path,
    out_tsv: str | Path,
    *,
    organism: str,
    germline_db_v: str | Path,
    germline_db_j: str | Path,
    germline_db_d: str | Path | None = None,
    auxiliary_data: str | Path | None = None,
    ig_seqtype: str = "TCR",
    num_threads: int = 1,
) -> Path:
    """Run ``igblastn -outfmt 19`` (AIRR rearrangement TSV)."""
    cmd = [
        tool("igblastn"),
        "-germline_db_V", str(germline_db_v),
        "-germline_db_J", str(germline_db_j),
        "-organism", organism,
        "-ig_seqtype", ig_seqtype,
        "-query", str(query_fasta),
        "-outfmt", "19",
        "-num_threads", str(num_threads),
    ]
    if germline_db_d is not None:
        cmd += ["-germline_db_D", str(germline_db_d)]
    else:
        # igblastn still wants a D db arg for VDJ; callers pass one for VDJ chains.
        pass
    if auxiliary_data is not None:
        cmd += ["-auxiliary_data", str(auxiliary_data)]
    _run(cmd, stdout_path=Path(out_tsv))
    return Path(out_tsv)
