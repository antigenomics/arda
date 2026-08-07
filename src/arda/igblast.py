"""Wrapper around an IgBLAST release, fetched on demand.

Used at *build time* (Phase 1) to construct the curated reference DB, and by ``arda igblast``,
which is how every gold-standard comparison in the benchmark is produced. The runtime annotator
does not depend on IgBLAST.

The release is a flat directory -- executables plus the ``internal_data`` and ``optional_file``
trees, with ``$IGDATA`` pointed at it::

    igblastn  igblastp  makeblastdb  edit_imgt_file.pl
    internal_data/   optional_file/

Resolved by :func:`igblast_root`, in order:

1. ``$ARDA_IGBLAST`` -- an explicit directory, never fetched;
2. ``<project>/bin`` if it already holds an IgBLAST (what ``setup.sh`` produces in a checkout);
3. ``<project>/igblast``, **auto-fetched from NCBI on first use**.

Step 3 is why this module changed shape. It used to resolve only through ``<project>/bin``, so a
plain ``pip install arda-mapper`` -- which has no checkout and never runs ``setup.sh`` -- could
not run ``arda igblast`` at all. The failure was also misleading: it surfaced as
``IgBlastError: IgBLAST ships no internal annotation for organism 'human'``, which names a
missing data file rather than a missing install, and so reads as a broken reference rather than
as "IgBLAST was never installed here".
"""

from __future__ import annotations

import os
import subprocess
from functools import lru_cache
from pathlib import Path

from ._igblast_fetch import fetch, installed_version
from .paths import bin_dir, project_root

__all__ = [
    "IgBlastError",
    "igblast_root",
    "igblast_version",
    "igdata_env",
    "tool",
    "edit_imgt_file",
    "makeblastdb",
    "igblastn_airr",
    "auxiliary_data",
    "SUPPORTED_ORGANISMS",
]

# Organisms shipping with IgBLAST internal_data that arda builds DBs for.
SUPPORTED_ORGANISMS = ("human", "mouse", "rat", "rabbit", "rhesus_monkey")


class IgBlastError(RuntimeError):
    """Raised when an IgBLAST tool invocation fails or is missing."""


def auxiliary_data(organism: str) -> Path:
    """``optional_file/<organism>_gl.aux`` — IgBLAST's J-gene coding frames.

    ⛔ **Without this file IgBLAST silently emits no CDR3 and no junction.** It is what tells
    igblastn each J allele's reading frame, and with no frame there is nothing to place the
    Phe/Trp 118 anchor against. Everything else still works: V and J are called, `v_score` is
    normal, the process exits 0 — only `cdr3*`, `junction` and `junction_aa` come back empty, on
    every read.

    It lives beside the executables, under :func:`igblast_root`. Both callers used to look under
    ``paths.bin_dir()`` and then fall back to passing nothing when the file was not there. Those
    two are THE SAME directory in a source checkout (``setup.sh`` installs IgBLAST into
    ``<repo>/bin``), so it worked everywhere it was developed and failed on every auto-fetched
    install, where the root is ``$XDG_CACHE_HOME/arda/igblast`` while ``bin_dir()`` is
    ``$XDG_CACHE_HOME/arda/bin``. Measured cost: a 10,000-read amplicon IgBLAST truth carrying
    `j_call` on 9,070 of 9,300 reads and `junction_aa` on **zero** — written up as an IgBLAST
    limitation at 151 bp before it was traced to here.

    Raises:
        IgBlastError: if the file is absent. A missing frame table must not degrade to "no
            junctions": that is indistinguishable from a truth set which genuinely has none, which
            is exactly how this went unnoticed.
    """
    aux = igblast_root() / "optional_file" / f"{organism}_gl.aux"
    if not aux.exists():
        raise IgBlastError(
            f"IgBLAST auxiliary data not found at {aux}. Without it igblastn emits no CDR3 and no "
            f"junction for any read, and reports no error. It ships inside the IgBLAST release: "
            f"re-fetch it, or point $ARDA_IGBLAST at a directory containing optional_file/."
        )
    return aux


@lru_cache(maxsize=1)
def igblast_root() -> Path:
    """The directory holding the IgBLAST executables and ``internal_data``.

    Fetches a release on first use; ``$ARDA_NO_AUTO_FETCH`` turns that into an error instead.
    """
    override = os.environ.get("ARDA_IGBLAST")
    if override:
        root = Path(override)
        if not (root / "igblastn").exists():
            raise IgBlastError(
                f"$ARDA_IGBLAST={override} does not contain igblastn. Point it at the directory "
                "holding the IgBLAST executables and internal_data/."
            )
        return root
    if (bin_dir() / "igblastn").exists():  # setup.sh layout, in a source checkout
        return bin_dir()
    return fetch(project_root() / "igblast")  # returns immediately if already complete


def igblast_version() -> str | None:
    """The auto-fetched IgBLAST release version, or None if it was not fetched by arda."""
    return installed_version(igblast_root())


def tool(name: str) -> Path:
    """Resolve an IgBLAST executable, fetching the release if it is not installed yet."""
    p = igblast_root() / name
    if not p.exists():
        raise IgBlastError(
            f"IgBLAST tool {name!r} not found at {p}. The IgBLAST release at "
            f"{igblast_root()} looks incomplete; remove it to force a re-fetch, or point "
            "$ARDA_IGBLAST at a complete one."
        )
    return p


def has_internal_annotation(organism: str, group: str) -> bool:
    """Whether IgBLAST ships V-region annotation for this organism + group.

    IG uses the generic ``<org>_V`` database; TR needs ``<org>_TR_V``, which only
    human and mouse ship. Missing annotation means IgBLAST cannot assign FR/CDR
    regions for that group, so the locus must be skipped during the build.
    """
    stem = f"{organism}_TR_V" if group == "TR" else f"{organism}_V"
    return (igblast_root() / "internal_data" / organism / f"{stem}.nin").exists()


def igdata_env() -> dict[str, str]:
    """Environment with ``IGDATA`` pointing at the IgBLAST data root."""
    env = dict(os.environ)
    env["IGDATA"] = str(igblast_root())
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
