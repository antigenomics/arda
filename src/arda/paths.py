"""Filesystem layout discovery for arda.

Resolves the locations of the project's ``bin/`` (downloaded IgBLAST release),
``data/`` (scratch/downloads), and ``database/`` (committed curated references).

Resolution order for the project root:

1. ``$ARDA_HOME`` if set.
2. Walking up from this file to find a directory containing both ``database``
   and one of ``pyproject.toml`` / ``setup.sh`` (the editable-install / source case).
"""

from __future__ import annotations

import os
from pathlib import Path
from functools import lru_cache

__all__ = [
    "project_root",
    "bin_dir",
    "data_dir",
    "database_dir",
    "vdj_dir",
    "mhc_dir",
]

_MARKERS = ("pyproject.toml", "setup.sh")


@lru_cache(maxsize=1)
def project_root() -> Path:
    """Return the arda project root directory.

    Honors ``$ARDA_HOME``; otherwise searches parent directories for a source
    checkout. Raises ``RuntimeError`` if neither resolves.
    """
    env = os.environ.get("ARDA_HOME")
    if env:
        root = Path(env).expanduser().resolve()
        if not root.is_dir():
            raise RuntimeError(f"ARDA_HOME={env!r} is not a directory")
        return root

    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "database").is_dir() and any(
            (parent / m).exists() for m in _MARKERS
        ):
            return parent
    raise RuntimeError(
        "Could not locate the arda project root. Set $ARDA_HOME to the checkout "
        "directory (the one containing 'database/')."
    )


def bin_dir() -> Path:
    """Directory holding the downloaded IgBLAST release (gitignored)."""
    return project_root() / "bin"


def data_dir() -> Path:
    """Scratch directory for downloads and intermediate artifacts (gitignored)."""
    return project_root() / "data"


def database_dir() -> Path:
    """Committed curated reference database root."""
    return project_root() / "database"


def vdj_dir(species: str | None = None) -> Path:
    """``database/vdj`` (or a single species subdirectory)."""
    base = database_dir() / "vdj"
    return base / species if species else base


def mhc_dir(species: str | None = None) -> Path:
    """``database/mhc`` (or a single species subdirectory)."""
    base = database_dir() / "mhc"
    return base / species if species else base
