"""Filesystem layout discovery for arda.

Two modes, resolved transparently:

* **Source checkout** (``$ARDA_HOME`` or an editable/development install) — ``bin/``,
  ``data/`` and ``database/`` live next to ``pyproject.toml``, exactly as committed.
* **PyPI install** (``pip install arda-mapper``, no source tree) — there is no bundled
  ``database/`` in the wheel (it is 50+ MB of curated references), so everything lives under
  a per-user cache (``$XDG_CACHE_HOME/arda`` or ``~/.cache/arda``). The curated reference is
  **auto-fetched once** into the cache on first use (see :mod:`arda._database_fetch`), and
  mmseqs target DBs are built there on demand. No ``ARDA_HOME`` needed.

Set ``$ARDA_NO_AUTO_FETCH`` to disable the one-time reference download (e.g. air-gapped
runs where the cache was pre-populated).
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

__all__ = [
    "project_root",
    "cache_root",
    "bin_dir",
    "data_dir",
    "database_dir",
    "vdj_dir",
]

_MARKERS = ("pyproject.toml", "setup.sh")
_REFERENCE_MARKER = "vdj"  # database_dir()/vdj must exist for the reference to be usable


@lru_cache(maxsize=1)
def _source_root() -> Path | None:
    """The source checkout root if this is a dev/editable install, else ``None``.

    Honors ``$ARDA_HOME`` (must be a directory); otherwise walks up from this file looking
    for a directory that has both ``database/`` and a project marker. Never raises for a
    plain PyPI install — it simply returns ``None`` and the caller falls back to the cache.
    """
    env = os.environ.get("ARDA_HOME")
    if env:
        root = Path(env).expanduser().resolve()
        if not root.is_dir():
            raise RuntimeError(f"ARDA_HOME={env!r} is not a directory")
        return root
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "database").is_dir() and any((parent / m).exists() for m in _MARKERS):
            return parent
    return None


@lru_cache(maxsize=1)
def cache_root() -> Path:
    """Per-user cache root (``$XDG_CACHE_HOME/arda`` or ``~/.cache/arda``)."""
    base = os.environ.get("XDG_CACHE_HOME")
    root = (Path(base) if base else Path.home() / ".cache") / "arda"
    root.mkdir(parents=True, exist_ok=True)
    return root


@lru_cache(maxsize=1)
def project_root() -> Path:
    """Source checkout root if present, else the per-user cache root.

    Unlike older arda, this never raises: a PyPI install with no source tree resolves to the
    cache dir, where the reference is auto-fetched and mmseqs DBs are built.
    """
    return _source_root() or cache_root()


def bin_dir() -> Path:
    """Directory holding the downloaded mmseqs/IgBLAST binaries (gitignored / cached)."""
    return project_root() / "bin"


def data_dir() -> Path:
    """Writable scratch directory for downloads and built mmseqs DBs (gitignored / cached)."""
    return project_root() / "data"


@lru_cache(maxsize=1)
def database_dir() -> Path:
    """Curated reference database root.

    Source checkout → the committed ``database/``. PyPI install → ``<cache>/database``, with
    the reference auto-fetched on first call unless ``$ARDA_NO_AUTO_FETCH`` is set.
    """
    src = _source_root()
    if src is not None:
        return src / "database"
    db = cache_root() / "database"
    if not (db / _REFERENCE_MARKER).is_dir() and "ARDA_NO_AUTO_FETCH" not in os.environ:
        from ._database_fetch import fetch_database

        fetch_database(db)
    return db


def vdj_dir(species: str | None = None) -> Path:
    """``database/vdj`` (or a single species subdirectory)."""
    base = database_dir() / "vdj"
    return base / species if species else base
