"""Auto-fetch the curated arda reference database (a GitHub release asset).

The wheel ships code only; the ~50 MB curated ``vdj/`` reference (allele FASTAs + region
markup, per species, AA + NT) is published as the ``arda-reference-vdj.tar.gz`` asset on the
matching ``vX.Y.Z`` GitHub release and downloaded once into the per-user cache on first use.
Version-sensitive precompiled mmseqs DBs are *not* shipped — they are built on demand from
the fetched FASTAs into ``<cache>/data``.
"""

from __future__ import annotations

import os
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

from . import __version__
from ._mmseqs_fetch import _download  # reuse the hardened UA/retry downloader

_ASSET = "arda-reference-vdj.tar.gz"
_URL = "https://github.com/antigenomics/arda/releases/download/v{version}/" + _ASSET


def reference_url(version: str | None = None) -> str:
    """Release-asset URL for the reference tarball at ``version`` (default: this build)."""
    return _URL.format(version=version or __version__)


def fetch_database(dest: Path, *, force: bool = False, version: str | None = None) -> Path:
    """Download + extract the reference (``vdj/<org>/...``) into ``dest``; return ``dest``.

    Skips the download if ``dest/vdj`` already exists (unless ``force``). Extraction rejects
    symlinks/hardlinks and any path that escapes the temp dir (same guards as the mmseqs
    fetch), then moves the ``vdj`` tree into place atomically-ish.
    """
    dest = Path(dest)
    if (dest / "vdj").is_dir() and not force:
        return dest
    dest.mkdir(parents=True, exist_ok=True)
    url = reference_url(version)
    print(f"[arda] fetching reference database (one-time): {url}", file=sys.stderr)
    with tempfile.TemporaryDirectory(prefix="arda_db_") as td:
        tmp = Path(td)
        tarball = tmp / _ASSET
        _download(url, tarball)
        with tarfile.open(tarball) as tf:
            tmp_resolved = tmp.resolve()
            for member in tf.getmembers():
                if member.issym() or member.islnk():
                    raise RuntimeError(f"Refusing to extract link from archive: {member.name}")
                member_path = (tmp / member.name).resolve()
                if not str(member_path).startswith(str(tmp_resolved) + os.sep):
                    raise RuntimeError(f"Refusing to extract path outside temp dir: {member.name}")
            tf.extractall(tmp)
        # tarball root holds vdj/ (created with `tar -C database ... vdj`); tolerate a nested layout.
        src = tmp / "vdj"
        if not src.is_dir():
            src = next((p for p in tmp.rglob("vdj") if p.is_dir()), None)
        if src is None:
            raise RuntimeError(f"Unexpected reference archive layout in {_ASSET} (no vdj/)")
        target = dest / "vdj"
        if target.exists():
            shutil.rmtree(target)
        shutil.move(str(src), str(target))
    print(f"[arda] installed reference -> {dest}", file=sys.stderr)
    return dest
