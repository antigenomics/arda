"""Auto-fetch the curated arda reference database (a GitHub release asset).

The wheel ships code only; the curated ``vdj/`` reference (allele FASTAs + region markup, per
species, AA + NT; ~50 MB on disk, ~3 MB compressed) is published as the
``arda-reference-vdj.tar.gz`` asset on the
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
from ._locking import build_lock
from ._mmseqs_fetch import _download  # reuse the hardened UA/retry downloader

_ASSET = "arda-reference-vdj.tar.gz"
_URL = "https://github.com/antigenomics/arda/releases/download/v{version}/" + _ASSET


def reference_url(version: str | None = None) -> str:
    """Release-asset URL for the reference tarball at ``version`` (default: this build)."""
    return _URL.format(version=version or __version__)


def fetch_database(dest: Path, *, force: bool = False, version: str | None = None) -> Path:
    """Download + extract the reference (``vdj/<org>/...``) into ``dest``; return ``dest``.

    Skips the download if ``dest/vdj`` already exists (unless ``force``).

    **``vdj/`` must never be visible in a partial state**, because its mere existence is the gate
    every other arda process tests (``paths.database_dir``) — a half-populated tree is not a slow
    download, it is a silent wrong answer. Two things guarantee that:

    * a :func:`~arda._locking.build_lock`, because arda is routinely run concurrently against the
      same cache (one Nextflow process per sample, one SLURM task per array index), and on first use
      in a fresh environment all of them would otherwise fetch into the same path at once;
    * extraction into ``dest`` itself, then a single :func:`os.replace`. The extraction directory has
      to be a *sibling* of the target for that to be a rename: the previous code extracted into
      ``/tmp`` and called ``shutil.move``, which silently degrades to a recursive **copy** across
      filesystems (``/tmp`` and ``~/.cache`` usually are different ones) — populating ``vdj/`` file
      by file, in full view of every other process.

    Extraction still rejects symlinks/hardlinks and any path escaping the staging dir.
    """
    dest = Path(dest)
    target = dest / "vdj"
    if target.is_dir() and not force:
        return dest
    dest.mkdir(parents=True, exist_ok=True)

    # Under --force the work is ours by definition; otherwise whoever installed vdj/ has done it.
    done = (lambda: False) if force else target.is_dir
    with build_lock(dest / ".vdj.lock", done=done) as ours:
        if not ours:
            return dest
        _fetch_into(dest, target, version)
    return dest


def _fetch_into(dest: Path, target: Path, version: str | None) -> None:
    """Download, verify and stage the reference beside ``target``, then swap it in atomically."""
    url = reference_url(version)
    print(f"[arda] fetching reference database (one-time): {url}", file=sys.stderr)
    # dir=dest, not /tmp: os.replace below is only a rename -- and only atomic -- within one filesystem.
    staging = Path(tempfile.mkdtemp(dir=dest, prefix=".arda_db_"))
    try:
        tarball = staging / _ASSET
        _download(url, tarball)
        with tarfile.open(tarball) as tf:
            staged = staging.resolve()
            for member in tf.getmembers():
                if member.issym() or member.islnk():
                    raise RuntimeError(f"Refusing to extract link from archive: {member.name}")
                member_path = (staging / member.name).resolve()
                if not str(member_path).startswith(str(staged) + os.sep):
                    raise RuntimeError(f"Refusing to extract path outside temp dir: {member.name}")
            tf.extractall(staging)
        # tarball root holds vdj/ (created with `tar -C database ... vdj`); tolerate a nested layout.
        src = staging / "vdj"
        if not src.is_dir():
            src = next((p for p in staging.rglob("vdj") if p.is_dir()), None)
        if src is None:
            raise RuntimeError(f"Unexpected reference archive layout in {_ASSET} (no vdj/)")

        if target.exists():
            # os.replace refuses a non-empty destination directory, so retire the old tree first.
            # Renaming it (rather than deleting it) keeps a reader that is mid-search on the old
            # files alive: it holds inodes, and those survive until it closes them.
            retired = dest / f".vdj.old.{os.getpid()}"
            shutil.rmtree(retired, ignore_errors=True)
            os.replace(target, retired)
            try:
                os.replace(src, target)
            finally:
                shutil.rmtree(retired, ignore_errors=True)
        else:
            os.replace(src, target)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    print(f"[arda] installed reference -> {dest}", file=sys.stderr)
