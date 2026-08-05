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

# The release tag that carries the CURRENT reference asset — deliberately NOT `__version__`.
#
# The reference is data and changes on its own schedule; the package version changes every
# release. Deriving the URL from `__version__` silently coupled them, so the moment the version
# was bumped ahead of a GitHub release every cold-cache install died on:
#     RuntimeError: failed to download .../v2.6.0/arda-reference-vdj.tar.gz: HTTP Error 404
# Found by running eight concurrent arda processes against a fresh cache on a cluster; it would
# have hit the first real `pip install` of any release whose asset was not up yet.
#
# Bump this ONLY when a release actually publishes a new reference tarball. `$ARDA_REFERENCE_TAG`
# overrides it for testing a candidate asset.
_REFERENCE_TAG = "2.5.7"


def reference_url(version: str | None = None) -> str:
    """Release-asset URL for the reference tarball.

    Defaults to :data:`_REFERENCE_TAG` (the release that published the current reference), not
    to the running package version — see the note there.
    """
    tag = version or os.environ.get("ARDA_REFERENCE_TAG") or _REFERENCE_TAG
    return _URL.format(version=tag.lstrip("v"))


def _candidate_urls(version: str | None = None) -> list[str]:
    """Reference URLs to try, in order.

    The pinned tag first. If an explicit ``version`` was asked for it is the only candidate —
    a caller naming a version must not silently get a different one. Otherwise fall back to the
    running package version, which covers a release that *did* ship its own asset before
    ``_REFERENCE_TAG`` was updated.
    """
    if version:
        return [reference_url(version)]
    urls = [reference_url()]
    fallback = _URL.format(version=__version__)
    if fallback not in urls:
        urls.append(fallback)
    return urls


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
    urls = _candidate_urls(version)
    print(f"[arda] fetching reference database (one-time): {urls[0]}", file=sys.stderr)
    # dir=dest, not /tmp: os.replace below is only a rename -- and only atomic -- within one filesystem.
    staging = Path(tempfile.mkdtemp(dir=dest, prefix=".arda_db_"))
    try:
        tarball = staging / _ASSET
        last: Exception | None = None
        for i, url in enumerate(urls):
            try:
                _download(url, tarball)
                last = None
                break
            except Exception as exc:  # noqa: BLE001 — try the next candidate tag
                last = exc
                if i + 1 < len(urls):
                    print(f"[arda] {url} unavailable; trying {urls[i + 1]}", file=sys.stderr)
        if last is not None:
            raise last
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
