"""Auto-fetch the curated arda reference database (a GitHub release asset).

The wheel ships code only; the curated ``vdj/`` reference (allele FASTAs + region markup, per
species, AA + NT; ~50 MB on disk, ~3 MB compressed) is published as the
``arda-reference-vdj.tar.gz`` asset on the
versioned data release (independent of the package version) and downloaded once into the
per-user cache on first use. Readers pin a verified immutable generation.
Version-sensitive precompiled mmseqs DBs are *not* shipped — they are built on demand from
the fetched FASTAs into ``<cache>/data``.
"""

from __future__ import annotations

import os
import hashlib
import json
import uuid
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
_PINNED_FILES = json.loads(Path(__file__).with_name("_reference_manifest.json").read_text())["files_sha256"]


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
    if version or os.environ.get("ARDA_REFERENCE_TAG"):
        return [reference_url(version)]
    urls = [reference_url()]
    fallback = _URL.format(version=__version__)
    if fallback not in urls:
        urls.append(fallback)
    return urls


_ORGANISMS = ("human", "mouse", "rat", "rabbit", "rhesus_monkey")
_REQUIRED = ("alleles.fasta", "alleles.aa.fasta", "markup.tsv", "markup.aa.tsv",
             "cdr3_anchors.tsv", "combinations.tsv", "d_germlines.fasta")


def _validate_required(root: Path) -> None:
    for organism in _ORGANISMS:
        files = _REQUIRED + (("d_prior.tsv",) if organism in ("human", "mouse") else ())
        for name in files:
            path = root / "vdj" / organism / name
            minimum = 1000 if name in _REQUIRED[:5] else 0
            if not path.is_file() or path.stat().st_size <= minimum:
                raise RuntimeError(f"Incomplete reference: missing or empty {organism}/{name}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def installed_database(dest: Path, *, version: str | None = None) -> Path | None:
    """Resolve one complete immutable generation, rejecting missing/corrupt cache files."""
    try:
        name = json.loads((dest / "REFERENCE.json").read_text())["generation"]
        if not isinstance(name, str) or len(name) != 32 or any(c not in "0123456789abcdef" for c in name):
            return None
        root = dest / ".references" / name
        manifest = json.loads((root / "MANIFEST.json").read_text())
        files = manifest["files_sha256"]
        if manifest["source_url"] not in _candidate_urls(version):
            return None
        _validate_required(root)
        required = {f"vdj/{o}/{f}" for o in _ORGANISMS for f in _REQUIRED}
        required |= {f"vdj/{o}/d_prior.tsv" for o in ("human", "mouse")}
        if not isinstance(files, dict) or not required <= files.keys():
            return None
        if (manifest["source_url"] == reference_url(_REFERENCE_TAG) and
                any(files.get(name) != digest for name, digest in _PINNED_FILES.items())):
            return None
        for relative, expected in files.items():
            path = root / relative
            if not path.resolve().is_relative_to(root.resolve()) or _sha256(path) != expected:
                return None
        return root
    except (OSError, ValueError, KeyError, TypeError, RuntimeError):
        return None


def legacy_database(dest: Path, *, version: str | None = None) -> Path | None:
    """Accept a flat legacy cache only when all canonical data hashes match."""
    if not _PINNED_FILES or reference_url(_REFERENCE_TAG) not in _candidate_urls(version):
        return None
    try:
        _validate_required(dest)
        return dest if all(_sha256(dest / name) == digest
                           for name, digest in _PINNED_FILES.items()) else None
    except (OSError, RuntimeError):
        return None


def fetch_database(dest: Path, *, force: bool = False, version: str | None = None) -> Path:
    """Fetch a complete reference and return its immutable root (containing ``vdj/``).

    Publication replaces a small pointer atomically. Existing readers keep their old root;
    failed or interrupted refreshes cannot remove it. A warm valid cache needs no network.
    """
    dest = Path(dest)
    current = installed_database(dest, version=version)
    if current is not None and not force:
        return current
    dest.mkdir(parents=True, exist_ok=True)
    with build_lock(dest / ".vdj.lock",
                    done=lambda: not force and installed_database(dest, version=version) is not None) as ours:
        if not ours:
            current = installed_database(dest, version=version)
            if current is None:
                raise RuntimeError("Reference changed while resolving its installed generation")
            return current
        if not force and legacy_database(dest, version=version) is not None:
            with tempfile.TemporaryDirectory(dir=dest, prefix=".arda_legacy_") as td:
                stage = Path(td)
                for name in _PINNED_FILES:
                    target = stage / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(dest / name, target)
                return _publish_reference(dest, stage / "vdj", reference_url(_REFERENCE_TAG))
        return _fetch_into(dest, version)


def _fetch_into(dest: Path, version: str | None) -> Path:
    urls = _candidate_urls(version)
    print(f"[arda] fetching reference database (one-time): {urls[0]}", file=sys.stderr)
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
                if not (member.isfile() or member.isdir()):
                    raise RuntimeError(f"Refusing to extract link from archive: {member.name}")
                if not (staging / member.name).resolve().is_relative_to(staged):
                    raise RuntimeError(f"Refusing to extract path outside temp dir: {member.name}")
            tf.extractall(staging)
        src = staging / "vdj"
        if not src.is_dir():
            candidates = [p for p in staging.rglob("vdj") if p.is_dir()]
            if len(candidates) != 1:
                raise RuntimeError(f"Unexpected reference archive layout in {_ASSET}")
            src = candidates[0]
        _validate_required(src.parent)
        if url == reference_url(_REFERENCE_TAG):
            for name, expected in _PINNED_FILES.items():
                if _sha256(src.parent / name) != expected:
                    raise RuntimeError(f"Reference checksum mismatch: {name}")
        return _publish_reference(dest, src, url)
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _publish_reference(dest: Path, src: Path, source_url: str) -> Path:
    """Install one private tree and publish only its pointer; never retire an active root."""
    prepared = Path(tempfile.mkdtemp(dir=dest, prefix=".arda_ref_"))
    pointer = dest / (".reference-" + uuid.uuid4().hex + ".json")
    try:
        os.replace(src, prepared / "vdj")
        files = {p.relative_to(prepared).as_posix(): _sha256(p)
                 for p in (prepared / "vdj").rglob("*") if p.is_file()}
        if (source_url == reference_url(_REFERENCE_TAG) and
                any(files.get(name) != digest for name, digest in _PINNED_FILES.items())):
            raise RuntimeError("Reference checksum mismatch while staging the selected data release")
        manifest = {"files_sha256": files, "source_url": source_url}
        (prepared / "MANIFEST.json").write_text(json.dumps(manifest, sort_keys=True))
        generations = dest / ".references"
        generations.mkdir(exist_ok=True)
        name = uuid.uuid4().hex
        root = generations / name
        os.replace(prepared, root)
        pointer.write_text(json.dumps({"generation": name}))
        os.replace(pointer, dest / "REFERENCE.json")
    finally:
        pointer.unlink(missing_ok=True)
        shutil.rmtree(prepared, ignore_errors=True)
    print(f"[arda] installed reference -> {root}", file=sys.stderr)
    return root
