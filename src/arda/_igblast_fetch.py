"""Download an NCBI IgBLAST release into the arda cache (packaged, stdlib only).

Used by :mod:`arda.igblast` for transparent auto-fetch, and by the
``scripts/fetch_igblast.py`` CLI wrapper that ``setup.sh`` calls. Kept in the package -- not
only in ``scripts/`` -- so ``arda igblast`` works from a plain ``pip install``, which has no
source checkout and therefore never runs ``setup.sh``.

NCBI ships the release as ``bin/`` (executables) plus sibling ``internal_data/`` and
``optional_file/`` trees. arda flattens all three into one directory and points ``$IGDATA`` at
it, which is the layout ``setup.sh`` has always produced.
"""

from __future__ import annotations

import os
import json
import uuid
import platform
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

from ._locking import build_lock
from ._mmseqs_fetch import _download

LATEST_URL = "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/"

_UA = "Mozilla/5.0 (arda igblast-fetch)"

# Written last, inside the staging tree, so it also records WHICH release a result came from.
VERSION_FILE = "IGBLAST_VERSION"


def platform_suffix() -> str:
    """The NCBI asset suffix for this platform.

    Honors ``$ARDA_IGBLAST_ASSET`` for the same reason ``$ARDA_MMSEQS_ASSET`` exists: a
    cross-build, or a host whose CPU the default asset does not suit.
    """
    override = os.environ.get("ARDA_IGBLAST_ASSET")
    if override:
        return override
    system = platform.system()
    if system == "Linux":
        return "x64-linux"
    if system == "Darwin":
        # NCBI ships only an x64 macOS build; on Apple Silicon it runs under Rosetta 2.
        return "x64-macosx"
    if system == "Windows":
        return "x64-win64"
    raise RuntimeError(
        f"Unsupported platform for IgBLAST auto-fetch: {system}. Install IgBLAST manually "
        "and point $ARDA_IGBLAST at the directory holding igblastn and internal_data/."
    )


def find_tarball(suffix: str) -> str:
    """Scrape the LATEST listing for ``ncbi-igblast-<ver>-<suffix>.tar.gz``."""
    req = urllib.request.Request(LATEST_URL, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=60) as resp:
        html = resp.read().decode("utf-8", "replace")
    pat = re.compile(rf"ncbi-igblast-([0-9.]+)-{re.escape(suffix)}\.tar\.gz")
    versions = sorted(set(pat.findall(html)), key=lambda v: [int(x) for x in v.split(".")])
    if not versions:
        raise RuntimeError(f"No IgBLAST tarball for {suffix!r} at {LATEST_URL}")
    return f"ncbi-igblast-{versions[-1]}-{suffix}.tar.gz"




def _safe_extract(tarball: Path, into: Path) -> Path:
    """Extract, refusing anything that could write outside *into*; return the release root."""
    with tarfile.open(tarball) as tf:
        resolved = into.resolve()
        for member in tf.getmembers():
            if member.issym() or member.islnk():
                raise RuntimeError(f"Refusing to extract link from IgBLAST archive: {member.name}")
            if not str((into / member.name).resolve()).startswith(str(resolved) + os.sep):
                raise RuntimeError(f"Refusing to extract outside temp dir: {member.name}")
        tf.extractall(into)
    roots = [p for p in into.iterdir() if p.is_dir() and p.name.startswith("ncbi-igblast")]
    if not roots:
        raise RuntimeError(f"Unexpected IgBLAST archive layout in {tarball.name}")
    return roots[0]


def lay_out(release_root: Path, dest: Path) -> None:
    """Flatten ``bin/* + internal_data + optional_file`` from a release into *dest*."""
    dest.mkdir(parents=True, exist_ok=True)
    src_bin = release_root / "bin"
    if not src_bin.is_dir():
        raise RuntimeError(f"Unexpected IgBLAST layout: no bin/ under {release_root}")
    for item in src_bin.iterdir():
        shutil.copy2(item, dest / item.name)
        (dest / item.name).chmod(0o755)
    for tree in ("internal_data", "optional_file"):
        src = release_root / tree
        if src.is_dir():
            shutil.copytree(src, dest / tree, dirs_exist_ok=True)


def _flat_ready(root: Path) -> bool:
    required = (VERSION_FILE, "igblastn", "makeblastdb", "edit_imgt_file.pl",
                "internal_data/human/human.ndm.imgt", "optional_file/human_gl.aux")
    try:
        return all((root / name).is_file() and (root / name).stat().st_size > 0 for name in required)
    except OSError:
        return False


def installed_root(dest: Path) -> Path | None:
    """Resolve the selected immutable refresh, or an existing flat installation."""
    pointer = dest.parent / f".{dest.name}.current.json"
    if pointer.exists():
        try:
            name = json.loads(pointer.read_text())["generation"]
            if not isinstance(name, str) or len(name) != 32 or any(c not in "0123456789abcdef" for c in name):
                return None
            root = dest.parent / f".{dest.name}.generations" / name
            return root if _flat_ready(root) else None
        except (OSError, ValueError, KeyError, TypeError):
            return None
    return dest if _flat_ready(dest) else None


def is_complete(dest: Path) -> bool:
    return installed_root(dest) is not None


def fetch(dest: Path, *, force: bool = False) -> Path:
    """Install IgBLAST and return its resolved flat root. Refreshes retain old reader roots.

    Concurrency-safe by construction, because arda runs concurrently against one cache by
    design (a Nextflow process per sample, a SLURM array task per shard). The release is built
    into a **sibling** staging directory and moved into place with a single :func:`os.replace`,
    so no other process can ever observe a partially populated root -- and the staging directory
    is a sibling specifically so that the move is a rename, not a silent recursive copy across a
    filesystem boundary.
    """
    dest = Path(dest)
    current = installed_root(dest)
    if current is not None and not force:
        return current
    if os.environ.get("ARDA_NO_AUTO_FETCH"):
        raise RuntimeError(
            f"IgBLAST is not installed at {dest} and $ARDA_NO_AUTO_FETCH is set. Fetch it with "
            "`python -m arda._igblast_fetch --dest <dir>`, or point $ARDA_IGBLAST at an "
            "existing IgBLAST directory."
        )

    dest.parent.mkdir(parents=True, exist_ok=True)
    with build_lock(dest.parent / ".igblast.lock", done=lambda: is_complete(dest) and not force,
                    timeout=2100) as ours:
        if not ours:
            current = installed_root(dest)
            if current is None:
                raise RuntimeError("IgBLAST installation changed while resolving its root")
            return current
        suffix = platform_suffix()
        tarball = find_tarball(suffix)
        print(f"[arda] fetching IgBLAST (one-time): {LATEST_URL}{tarball}", file=sys.stderr)
        staging = dest.parent / f".{dest.name}.staging.{os.getpid()}"
        shutil.rmtree(staging, ignore_errors=True)
        try:
            # Scratch next to the destination, NOT in the system temp dir. The release is a
            # ~400 MB download that extracts to more, and on a cluster /tmp is routinely a small
            # node-local filesystem -- aldan3's is 2.0 GB with 29 MB free, where this died on
            # `OSError: [Errno 28] No space left on device`. Putting it here also means the whole
            # download-extract-lay out-replace sequence happens on one filesystem, so the final
            # move is a rename by construction rather than by assumption.
            with tempfile.TemporaryDirectory(prefix=".arda_igblast_", dir=dest.parent) as td:
                tmp = Path(td)
                tar_path = tmp / tarball
                _download(LATEST_URL + tarball, tar_path, max_seconds=1800)
                lay_out(_safe_extract(tar_path, tmp), staging)
            # Last, so `is_complete` cannot be true before everything above landed. The version
            # is read back out of the asset name rather than stripped by suffix -- a strip is
            # silently a no-op when the two disagree, and would then record the whole filename
            # as the version.
            m = re.search(r"ncbi-igblast-([0-9.]+)-", tarball)
            (staging / VERSION_FILE).write_text((m.group(1) if m else tarball) + "\n")
            if not _flat_ready(staging):
                raise RuntimeError("Incomplete IgBLAST runtime: required tools or anchor tables are missing")
            if dest.exists() or (dest.parent / f".{dest.name}.current.json").exists():
                # Keep old readers and sibling tools untouched, even on failed publication.
                generations = dest.parent / f".{dest.name}.generations"
                generations.mkdir(exist_ok=True)
                name = uuid.uuid4().hex
                root = generations / name
                os.replace(staging, root)
                pointer = dest.parent / f".{dest.name}.current.json"
                pending = dest.parent / (".igblast-pointer-" + uuid.uuid4().hex)
                try:
                    pending.write_text(json.dumps({"generation": name}))
                    os.replace(pending, pointer)
                finally:
                    pending.unlink(missing_ok=True)
            else:
                os.replace(staging, dest)
                root = dest
        finally:
            shutil.rmtree(staging, ignore_errors=True)
    print(f"[arda] installed IgBLAST -> {root}", file=sys.stderr)
    return root


def installed_version(dest: Path) -> str | None:
    """The IgBLAST release version in *dest*, or None. Provenance for a benchmark record."""
    marker = (installed_root(Path(dest)) or Path(dest)) / VERSION_FILE
    return marker.read_text().strip() if marker.exists() else None


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description="Download an NCBI IgBLAST release for arda.")
    ap.add_argument("--dest", required=True, help="Install location; the resolved tool root is printed.")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    fetch(Path(args.dest), force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
