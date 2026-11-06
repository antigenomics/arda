"""Download a static MMseqs2 binary into ``bin/`` (packaged, stdlib only).

Used by :func:`arda.mmseqs.mmseqs_binary` for transparent auto-fetch and by the
``scripts/fetch_mmseqs.py`` CLI wrapper. Kept in the package (not just scripts/)
so auto-fetch works for wheel installs, not only editable checkouts.
"""

from __future__ import annotations

import os
import hashlib
import subprocess
import zipfile
import http.client
import platform
import shutil
import stat
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

from ._locking import build_lock

RELEASE_URL = "https://github.com/soedinglab/MMseqs2/releases/latest/download/{asset}"

# GitHub rejects the default Python-urllib User-Agent on some networks; mimic a
# real client so the release redirect chain resolves.
_UA = "Mozilla/5.0 (arda mmseqs-fetch)"


def _download(url: str, dest: Path, *, retries: int = 3, max_seconds: float = 120) -> None:
    """Bound network waits, retry transient errors and verify the received byte count."""
    last: Exception | None = None
    deadline = time.monotonic() + max_seconds
    for attempt in range(retries):
        try:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f"Download exceeded {max_seconds:g}s")
            req = urllib.request.Request(url, headers={"User-Agent": _UA})
            with urllib.request.urlopen(req, timeout=min(30, remaining)) as resp, open(dest, "wb") as fh:
                count = 0
                read = getattr(resp, "read1", resp.read)
                while True:
                    if time.monotonic() >= deadline:
                        raise TimeoutError(f"Download exceeded {max_seconds:g}s")
                    block = read(1024 * 1024)
                    if not block:
                        break
                    fh.write(block)
                    count += len(block)
                length = resp.headers.get("Content-Length")
                if length is not None and count != int(length):
                    raise urllib.error.URLError(f"Incomplete download: {count}/{length} bytes")
            return
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException) as exc:
            last = exc
            if (isinstance(exc, urllib.error.HTTPError) and
                    400 <= exc.code < 500 and exc.code not in (408, 429)):
                break
            if attempt < retries - 1 and time.monotonic() < deadline:
                print(f"[arda] retrying download ({attempt + 2}/{retries}): {exc}", file=sys.stderr)
                time.sleep(min(2 * (attempt + 1), max(0, deadline - time.monotonic())))
    raise RuntimeError(f"failed to download {url}: {last}")


def default_asset() -> str:
    """Pick the MMseqs2 release asset for this platform.

    Honors ``$ARDA_MMSEQS_ASSET``. Linux x86-64 defaults to the AVX2 build
    (universal on ~2013+ hardware); set ``$ARDA_MMSEQS_ASSET`` to
    ``mmseqs-linux-sse41.tar.gz`` for older CPUs.
    """
    override = os.environ.get("ARDA_MMSEQS_ASSET")
    if override:
        return override
    system = platform.system()
    machine = platform.machine().lower()
    if system == "Darwin":
        return "mmseqs-osx-universal.tar.gz"
    if system == "Windows" and machine in ("amd64", "x86_64"):
        return "mmseqs-win64.zip"
    if system == "Linux":
        if machine in ("aarch64", "arm64"):
            return "mmseqs-linux-arm64.tar.gz"
        return "mmseqs-linux-avx2.tar.gz"
    raise RuntimeError(
        f"Unsupported platform for MMseqs2 auto-fetch: {system}/{machine}. "
        "Install "
        "a compatible MMseqs2 and set $ARDA_MMSEQS, or select a supported tar asset with "
        "$ARDA_MMSEQS_ASSET."
    )


def installed_binary(bin_dir: Path) -> Path:
    """Resolve an atomically published Windows runtime, or the Unix executable."""
    marker = Path(bin_dir) / "mmseqs-windows" / "CURRENT"
    if platform.system() == "Windows" and marker.is_file():
        generation = marker.read_text().strip()
        if len(generation) == 64 and all(c in "0123456789abcdef" for c in generation):
            return marker.parent / generation / "bin" / "mmseqs.exe"
    return Path(bin_dir) / ("mmseqs.exe" if platform.system() == "Windows" else "mmseqs")


def _install_windows(archive: Path, bin_dir: Path, tmp: Path) -> Path:
    """Stage the official DLL/BusyBox runtime; publish a complete immutable generation."""
    generation = hashlib.sha256(archive.read_bytes()).hexdigest()
    with zipfile.ZipFile(archive) as zf:
        for member in zf.infolist():
            name = member.filename.replace("\\", "/")
            path = (tmp / name).resolve()
            if (not path.is_relative_to(tmp.resolve()) or ":" in name or
                    stat.S_ISLNK(member.external_attr >> 16)):
                raise RuntimeError(f"Unsafe MMseqs2 zip member: {member.filename}")
        zf.extractall(tmp)
    source = tmp / "mmseqs"
    runtime = source / "bin"
    if not all((runtime / name).is_file() for name in ("mmseqs.exe", "busybox.exe", "cygwin1.dll")):
        raise RuntimeError("Unexpected MMseqs2 Windows archive layout")
    # Hard links need no Administrator/developer-mode permission, unlike upstream's -s.
    subprocess.run([str(runtime / "busybox.exe"), "--install", str(runtime)],
                   check=True, capture_output=True, timeout=60)
    subprocess.run([str(runtime / "mmseqs.exe"), "version"],
                   check=True, capture_output=True, timeout=60)
    root = bin_dir / "mmseqs-windows"
    root.mkdir(exist_ok=True)
    target = root / generation
    if not target.exists():
        os.replace(source, target)
    marker = tmp / "CURRENT"
    marker.write_text(generation)
    os.replace(marker, root / "CURRENT")
    return target / "bin" / "mmseqs.exe"


def fetch(bin_dir: Path, *, force: bool = False) -> Path:
    """Download + install the mmseqs binary into *bin_dir*; return its path.

    Concurrency-safe, because arda runs concurrently against one cache by design -- a Nextflow
    process per sample, a SLURM array task per shard -- and on a cold cache every one of them
    calls this at the same moment.

    It used to `shutil.copy2` straight onto ``bin/mmseqs``, which is this repo's signature bug in
    its purest form. Two failures, both silent:

    * copy2 writes **in place**, so two processes doing it at once interleave their bytes and the
      result is a corrupt binary that nobody reported an error about;
    * ``dest.exists()`` -- the gate that decides the work is already done -- becomes true on the
      **first byte** copy2 writes. A third process therefore sees an mmseqs, returns it, and
      executes a truncated file. (The chmod came after the copy too, so there was a second window
      where the file existed but was not executable.)

    Now: one downloader under the build lock, the binary made executable while still under a
    private name, and a single ``os.replace`` to publish it. A rename is atomic within a
    filesystem, so ``dest`` only ever exists as a complete, executable binary -- which is what
    makes ``dest.exists()`` a legitimate gate rather than a bug.
    """
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    dest = installed_binary(bin_dir)
    if dest.exists() and not force:
        return dest

    # Under --force the work is ours by definition; otherwise whoever we queued behind may have
    # finished it while we waited, and there is then nothing left to do.
    done = (lambda: False) if force else lambda: installed_binary(bin_dir).exists()
    with build_lock(bin_dir / ".mmseqs.lock", done=done) as ours:
        if not ours:
            return installed_binary(bin_dir)
        asset = default_asset()
        url = RELEASE_URL.format(asset=asset)
        print(f"[arda] fetching mmseqs: {url}", file=sys.stderr)
        # Scratch next to the destination, NOT in the system temp dir: on a cluster /tmp is
        # routinely a small node-local filesystem (aldan3's is 2.0 GB with 29 MB free) and this
        # archive plus its extraction does not fit. The IgBLAST fetch died there first, with
        # ENOSPC. It also puts the staged binary on the same filesystem as `dest`, which is what
        # makes the `os.replace` below a rename rather than a copy.
        with tempfile.TemporaryDirectory(prefix=".arda_mmseqs_", dir=bin_dir) as td:
            tmp = Path(td)
            tarball = tmp / asset
            _download(url, tarball)
            if asset.endswith(".zip"):
                return _install_windows(tarball, bin_dir, tmp)
            with tarfile.open(tarball) as tf:
                tmp_resolved = tmp.resolve()
                for member in tf.getmembers():
                    if member.issym() or member.islnk():
                        raise RuntimeError(
                            f"Refusing to extract link from MMseqs2 archive: {member.name}"
                        )
                    member_path = (tmp / member.name).resolve()
                    if not str(member_path).startswith(str(tmp_resolved) + os.sep):
                        raise RuntimeError(
                            f"Refusing to extract path outside temp dir: {member.name}"
                        )
                tf.extractall(tmp, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
            extracted = next(tmp.rglob("bin/mmseqs"), None)  # static archives: mmseqs/bin/mmseqs
            if extracted is None:
                raise RuntimeError(f"Unexpected MMseqs2 archive layout in {asset}")
            staged = tmp / "mmseqs.staged"
            shutil.copy2(extracted, staged)
            # Executable BEFORE it is published, never after: a rename cannot be observed
            # half-done, but a chmod that follows the publish can.
            staged.chmod(staged.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
            os.replace(staged, dest)
    print(f"[arda] installed mmseqs -> {dest}", file=sys.stderr)
    return dest
