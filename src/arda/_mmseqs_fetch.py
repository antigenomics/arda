"""Download a static MMseqs2 binary into ``bin/`` (packaged, stdlib only).

Used by :func:`arda.mmseqs.mmseqs_binary` for transparent auto-fetch and by the
``scripts/fetch_mmseqs.py`` CLI wrapper. Kept in the package (not just scripts/)
so auto-fetch works for wheel installs, not only editable checkouts.
"""

from __future__ import annotations

import os
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

RELEASE_URL = "https://github.com/soedinglab/MMseqs2/releases/latest/download/{asset}"

# GitHub rejects the default Python-urllib User-Agent on some networks; mimic a
# real client so the release redirect chain resolves.
_UA = "Mozilla/5.0 (arda mmseqs-fetch)"


def _download(url: str, dest: Path, *, retries: int = 3) -> None:
    """Download *url* to *dest* with a browser UA and simple retry/backoff."""
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": _UA})
            with urllib.request.urlopen(req, timeout=600) as resp, open(dest, "wb") as fh:
                shutil.copyfileobj(resp, fh)
            return
        except (urllib.error.URLError, TimeoutError) as exc:  # incl. HTTPError
            last = exc
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
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
    if system == "Linux":
        if machine in ("aarch64", "arm64"):
            return "mmseqs-linux-arm64.tar.gz"
        return "mmseqs-linux-avx2.tar.gz"
    raise RuntimeError(
        f"Unsupported platform for MMseqs2 auto-fetch: {system}/{machine}. "
        "Install mmseqs2 manually (conda install -c bioconda mmseqs2) and set "
        "$ARDA_MMSEQS, or set $ARDA_MMSEQS_ASSET."
    )


def fetch(bin_dir: Path, *, force: bool = False) -> Path:
    """Download + install the mmseqs binary into *bin_dir*; return its path."""
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    dest = bin_dir / "mmseqs"
    if dest.exists() and not force:
        return dest

    asset = default_asset()
    url = RELEASE_URL.format(asset=asset)
    print(f"[arda] fetching mmseqs: {url}", file=sys.stderr)
    with tempfile.TemporaryDirectory(prefix="arda_mmseqs_") as td:
        tmp = Path(td)
        tarball = tmp / asset
        _download(url, tarball)
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
            tf.extractall(tmp)
        extracted = next(tmp.rglob("bin/mmseqs"), None)  # static archives: mmseqs/bin/mmseqs
        if extracted is None:
            raise RuntimeError(f"Unexpected MMseqs2 archive layout in {asset}")
        shutil.copy2(extracted, dest)
        dest.chmod(dest.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    print(f"[arda] installed mmseqs -> {dest}", file=sys.stderr)
    return dest
