#!/usr/bin/env python3
"""Download the latest IgBLAST release and lay it out under ``bin/``.

Standalone (stdlib only) so it can run before ``arda`` is installed. Resolves
the platform-appropriate tarball from the NCBI ``LATEST`` directory, extracts
it, and arranges ``bin/`` so that it contains the executables *and* the
``internal_data`` / ``optional_file`` trees (IgBLAST reads them via ``$IGDATA``).

Idempotent: a no-op if ``bin/igblastn`` already exists (use ``--force`` to redo).
"""

from __future__ import annotations

import argparse
import platform
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

LATEST_URL = "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/"


def platform_suffix() -> str:
    system = platform.system()
    if system == "Linux":
        return "x64-linux"
    if system == "Darwin":
        # NCBI ships an x64 macOS build; on Apple Silicon it runs via Rosetta 2.
        return "x64-macosx"
    raise SystemExit(f"Unsupported platform for IgBLAST: {system}")


def find_tarball(suffix: str) -> str:
    """Scrape the LATEST listing for ``ncbi-igblast-<ver>-<suffix>.tar.gz``."""
    with urllib.request.urlopen(LATEST_URL, timeout=60) as resp:
        html = resp.read().decode("utf-8", "replace")
    pat = re.compile(rf"ncbi-igblast-([0-9.]+)-{re.escape(suffix)}\.tar\.gz")
    names = sorted(set(pat.findall(html)))
    if not names:
        raise SystemExit(
            f"No IgBLAST tarball for suffix {suffix!r} found at {LATEST_URL}"
        )
    version = names[-1]
    return f"ncbi-igblast-{version}-{suffix}.tar.gz"


def download(url: str, dest: Path) -> None:
    print(f"[fetch_igblast] downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=600) as resp, open(dest, "wb") as fh:
        shutil.copyfileobj(resp, fh)


def lay_out(extracted_root: Path, bin_dir: Path) -> None:
    """Copy bin/* + internal_data + optional_file into ``bin_dir``."""
    bin_dir.mkdir(parents=True, exist_ok=True)
    src_bin = extracted_root / "bin"
    if not src_bin.is_dir():
        raise SystemExit(f"Unexpected IgBLAST layout: no bin/ under {extracted_root}")
    for item in src_bin.iterdir():
        shutil.copy2(item, bin_dir / item.name)
        (bin_dir / item.name).chmod(0o755)
    for tree in ("internal_data", "optional_file"):
        src = extracted_root / tree
        if src.is_dir():
            shutil.copytree(src, bin_dir / tree, dirs_exist_ok=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", required=True, help="Target bin/ directory.")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    bin_dir = Path(args.dest)
    if (bin_dir / "igblastn").exists() and not args.force:
        print("[fetch_igblast] bin/igblastn present — skipping (use --force).",
              file=sys.stderr)
        return 0

    suffix = platform_suffix()
    tarball = find_tarball(suffix)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        tar_path = tmp / tarball
        download(LATEST_URL + tarball, tar_path)
        with tarfile.open(tar_path) as tf:
            tf.extractall(tmp)  # noqa: S202 — trusted NCBI source
        roots = [p for p in tmp.iterdir() if p.is_dir() and p.name.startswith("ncbi-igblast")]
        if not roots:
            raise SystemExit("Could not find extracted ncbi-igblast directory.")
        lay_out(roots[0], bin_dir)
    print(f"[fetch_igblast] installed IgBLAST into {bin_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
