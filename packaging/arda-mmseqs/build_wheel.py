#!/usr/bin/env python
# Build one platform-tagged `arda-mmseqs` wheel.
#
# Run from this directory, in any venv with `build` installed:
#     python build_wheel.py                       # this platform
#     python build_wheel.py --plat manylinux_2_17_x86_64 --asset mmseqs-linux-avx2.tar.gz
#
# Reuses arda's own platform->asset table (`arda._mmseqs_fetch`) so the wheel and the
# auto-fetch path can never disagree about which build is the right one.
#
# 2026-08-04
"""Build one platform-tagged ``arda-mmseqs`` wheel."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIN = HERE / "src" / "arda_mmseqs" / "bin"
# Repo root: packaging/arda-mmseqs -> repo. Lets this run without arda installed.
sys.path.insert(0, str(HERE.parents[1] / "src"))

from arda._mmseqs_fetch import default_asset, fetch  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--plat", default=sysconfig.get_platform().replace("-", "_").replace(".", "_"),
                    help="wheel platform tag (e.g. manylinux_2_17_x86_64, macosx_11_0_arm64)")
    ap.add_argument("--asset", default=None, help="MMseqs2 release asset (default: this platform)")
    ap.add_argument("--outdir", default=str(HERE / "dist"))
    args = ap.parse_args()

    # `fetch` selects via default_asset(), which honours $ARDA_MMSEQS_ASSET — that env var
    # is the cross-build knob, so there is nothing to plumb through.
    if args.asset:
        os.environ["ARDA_MMSEQS_ASSET"] = args.asset
    print(f"[arda-mmseqs] asset={default_asset()} plat={args.plat}")

    if BIN.exists():
        shutil.rmtree(BIN)
    BIN.mkdir(parents=True)

    fetch(BIN)  # downloads, extracts, chmods; lands at BIN/mmseqs
    print(f"[arda-mmseqs] bundled {(BIN / 'mmseqs').stat().st_size / 1e6:.1f} MB")

    # Ship the upstream licence with the binary — MIT permits redistribution WITH the notice.
    lic = HERE / "LICENSE.mmseqs2"
    if not lic.exists():
        print("[arda-mmseqs] WARNING: LICENSE.mmseqs2 missing; MIT requires shipping it", file=sys.stderr)

    subprocess.run([sys.executable, "-m", "build", "--wheel", "--outdir", args.outdir],
                   cwd=HERE, check=True)
    # Retag to py3-none-<plat>. `has_ext_modules()` is a lie we tell setuptools to stop it
    # producing a `py3-none-any` wheel (there is no extension module, just a binary in
    # package data) -- but it also makes setuptools stamp the building interpreter's ABI,
    # which would mean one wheel per Python version for no reason. There is no C extension,
    # so the wheel is ABI-independent: one per platform is enough.
    for whl in Path(args.outdir).glob("arda_mmseqs-*.whl"):
        if not whl.name.endswith(f"py3-none-{args.plat}.whl"):
            subprocess.run([sys.executable, "-m", "wheel", "tags",
                            "--python-tag", "py3", "--abi-tag", "none",
                            "--platform-tag", args.plat, "--remove", str(whl)], check=True)
    print(f"[arda-mmseqs] wheels in {args.outdir}")
    return 0



if __name__ == "__main__":
    raise SystemExit(main())
