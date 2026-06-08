#!/usr/bin/env python3
"""CLI to download a static MMseqs2 binary into ``bin/`` (see arda._mmseqs_fetch).

Usage::

    python scripts/fetch_mmseqs.py [--dest bin] [--force]

The download logic lives in the packaged module :mod:`arda._mmseqs_fetch` so the
same code path serves both this CLI and the transparent auto-fetch in
:func:`arda.mmseqs.mmseqs_binary`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from arda._mmseqs_fetch import fetch
from arda.paths import bin_dir


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dest", type=Path, default=None,
                    help="bin directory (default: arda project bin/)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    fetch(args.dest if args.dest is not None else bin_dir(), force=args.force)
