#!/usr/bin/env python3
"""Download an NCBI IgBLAST release into a directory. Thin wrapper; the logic is packaged.

``setup.sh`` calls this before arda is necessarily installed, so it puts ``src/`` on the path
rather than importing an installed ``arda``. The logic itself lives in
:mod:`arda._igblast_fetch`, because the *runtime* needs it too: a plain ``pip install`` never
runs ``setup.sh``, and ``arda igblast`` has to work there. Keeping a second copy here is what
let the two drift apart in the first place.

    python scripts/fetch_igblast.py --dest bin [--force]
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from arda._igblast_fetch import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
