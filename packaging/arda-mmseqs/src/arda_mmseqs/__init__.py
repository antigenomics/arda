"""A static MMseqs2 binary shipped as a wheel, so ``arda`` needs nothing installed.

``pip install arda-mapper`` fetches mmseqs over the network on first use.
``pip install 'arda-mapper[mmseqs]'`` pulls this package instead and the binary is simply
there — no conda channel, no PATH surgery, no first-run download.

Two distributions rather than two variants of one, because pip cannot be asked to *prefer*
a wheel: build tags are ranked, not chosen. This is the same shape ``cmake``, ``ninja``,
``ruff`` and ``patchelf`` use on PyPI.

MMseqs2 is MIT-licensed (Steinegger & Söding); ``LICENSE.mmseqs2`` travels with the wheel.
"""

from __future__ import annotations

import os
import stat
import subprocess
from functools import lru_cache
from pathlib import Path

__all__ = ["mmseqs_version", "binary"]

_EXE = "mmseqs.exe" if os.name == "nt" else "mmseqs"


@lru_cache(maxsize=1)
def mmseqs_version() -> str | None:
    """What the bundled binary reports for ``mmseqs version``.

    Asked, not asserted. A hard-coded constant here would be a lie on at least one platform:
    the same MMseqs2 build spells itself differently per asset -- the osx-universal release
    prints a full commit hash (``8cc5ce367b...``) where bioconda prints ``18.8cc5c``. Running
    the binary is the only answer that is right everywhere.
    """
    try:
        p = subprocess.run([str(binary()), "version"], capture_output=True, text=True, timeout=60)
    except Exception:  # noqa: BLE001 — a missing/!exec binary is reported by binary() itself
        return None
    return p.stdout.strip() or None if p.returncode == 0 else None


def binary() -> Path:
    """Absolute path to the bundled ``mmseqs`` executable.

    Ensures the execute bit: wheels are zip archives and pip does not reliably preserve
    the mode of files outside ``scripts/``.

    Raises:
        FileNotFoundError: if the wheel was built without a binary for this platform.
    """
    path = Path(__file__).resolve().parent / "bin" / _EXE
    if not path.exists():
        raise FileNotFoundError(
            f"arda-mmseqs carries no binary at {path}. This wheel was built for a "
            f"different platform; reinstall, or use plain `arda-mapper` and let it fetch one."
        )
    mode = path.stat().st_mode
    if not mode & stat.S_IXUSR:
        path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path
