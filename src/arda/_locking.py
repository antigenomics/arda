"""A cross-process build lock, for the two places arda races with itself.

arda is routinely run concurrently against the **same** cache: the Nextflow module launches one
process per sample and a SLURM array one per task. On first use in a fresh environment every one of
them finds no reference and no mmseqs index, and starts building the same thing into the same path.

Both races have the same shape and the same silent failure: an artifact whose mere *presence* is the
"already built" gate becomes visible before it is complete, so every other process happily uses a
half-written one and reports success over nothing.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

LOCK_TIMEOUT_S = 900
"""These builds take seconds. Fifteen minutes means a dead lock, not a slow one."""


@contextmanager
def build_lock(lock: Path, *, done: Callable[[], bool]) -> Iterator[bool]:
    """Serialize a build at ``lock``. Yields ``True`` if the work is ours, ``False`` if it is done.

    ``mkdir`` is the lock: it is atomic on POSIX and Windows, needs no dependency, and leaves nothing
    to clean up but an empty directory.

    ``done()`` is checked while waiting *and* again after acquiring, because the common case is not
    contention over the work — it is that the process we queued behind already finished it, and there
    is nothing left for us to do.
    """
    deadline = time.monotonic() + LOCK_TIMEOUT_S
    lock.parent.mkdir(parents=True, exist_ok=True)
    while True:
        try:
            lock.mkdir()
            break
        except FileExistsError:
            if done():  # the holder finished while we waited
                yield False
                return
            if time.monotonic() > deadline:
                raise TimeoutError(
                    f"waited {LOCK_TIMEOUT_S}s for another arda process to finish {lock.parent}. "
                    f"If no other arda run is active, remove {lock} and retry."
                )
            time.sleep(0.5)
    try:
        yield not done()  # won the lock, but someone may have finished just before we took it
    finally:
        lock.rmdir()
