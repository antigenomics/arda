"""Crash-releasing cross-process and cross-thread coordination for cache builds."""
from __future__ import annotations

import errno
import os
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

LOCK_TIMEOUT_S = 900


@contextmanager
def build_lock(lock: Path, *, done: Callable[[], bool], timeout: float | None = None) -> Iterator[bool]:
    """Yield whether this caller must build; the kernel releases ownership on exit or death.

    A persistent guard file replaces the old mkdir lock. Legacy empty lock directories
    do not block recovery. Never unlink the guard: waiters may already hold its inode.
    """
    if done():
        yield False
        return
    lock.parent.mkdir(parents=True, exist_ok=True)
    guard = lock.with_name(lock.name + ".guard")
    timeout = LOCK_TIMEOUT_S if timeout is None else timeout
    deadline = time.monotonic() + timeout
    with guard.open("a+b") as handle:
        if os.name == "nt":
            import msvcrt
            if handle.seek(0, 2) == 0:
                handle.write(b"\0")
                handle.flush()

            def acquire():
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)

            def release():
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            def acquire():
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)

            def release():
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

        while True:
            try:
                acquire()
                break
            except OSError as exc:
                if exc.errno not in (errno.EACCES, errno.EAGAIN, errno.EDEADLK):
                    raise
                if done():
                    yield False
                    return
                if time.monotonic() >= deadline:
                    raise TimeoutError(f"Another active arda build still owns {guard} after "
                                       f"{timeout:g}s") from exc
                time.sleep(0.05)
        try:
            yield not done()
        finally:
            release()
