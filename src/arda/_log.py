"""Logging, progress and resource reporting.

One ``arda`` logger, configured once from the CLI callback, that every module already feeds:
``cdr3fix``, ``prefilter``, ``segmap``, ``refbuild.*`` all call ``logging.getLogger(__name__)``,
so they are children of ``arda`` and inherit whatever ``setup`` installs. Nothing else has to
know about verbosity.

⛔ **Progress goes to stderr, results go to stdout.** The stage lines used to be ``typer.echo``
on stdout, which means ``arda export-ref`` piped into a file interleaved a progress line with the
data. Everything informational is a log record now; only paths and the export payload stay on
stdout.

``--log-file`` is always DEBUG whatever the console level is, and its format carries a timestamp
and the process peak RSS. That is the artifact a cluster job leaves behind, and re-running a
10-hour bulk sample because the console was at the default level is not a thing anyone should
have to do.
"""

from __future__ import annotations

import logging
import os
import sys
import time
from pathlib import Path

__all__ = ["logger", "setup", "Throttle", "peak_rss_mb"]

logger = logging.getLogger("arda")


def peak_rss_mb() -> float:
    """Peak RSS of this process AND its children, in MB; 0.0 where ``resource`` is absent.

    ``RUSAGE_SELF`` alone is wrong and was: 92 % of a `map` run's wall time is spent inside the
    `mmseqs` **subprocess**, whose nucleotide prefilter allocates a ``4**k`` index table that
    dominates the footprint. Reporting only the Python process understated peak RSS by roughly an
    order of magnitude. ``ru_maxrss`` is bytes on macOS, KB on Linux.

    ⛔ Lives HERE, in the module with no arda imports, and is re-exported by
    :mod:`arda.rnaseq._res`. The other direction is an import cycle: ``arda.rnaseq.__init__``
    imports ``map``, which needs :class:`Throttle`.
    """
    try:
        import resource  # noqa: PLC0415 — POSIX-only; wheels ship for Windows too
    except ImportError:  # pragma: no cover
        return 0.0
    scale = 1024 * 1024 if sys.platform == "darwin" else 1024
    return max(resource.getrusage(who).ru_maxrss
               for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)) / scale

_CONSOLE = "[arda] %(message)s"
_CONSOLE_V = "[arda] %(levelname)s %(name)s: %(message)s"
_FILE = "%(asctime)s %(levelname)-7s %(name)-24s %(rss)7.1f MB  %(message)s"


class _Rss(logging.Filter):
    """Injects ``%(rss)s`` -- process+children peak RSS in MB -- into every record.

    A filter rather than a call at each log site: memory is the number a cluster operator sizes
    a job from, and it is only useful if EVERY line carries it, including the ones written by
    modules that know nothing about this.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        record.rss = peak_rss_mb()
        return True


class Throttle:
    """One-every-``seconds`` gate for a progress line inside a hot loop.

    A bulk `map` run flushes hundreds of chunks; logging each at INFO turns a 100 M-read sample
    into 500 lines of noise, and logging none leaves a multi-hour job with no sign of life. Time
    is the right axis, not chunk count -- chunk wall time varies ~50x with the receptor fraction.
    """

    def __init__(self, seconds: float = 30.0) -> None:
        self.seconds = seconds
        self._last = time.monotonic()

    def ready(self) -> bool:
        now = time.monotonic()
        if now - self._last < self.seconds:
            return False
        self._last = now
        return True


def setup(verbosity: int = 0, quiet: bool = False, log_file: str | Path | None = None) -> None:
    """Configure the ``arda`` logger. Idempotent -- handlers are replaced, never stacked.

    Args:
        verbosity: 0 = INFO (the stage/progress lines), 1+ = DEBUG with level and module names.
        quiet: WARNING and above only. Loses to ``--log-file``, which stays at DEBUG.
        log_file: also write DEBUG records here, with timestamps and peak RSS per line.
    """
    logger.handlers.clear()
    logger.propagate = False
    logger.setLevel(logging.DEBUG)  # the handlers decide; the logger must not gate them

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(logging.WARNING if quiet
                     else logging.INFO if verbosity <= 0 else logging.DEBUG)
    console.setFormatter(logging.Formatter(_CONSOLE if verbosity <= 0 else _CONSOLE_V))
    logger.addHandler(console)

    if log_file is not None:
        handler = logging.FileHandler(log_file, mode="w")
        handler.setLevel(logging.DEBUG)
        handler.addFilter(_Rss())
        handler.setFormatter(logging.Formatter(_FILE, datefmt="%Y-%m-%d %H:%M:%S"))
        logger.addHandler(handler)

    from . import __version__

    logger.debug("arda %s | python %s | %s | %d cores | pid %d",
                 __version__, sys.version.split()[0], sys.platform,
                 os.cpu_count() or 0, os.getpid())
    logger.debug("argv: %s", " ".join(sys.argv))
