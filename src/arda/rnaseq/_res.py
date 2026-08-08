"""Wall time and peak RSS for one pipeline stage.

Only `map` used to report resources, which left the expensive stage unmeasured: mapping is
flat (~300-650 MB at any depth) while Stage 3 holds the clone set, and on a B-cell-rich tumour
(28,444 clonotypes from 105M reads) `correct` peaked at 2,071.7 MB. Anyone sizing a SLURM
`--mem` or Nextflow `memory` directive from the mapping number alone would be OOM-killed.

What the numbers mean, exactly -- because `resource.getrusage` offers no way to be more
precise, and a vaguer definition here would be a lie rather than a simplification:

``peak_rss_mb``
    The **whole process** (plus reaped children) high-water mark **as of the end of this
    stage**. Monotone across stages: it can only rise. `getrusage` reports high-water marks
    only -- `RUSAGE_SELF` since process start, `RUSAGE_CHILDREN` cumulative over children --
    and there is no per-stage reset, so a stage cannot be charged its own peak in isolation
    when all three run in one process (`arda rnaseq run`).
``rss_gain_mb``
    How much *this* stage raised that mark; 0 if it stayed under an earlier stage's peak.

The monotone number is the one an operator actually needs: it is what the process required at
that point, which is what a memory directive has to cover. For per-stage attribution, run the
stage in its own process (`arda rnaseq map|assemble|correct` separately) -- then
``peak_rss_mb`` is that stage alone.
"""

from __future__ import annotations

import resource
import sys
import time

__all__ = ["peak_rss_mb", "Stage"]


def peak_rss_mb() -> float:
    """Peak RSS of this process AND its children, in MB.

    ``RUSAGE_SELF`` alone is wrong and was: 92 % of a `map` run's wall time is spent inside the
    `mmseqs` **subprocess**, whose nucleotide prefilter allocates a ``4**k`` index table that
    dominates the footprint. Reporting only the Python process understated peak RSS by roughly
    an order of magnitude. ``ru_maxrss`` is bytes on macOS, KB on Linux.
    """
    scale = 1024 * 1024 if sys.platform == "darwin" else 1024
    return max(resource.getrusage(who).ru_maxrss
               for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)) / scale


class Stage:
    """Times a stage and records its resource footprint onto a report dataclass.

    Not a context manager: `correct_airr` and `assemble_contigs` write their JSON report
    *inside* the function, so measurement has to close before that write, not after the block.

        stage = Stage()
        ...
        stage.finish(report)          # sets wall_seconds / peak_rss_mb / rss_gain_mb
        if report_path: ...
    """

    def __init__(self) -> None:
        self._t0 = time.perf_counter()
        self._rss0 = peak_rss_mb()

    @property
    def wall_seconds(self) -> float:
        return time.perf_counter() - self._t0

    def finish(self, report=None):
        """Stamp the elapsed time and RSS onto *report* (if it has those fields)."""
        wall = self.wall_seconds
        peak = peak_rss_mb()
        if report is not None:
            for name, value in (("wall_seconds", round(wall, 3)),
                                ("peak_rss_mb", round(peak, 1)),
                                ("rss_gain_mb", round(max(0.0, peak - self._rss0), 1))):
                if hasattr(report, name):
                    setattr(report, name, value)
        return report
