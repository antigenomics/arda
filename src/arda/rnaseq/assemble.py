"""Stage 3 — contig assembly (DEFERRED — interface stub only).

Reconstruct full-length V(D)J contigs from the candidate reads that Stage 1 filtered,
using the per-read anchors (V/J/junction keyed by ``sequence_id``) it already emits —
the role assembly-based extractors play with de-novo assembly. Not implemented yet; see the plan and
``ROADMAP.md``. The I/O contract is fixed here so callers/pipelines can wire it early.

# ponytail: interface only until Stage 1+2 + the benchmark validate; assembly is a much
# larger build (a de-novo assembler), explicitly out of scope for the filter-first goal.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["assemble_contigs"]


def assemble_contigs(
    reads_fasta: str | Path,
    airr_tsv: str | Path,
    output: str | Path,
    *,
    organism: str = "human",
) -> Path:
    """Assemble full-length V(D)J contigs from Stage-1 candidate reads (NOT IMPLEMENTED).

    Args:
        reads_fasta: candidate reads from ``arda rnaseq map --emit-reads``.
        airr_tsv: the Stage-1 mapped-reads AIRR TSV (read-id → V/J/junction anchors).
        output: contig FASTA / AIRR TSV to write.

    Raises:
        NotImplementedError: always — Stage 3 is deferred (see module docstring).
    """
    raise NotImplementedError(
        "arda rnaseq assemble (Stage 3, contig assembly) is not implemented yet — "
        "see ROADMAP.md. Use `map` (filter+annotate) and `correct` (CDR3 error correction).")
