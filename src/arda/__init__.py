"""arda — Antigen Receptor Domain Annotation.

Fast FR/CDR region annotation for TCR/BCR nucleotide and amino acid sequences,
via offline IgBLAST-built references mapped at runtime with MMseqs2.
"""

from __future__ import annotations

__version__ = "2.38.1"

from .adapter import annotate_sequences  # noqa: E402


def __getattr__(name):
    if name == "alignment_candidates":
        from .annotate.candidates import alignment_candidates

        return alignment_candidates
    if name == "reference_index_identity":
        from .annotate.mapper import reference_index_identity

        return reference_index_identity
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["alignment_candidates", "annotate_sequences", "reference_index_identity", "__version__"]
