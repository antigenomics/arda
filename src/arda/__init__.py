"""arda — Antigen Receptor Domain Annotation.

Fast FR/CDR region annotation for TCR/BCR nucleotide and amino acid sequences,
via offline IgBLAST-built references mapped at runtime with MMseqs2.
"""

from __future__ import annotations

__version__ = "2.0.1"

from .adapter import annotate_sequences  # noqa: E402

__all__ = ["annotate_sequences", "__version__"]
