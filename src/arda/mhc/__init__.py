"""MHC domain annotation — scaffold only (deferred).

The TCR/BCR pipeline annotates FR/CDR regions. The analogous MHC task is to
annotate the peptide-binding groove (β-sheet floor) and the two α-helices for
class I (α1/α2 on a single heavy chain + β2-microglobulin) and class II (α1 on
the α chain, β1 on the β chain).

This module pre-fetches MHC allele references (see ``scripts/fetch_mhc.py`` →
``database/mhc/<organism>/``) but does **not** yet map the groove/helix domains.
That work will derive domain coordinates from solved structures (PDB) and project
them onto query MHC sequences with the same alignment-transfer machinery used for
antigen receptors. See ``memory/`` for the plan.
"""

from __future__ import annotations

__all__ = ["annotate_mhc"]


def annotate_mhc(sequences, mhc_class: str = "I", organism: str = "human"):
    """Annotate MHC groove/helix domains (not yet implemented).

    Args:
        sequences: Iterable of MHC protein sequences (or ``(id, seq)`` pairs).
        mhc_class: ``"I"`` or ``"II"``.
        organism: Reference organism.

    Raises:
        NotImplementedError: Always — MHC domain mapping is deferred. The allele
            references are available under ``database/mhc/`` for development.
    """
    raise NotImplementedError(
        "MHC groove/helix annotation is not implemented yet. Allele references are "
        "pre-fetched under database/mhc/; domain mapping (from PDB-derived "
        "coordinates) is a future phase."
    )
