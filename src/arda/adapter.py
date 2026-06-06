"""Library-facing API.

A stable, import-friendly surface for embedding arda in other Python tools.
The heavy lifting lands in Phase 2 (``arda.annotate``); this module keeps the
public signature stable.
"""

from __future__ import annotations

from typing import Iterable, Literal

SeqType = Literal["nt", "aa"]


def annotate_sequences(
    sequences: Iterable[str] | Iterable[tuple[str, str]],
    seqtype: SeqType = "nt",
    organism: str = "human",
):
    """Annotate FR/CDR regions for a batch of sequences.

    Args:
        sequences: Either raw sequence strings or ``(id, sequence)`` pairs.
        seqtype: ``"nt"`` for nucleotide input, ``"aa"`` for amino acid.
        organism: One of the supported organisms (human, mouse, rat, rabbit,
            rhesus_monkey).

    Returns:
        A list of AIRR-style annotation record dicts (one per input sequence).
    """
    from .annotate.mapper import annotate_records

    pairs: list[tuple[str, str]] = []
    for i, item in enumerate(sequences):
        if isinstance(item, str):
            pairs.append((f"seq{i}", item))
        else:
            pairs.append((str(item[0]), str(item[1])))
    return annotate_records(pairs, organism=organism, seqtype=seqtype)
