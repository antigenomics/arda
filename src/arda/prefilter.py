"""Exact k-mer prefilter — drop reads that cannot align before MMseqs2 sees them.

On bulk RNA-seq, ``mmseqs search`` spends nearly all of its time proving that reads are *not*
receptor reads: 4 M reads of SRR10611239 take 48.9 s to find 947 hits (0.024 %). The fitted cost
model says why — ``wall ~ reads/46,353 + hits/350``, so the dominant term is set by the read count
and not by the answer. A read can only align to a V(D)J scaffold if it shares an exact k-mer with
one; testing that is a lookup, proving it is Smith-Waterman.

The index is built from :attr:`~arda.annotate.reference.Reference.target_fasta` — **the same FASTA
MMseqs2 searches**. That is deliberate: the design's largest single finding is that a prefilter
built over ``V+pad+J`` alone loses 16.29 % of real reads, **69.27 % of them J->C reads**, and that
indexing the constant region takes the loss to 0.53 % (OPTIMIZATION.md §3.3). Deriving the index
from a hand-listed set of segments is exactly how that hole would come back; deriving it from the
search target makes it structurally impossible for the two to disagree about what is indexable.

The filter is **off by default**. It trades a measured ~0.5 % of real reads for ~6x on bulk, and
arda's near-zero Stage-1 false-negative rate is the thing it is not allowed to trade silently.
"""
from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

logger = logging.getLogger(__name__)

#: Seed length. k=12 passes 62-65 % of reads (no filtering); k>=18 adds nothing over 16.
K = 16
#: Windows a read must share with the reference to survive. A real read sharing one exact 16-mer
#: usually shares many, so >=2 barely moves the pass rate (4.64 % -> 3.91 %) while FN climbs.
MIN_HITS = 1
#: Above this pass rate the filter costs more than it saves — amplicon libraries run 46-49 %
#: receptor, where MMseqs2 has to look at nearly every read anyway.
MAX_USEFUL_PASS_RATE = 0.30


def _read_fasta_seqs(path: Path) -> list[str]:
    """Sequences only — the prefilter indexes bases, and the ids are never used."""
    seqs: list[str] = []
    cur: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


@lru_cache(maxsize=2)
def build(target_fasta: Path, k: int = K):
    """Index every k-mer of ``target_fasta`` and its reverse complement.

    Cached per FASTA: the index is ~1-3 MB and takes well under a second to build, but a
    per-chunk rebuild would charge every chunk of a 4 M-read run for it.
    """
    from ._prefilter import Prefilter  # noqa: PLC0415 — optional native extension

    seqs = _read_fasta_seqs(Path(target_fasta))
    idx = Prefilter(seqs, k)
    logger.debug("prefilter: %d targets -> %d distinct %d-mers", len(seqs), idx.size, k)
    return idx


def available() -> bool:
    """Is the native extension importable? A source tree without a built ext is a normal state."""
    try:
        from . import _prefilter  # noqa: F401,PLC0415
    except ImportError:
        return False
    return True


def keep_mask(records: list[tuple[str, str]], target_fasta: Path, *,
              threads: int = 1, min_hits: int = MIN_HITS, k: int = K) -> list[int]:
    """1 for each ``(id, sequence)`` record worth searching, 0 for the rest."""
    idx = build(Path(target_fasta), k)
    return list(idx.mask([seq for _sid, seq in records], min_hits, threads))
