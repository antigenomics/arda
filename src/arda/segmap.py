"""Segment assignment without a homology search.

The two-pass path searches a 924-target segment reference purely to learn, per read, its best V
allele and its best J allele with coordinates. That is a *structural* question about a fixed 236 kb
germline reference, not a homology search, and answering it structurally is dramatically cheaper.
Measured on 100,000 amplicon reads against the shipped reference, 8 threads:

===========================  =========  ===============================================
step                         wall       agreement with `mmseqs search`
===========================  =========  ===============================================
``mmseqs search``            2770 ms    —
:func:`segment_rows`           74 ms    V allele .9997, J allele .9998, C allele 1.0000
===========================  =========  ===============================================

37x, with 53,048 reads getting a hit against mmseqs' 53,121 (-0.14 %). The residual disagreements
are almost entirely allele-level inside one gene (``TRAV36/DV7*01`` vs ``*04``) — the degeneracy
:data:`arda.annotate.mapper._MAX_TIED_V` exists for.

⛔ **This nominates candidates; it does not decide anything.** The winner is still aligned against
the full V+pad+J scaffold by MMseqs2 and scored there, so the contract is that arda's AIRR output
must not move — not that these scores match mmseqs' bit scores. They are on a different scale by
design (see ``MATCH``/``MISMATCH`` in ``src/_segmap/segmap.cpp``).

Off by default. Enable with ``arda map --two-pass --fast-segments``.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

__all__ = ["available", "build", "segment_rows", "K", "MIN_SCORE", "SIDE_GROUP"]

logger = logging.getLogger(__name__)

try:  # pragma: no cover - import guard
    from . import _segmap
except ImportError:  # pragma: no cover
    _segmap = None  # type: ignore[assignment]

#: Seed length. **Matches `mapper._KMER["nt"]`, the `-k` arda passes MMseqs2**, and that is the
#: point: the segment pass has to be as sensitive as the search it replaces, and seed length is
#: what sets sensitivity to mismatches. 16 was inherited from `prefilter`, where it is calibrated
#: for REJECTION (k=12 passes 62-65 % of reads there, which is useless as a filter) -- the opposite
#: problem. Measured on 100,000 amplicon reads, reads with a hit:
#:
#:   k=16  53,048     k=14  53,087     k=13  53,108     k=12  53,121     mmseqs  53,121
#:
#: k=12 reproduces mmseqs exactly. It costs 129 ms against 68 ms at k=16 -- still ~21x the 2,770 ms
#: search, so sensitivity is the right thing to spend it on.
K = 12

#: Significance floor, calibrated against mmseqs rather than chosen — see the C++ source. mmseqs
#: applies ``-e 1e-3``; this scheme has no e-value, so without a floor every seeded diagonal with a
#: positive extension is reported and 43,010 reads pick up a constant-region hit against mmseqs' 473.
MIN_SCORE = 40

#: Target prefix -> group index handed to the native mapper, which returns the best hit per
#: (read, group).
#:
#: ⛔ These are the same three sides as ``mapper._SEGMENT_SIDE`` and must stay that way: `C` is its
#: OWN group because a constant-region hit says what the isotype is and nothing about which J the
#: read carries, and `JC` (the pre-2.8.0 kind) stays J-side. ``test_segmap_wiring`` asserts the two
#: mappings agree rather than trusting this comment.
SIDE_GROUP = {"V": 0, "J": 1, "JC": 1, "C": 2}


def available() -> bool:
    """Is the native extension importable?"""
    return _segmap is not None


def _read_fasta(path: Path) -> tuple[list[str], list[str]]:
    names: list[str] = []
    seqs: list[str] = []
    cur: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if names:
                    seqs.append("".join(cur))
                cur = []
                names.append(line[1:].strip().split()[0])
            else:
                cur.append(line.strip())
    if names:
        seqs.append("".join(cur))
    return names, seqs


@lru_cache(maxsize=4)
def build(segments_fasta: Path, k: int = K):
    """Index ``segments.fasta``. Cached: the reference never changes within a run.

    Returns:
        ``(mapper, names)`` — ``names[i]`` is the target the mapper's ``target_index`` refers to.

    Raises:
        RuntimeError: if the native extension is not built.
    """
    if _segmap is None:
        raise RuntimeError("arda._segmap is not built; --fast-segments is unavailable")
    names, seqs = _read_fasta(Path(segments_fasta))
    unknown = {n.split("|", 1)[0] for n in names} - set(SIDE_GROUP)
    if unknown:
        raise ValueError(f"segments.fasta carries unknown target kinds: {sorted(unknown)}")
    groups = [SIDE_GROUP[n.split("|", 1)[0]] for n in names]
    mapper = _segmap.SegmentMapper(seqs, groups, k)
    logger.debug("segmap: %d targets, %d codes, %d postings",
                 len(names), mapper.size, mapper.postings)
    return mapper, names


#: Largest diagonal shift still read as one indel rather than a repeat or a stray seed, when indel
#: detection is enabled. IMGT V genes are ~300 nt and SHM indels are typically codon-sized, so this
#: is generous; the bound exists to stop two unrelated seeds on one target from routing a read to a
#: gapped realignment that would tell us nothing.
MAX_INDEL_NT = 30


def segment_rows(segments_fasta: Path, reads: dict[str, str], *, max_tied: int,
                 threads: int = 1, min_score: int = MIN_SCORE,
                 max_indel: int = 0) -> list[dict]:
    """Best V / best J / best C per read, in the shape ``mapper._segment_rows`` returns.

    Args:
        reads: ``{read_id: sequence}``. Insertion order is preserved and is the only thing tying a
            native row back to a read id.
        max_tied: exactly-tied targets to keep per side, i.e. ``_MAX_TIED``.
        max_indel: if > 0, flag rows whose target carries two well-supported diagonals at most this
            many nt apart — the signature of an indel, which a single ungapped extension cannot
            score. 0 disables the check and every ``split`` is 0.

    Returns:
        Row dicts with ``query``, ``target``, ``bits``, ``qstart``, ``qend``, ``tstart``, ``split``
        — the keys ``_segment_best_hits`` consumes, and **per read in descending score order**,
        which is the ordering its loop depends on.

    ⚠ ``bits`` here is an ungapped match/mismatch score, not an MMseqs2 bit score. Nothing
    downstream compares the two: the value is used to rank candidates within a read and to detect
    exact ties, both of which are internal to the segment pass.
    """
    mapper, names = build(Path(segments_fasta))
    ids = list(reads)
    hits = mapper.map([reads[i] for i in ids], max_tied, min_score, threads, max_indel)
    return [
        {"query": ids[qi], "target": names[ti], "bits": float(score),
         "qstart": qstart, "qend": qend, "tstart": tstart, "split": bool(split)}
        for qi, ti, score, qstart, qend, tstart, _rc, split in hits
    ]
