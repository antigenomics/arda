"""CDR3 repair by k-mer scan against the **flank-extended** germline.

A port of VDJdb's ``Cdr3Fixer``/``KmerScanner``, which :mod:`arda.cdr3fix` replaced and does not
reproduce on one class of input. The difference is the reference, not the algorithm.

:mod:`arda.cdr3fix` aligns against ``Anchor.templated_aa``, which runs Cys104 through [FW]118
**inclusive and no further**. A junction carrying framework past an anchor therefore has nothing to
align that framework *to*, so the boundary can only be inferred from the junction side, and arda
reports the finding without applying it: ``YFCASSQSPGGVAFFGQG`` on ``TRBV14``/``TRBJ1-1`` comes back
``V ins@0 extra 'YF'`` and ``J ins@13 extra 'GQG'``, both "reported, not repaired", and the shipped
junction keeps the framework.

Here the segment sequence is sliced from the full allele so that it *starts* (V) or *ends* (J)
outside the junction, exactly as VDJdb did:

* **V**: ``allele_nt[anchor_nt - 3:]`` translated, so one codon of FR3 precedes Cys104;
* **J**: ``allele_nt[:anchor_nt + 4]`` translated in the J frame, so the [FW]GXG motif follows 118.

With flanking context present the repair is positional and needs no scoring. One scan returns the
longest common substring and where it sits in each sequence, and the four outcomes follow from the
two offsets:

=========================  =====================  ===============================================
``start_in_segment``       ``start_in_cdr3``      outcome
=========================  =====================  ===============================================
0                          0                      ``NoFixNeeded``
0                          > 0                    ``FixTrim`` -- drop the flank before the anchor
> 0                        0                      ``FixAdd`` -- restore the germline the cut removed
> 0                        <= ``MAX_REPLACE``     ``FixReplace`` -- substitute at the anchor
> 0                        > ``MAX_REPLACE``      ``FailedReplace``
=========================  =====================  ===============================================

The J end runs the same code on reversed strings, which is why there is one scanner rather than two.

``NoFixNeeded`` is an **alignment** verdict and not an anchor-residue check: it requires a hit of at
least :data:`MIN_HIT` residues at offset zero in both sequences. That distinction is the reason this
module exists alongside the aligner - ``CASSQQQQQQQQQF`` on ``TRBJ1-1`` (``NTEAFF``) shares only its
terminal Phe, which is one residue of coincidence and no alignment at all.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

#: Shortest k-mer that counts as a hit. VDJdb's ``min_hit_size``. Two residues at the anchor is the
#: floor at which a match is evidence rather than coincidence.
MIN_HIT = 2

#: How many residues may be substituted at the anchor. VDJdb's ``max_replace_size``. Beyond this the
#: record is likelier assigned to the wrong allele than to hold that many typos, so it declines.
MAX_REPLACE = 1

#: One codon of flank on the V side and four nucleotides on the J side, which is what puts the
#: framework into the scanned sequence. VDJdb's slice offsets, kept verbatim.
V_FLANK_NT, J_FLANK_NT = 3, 4

_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str, *, from_end: bool = False) -> str:
    """Translate in one frame. ``from_end`` drops the leading remainder so the frame ends flush.

    A J germline is in frame with its 3' end, not its 5' one, which is what ``from_end`` is for.
    An unknown codon carrying an ``N`` becomes ``X`` and anything else ``?``, as VDJdb's
    ``translate_linear`` did: a reference defect should be visible, not silently dropped.
    """
    if from_end:
        seq = seq[len(seq) % 3:]
    out = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        out.append(_CODONS.get(codon, "X" if "N" in codon else "?"))
    return "".join(out)


@dataclass(frozen=True)
class Hit:
    """The longest common substring: where it starts in each sequence, and how long it is."""

    start_in_segment: int
    start_in_cdr3: int
    match_size: int


class KmerScanner:
    """Every substring of one germline segment, indexed by its first position.

    Built once per segment and reused across queries, which is what makes the scan cheap: the index
    is O(len^2) substrings of a sequence tens of residues long, and a query is a dictionary lookup
    per substring.
    """

    __slots__ = ("kmers", "min_hit")

    def __init__(self, seq: str, min_hit: int = MIN_HIT) -> None:
        self.min_hit = min_hit
        self.kmers: dict[str, int] = {}
        for size in range(min_hit, len(seq) + 1):
            for start in range(len(seq) - size + 1):
                self.kmers.setdefault(seq[start:start + size], start)

    def scan(self, seq: str) -> Hit | None:
        """The longest substring of ``seq`` that occurs in the segment, or ``None``.

        Ties keep the first position found, scanning shortest-to-longest and left-to-right, which is
        the order VDJdb's scanner used.
        """
        best: Hit | None = None
        for size in range(self.min_hit, len(seq) + 1):
            for start in range(len(seq) - size + 1):
                found = self.kmers.get(seq[start:start + size])
                if found is not None and (best is None or size > best.match_size):
                    best = Hit(found, start, size)
        return best
