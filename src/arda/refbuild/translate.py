"""Native nucleotide translation and reading-frame utilities.

No BioPython. The hot functions (``translate``, ``detect_coding_frame``,
``reverse_complement``, ``back_translate``) are implemented in the C++ extension
``arda._markup`` and re-exported here; a pure-Python fallback keeps the module
importable if the extension is unavailable. These mirror mirpy's mirseq API so
mirpy can later ``import arda`` and reuse them.
"""

from __future__ import annotations

__all__ = [
    "CODON_TABLE",
    "translate",
    "detect_coding_frame",
    "reverse_complement",
    "back_translate",
    "aa_coords_from_nt",
]

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# Build the standard genetic code: codons in TCAG×TCAG×TCAG order map to _AA.
CODON_TABLE: dict[str, str] = {}
_i = 0
for _a in _BASES:
    for _b in _BASES:
        for _c in _BASES:
            CODON_TABLE[_a + _b + _c] = _AA[_i]
            _i += 1


# Human (Kazusa) most-frequent codon per amino acid — used for back-translation.
_HUMAN_CODON = {
    "A": "GCC", "R": "AGG", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG",
    "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC", "T": "ACC", "W": "TGG",
    "Y": "TAC", "V": "GTG",
}

try:  # Fast path: C++ extension.
    from .. import _markup as _ext

    def translate(nt: str, frame: int = 0) -> str:
        """Translate a nucleotide string from ``frame`` (0/1/2)."""
        return _ext.translate(nt, frame)

    def detect_coding_frame(nt: str) -> int:
        """Return the reading frame (0/1/2) with the fewest stop codons."""
        return _ext.detect_coding_frame(nt)

    def reverse_complement(nt: str) -> str:
        """Reverse-complement a nucleotide string (non-ACGT -> ``N``)."""
        return _ext.reverse_complement(nt)

    def back_translate(aa: str, unknown: str = "NNN") -> str:
        """Mock back-translation via most-frequent human codons."""
        return _ext.back_translate(aa, unknown)

except Exception:  # pragma: no cover - pure-Python fallback
    _comp = str.maketrans("ACGTacgtN", "TGCATGCAN")

    def translate(nt: str, frame: int = 0) -> str:
        """Translate a nucleotide string from ``frame`` (0/1/2)."""
        s = nt.upper()
        table = CODON_TABLE
        return "".join(table.get(s[i : i + 3], "X") for i in range(frame, len(s) - 2, 3))

    def detect_coding_frame(nt: str) -> int:
        """Return the reading frame (0/1/2) with the fewest stop codons."""
        best_frame, best_stops = 0, None
        for f in (0, 1, 2):
            stops = translate(nt, f).count("*")
            if best_stops is None or stops < best_stops:
                best_stops, best_frame = stops, f
                if stops == 0:
                    break
        return best_frame

    def reverse_complement(nt: str) -> str:
        """Reverse-complement a nucleotide string (non-ACGT -> ``N``)."""
        return nt.upper().translate(_comp)[::-1]

    def back_translate(aa: str, unknown: str = "NNN") -> str:
        """Mock back-translation via most-frequent human codons."""
        return "".join(_HUMAN_CODON.get(c, unknown) for c in aa)


def aa_coords_from_nt(nt_start: int, nt_end: int, coding_start: int) -> tuple[int, int]:
    """Map a 1-based closed nt interval to 1-based closed aa coordinates.

    Args:
        nt_start: 1-based start of the region in the nucleotide sequence.
        nt_end: 1-based end (closed).
        coding_start: 1-based nt position where translation begins (frame origin).

    Returns:
        ``(aa_start, aa_end)`` 1-based closed, in the translated protein.
    """
    aa_start = (nt_start - coding_start) // 3 + 1
    aa_end = (nt_end - coding_start) // 3 + 1
    return aa_start, aa_end
