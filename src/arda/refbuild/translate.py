"""Native nucleotide translation and reading-frame utilities.

No BioPython: a plain dict codon table is faster and dependency-free. Used both
to normalize V germline to its coding frame when building scaffolds and to
derive protein markup from nucleotide coordinates.
"""

from __future__ import annotations

__all__ = ["CODON_TABLE", "translate", "detect_coding_frame", "aa_coords_from_nt"]

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


def translate(nt: str, frame: int = 0) -> str:
    """Translate a nucleotide string from ``frame`` (0/1/2).

    Codons containing any non-ACGT base (e.g. ``N``) translate to ``X``; stop
    codons are ``*``.
    """
    s = nt.upper()
    out = []
    table = CODON_TABLE
    for i in range(frame, len(s) - 2, 3):
        out.append(table.get(s[i : i + 3], "X"))
    return "".join(out)


def detect_coding_frame(nt: str) -> int:
    """Return the reading frame (0/1/2) with the fewest stop codons.

    For a genuine germline V-REGION exactly one frame is stop-free; ties resolve
    to the lowest frame.
    """
    best_frame = 0
    best_stops = None
    for f in (0, 1, 2):
        stops = translate(nt, f).count("*")
        if best_stops is None or stops < best_stops:
            best_stops = stops
            best_frame = f
            if stops == 0:
                break
    return best_frame


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
