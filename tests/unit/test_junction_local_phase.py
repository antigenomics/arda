"""Local junction rendering must not bridge an upstream coding-frame offset."""
import pytest

from arda.annotate.transfer import _junction_nt


@pytest.mark.parametrize("coding_start", [1, 2, 3])
def test_upstream_phase_does_not_change_local_junction_translation(coding_start):
    nucleotide, protein, cdr3, phase = _junction_nt("TGTAAATGG", 4, 7, coding_start, 3)
    assert nucleotide == "TGTAAATGG"
    assert protein == "CKW"
    assert cdr3 == "K"
    assert phase == (7 - coding_start) % 3


@pytest.mark.parametrize("inserted", ["A", "AA"])
@pytest.mark.parametrize("coding_start", [1, 2, 3])
def test_local_frameshift_still_bridges_into_the_terminal_anchor(inserted, coding_start):
    sequence = "TGTAAA" + inserted + "TGG"
    nucleotide, protein, cdr3, phase = _junction_nt(
        sequence, 4, len(sequence) - 2, coding_start, 6)
    assert nucleotide == sequence
    assert "_" in protein
    assert protein.startswith("C") and protein.endswith("W")
    assert cdr3 == protein[1:-1]
    assert phase == (len(sequence) - 2 - coding_start) % 3
