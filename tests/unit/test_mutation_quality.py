"""``map --mutation-quality``: the Phred behind each ``v_mutations`` / ``j_mutations`` entry.

A novel allele, somatic hypermutation and a base miscall are the same string in the mutation
list. What separates them is how often the mutation recurs across an allele's reads and how good
the base is; the second needs this column.

⛔ The failure this file exists to catch is a quality list of the RIGHT LENGTH holding the wrong
bases' scores. Two ways to get there, and both are asserted against:

* **Re-deriving the mutation list instead of reading it.** Walking the alignment and scoring every
  mismatch reproduces what ``_markup.segment_cigars`` found, which since 2.16.0 is a SUPERSET of
  the emitted columns -- ``arda.shm`` drops the junction-internal entries afterwards. On this
  repo's own real-read fixture that is 25 of 242 V rows.
* **Reading the wrong strand.** Quality belongs to the read as submitted; the alignment and every
  coordinate are on the coding strand.
"""

from __future__ import annotations

import pytest

from arda.rnaseq.map import MUTATION_QUALITY, mutation_quality
from arda.refbuild.translate import reverse_complement

_VQ, _JQ = MUTATION_QUALITY

# A 12 nt query against a scaffold whose V runs to target 12 and whose J starts at 20. Three V
# mismatches at germline 3, 7, 11; the query is the germline with those three bases swapped.
_GERM = "ACGTACGTACGT"
_QUERY = "ACATACCTACAT"                                  # positions 3, 7, 11 differ
_QUAL = "IIA!IIIIII+I"                                   # Phred 40,40,32,0,... ,10 at pos 11


def _rec(**kw) -> dict:
    rec = {
        "sequence": _QUERY, "rev_comp": "F",
        "sequence_alignment": _QUERY, "germline_alignment": _GERM,
        "mmseqs2_qstart": 1.0, "mmseqs2_tstart": 1.0,
        "mmseqs2_t_vend": 12, "mmseqs2_t_jstart": 20, "mmseqs2_t_vjend": 30,
        "v_mutations": "G3A,G7C,G11A", "j_mutations": "",
    }
    rec.update(kw)
    return rec


def _phred(qual: str, one_based: int) -> int:
    return ord(qual[one_based - 1]) - 33


def test_each_entry_gets_its_own_base_score():
    out = mutation_quality(_rec(), _QUAL)
    assert out[_VQ] == ",".join(str(_phred(_QUAL, p)) for p in (3, 7, 11))
    assert out[_JQ] == ""


def test_the_scores_are_not_all_the_same_number():
    """Guards the test itself: a walk that returned the FIRST base's score every time would pass a
    length check and a spot check on entry 0."""
    scores = mutation_quality(_rec(), _QUAL)[_VQ].split(",")
    assert len(set(scores)) == 3, scores


def test_a_subset_of_the_walk_is_scored_at_its_own_positions():
    """The emitted list drives the lookup. ``arda.shm`` scopes it AFTER `segment_cigars` ran, so a
    read whose column carries 1 of 3 mismatches must get that one's Phred, not the first one's."""
    out = mutation_quality(_rec(v_mutations="G11A"), _QUAL)
    assert out[_VQ] == str(_phred(_QUAL, 11))
    assert out[_VQ] != str(_phred(_QUAL, 3))


def test_reverse_complemented_reads_read_the_reversed_quality():
    rc_seq = reverse_complement(_QUERY)
    rc_qual = _QUAL[::-1]
    out = mutation_quality(_rec(sequence=rc_seq, rev_comp="T"), rc_qual)
    assert out[_VQ] == ",".join(str(_phred(_QUAL, p)) for p in (3, 7, 11))


def test_a_j_side_mutation_uses_j_germline_coordinates():
    """J germline position is ``target - t_jstart + 1``, so an entry at J position 1 is the base
    aligned to target 20 -- not to target 1."""
    germ = "ACGTACGTACGT" + "N" * 7 + "TTTT"
    query = "ACGTACGTACGT" + "N" * 7 + "ATTT"
    rec = _rec(sequence=query, sequence_alignment=query, germline_alignment=germ,
               v_mutations="", j_mutations="T1A")
    qual = "I" * 19 + "5" + "I" * 3
    out = mutation_quality(rec, qual)
    assert out[_VQ] == ""
    assert out[_JQ] == str(_phred(qual, 20))


@pytest.mark.parametrize("rec", [
    {"sequence_alignment": ""},                                    # no alignment
    {"v_mutations": "", "j_mutations": ""},                        # nothing to score
    {"v_mutations": "G99A"},                                       # position outside the alignment
    {"v_mutations": "Gxx A"},                                      # unparseable entry
    {"mmseqs2_qstart": ""},                                        # no coordinates
])
def test_it_returns_empty_rather_than_a_shorter_list(rec):
    """A short list would be read positionally against the mutation column and silently pair entry
    *i* with a different base's score. Refuse instead."""
    out = mutation_quality(_rec(**rec), _QUAL)
    assert out[_VQ] == "" and out[_JQ] == ""


def test_a_quality_string_of_the_wrong_length_is_refused():
    assert mutation_quality(_rec(), _QUAL[:-2]) == {_VQ: "", _JQ: ""}


def test_reconstruct_is_refused_rather_than_given_one_mates_quality(tmp_path):
    """A merged fragment's bases come from two reads, so no input quality string describes it."""
    from arda.rnaseq.map import map_rnaseq

    with pytest.raises(ValueError, match="mutation-quality"):
        map_rnaseq(tmp_path / "r1.fq", tmp_path / "out.tsv", r2=tmp_path / "r2.fq",
                   reconstruct=True, with_mutation_quality=True)
