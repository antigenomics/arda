"""A stop codon in FR4 is a stop codon in the rearrangement.

`productive` and `stop_codon` scanned the V-side regions and the junction, and nothing else. The
junction ends AT [FW]118, which is FR4's first residue — so residues 2..n of the J were looked at
by neither, and a read whose J carried a stop came out `productive=T`, `stop_codon=F`.

Measured before the fix, on the committed IgBLAST fixtures: every evaluable row carries an FR4
(mean 10.1 aa human, 8.6 mouse) and **none of 4,217** contains a stop — which is why closing the
hole left all 4,843 annotated rows byte-identical. Real, and unexercised by real data, which is
exactly the kind of hole that needs a test written for it deliberately.
"""

import random

import pytest

from arda.annotate.mapper import annotate_records
from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]


@pytest.fixture(scope="module")
def clean_scaffold(human_scaffolds):
    """A scaffold arda annotates as productive, with an FR4 long enough to mutate inside."""
    random.seed(7)
    for sid, seq in human_scaffolds:
        rec = annotate_records([(sid, seq)], "human", "nt", threads=4)[0]
        if (rec["productive"] == "T" and rec["stop_codon"] == "F"
                and len(rec.get("fwr4_aa") or "") >= 4 and rec.get("fwr4_start")):
            return sid, seq, rec
    pytest.skip("no productive scaffold with a usable FR4 in the sample")


def _stop_at(seq: str, start_nt: int, codon_index: int) -> str:
    """Replace one whole codon of FR4 with TAA, in FR4's own frame."""
    at = (start_nt - 1) + 3 * codon_index
    return seq[:at] + "TAA" + seq[at + 3:]


def test_a_stop_inside_fr4_makes_the_read_non_productive(clean_scaffold):
    sid, seq, before = clean_scaffold
    assert before["productive"] == "T" and before["stop_codon"] == "F"

    # Codon 1, not 0: codon 0 is [FW]118 itself, which the junction already covers — mutating it
    # would be caught by the old rule and prove nothing about the new one.
    mutated = _stop_at(seq, int(before["fwr4_start"]), 1)
    after = annotate_records([("fr4stop", mutated)], "human", "nt", threads=4)[0]

    assert "*" in after["fwr4_aa"], "the mutation did not land in FR4's reading frame"
    assert "*" not in (after.get("junction_aa") or ""), "the stop leaked into the junction"
    assert after["stop_codon"] == "T"
    assert after["productive"] == "F"


def test_the_stop_is_the_only_thing_that_changed(clean_scaffold):
    """Never: `productive` must flip because of the stop, not because the mutation moved a
    boundary. If `vj_in_frame` or the junction changed too, the test is proving something else."""
    sid, seq, before = clean_scaffold
    after = annotate_records(
        [("fr4stop", _stop_at(seq, int(before["fwr4_start"]), 1))], "human", "nt", threads=4)[0]

    assert after["junction_aa"] == before["junction_aa"]
    assert after["vj_in_frame"] == before["vj_in_frame"] == "T"
    assert after["v_call"] == before["v_call"]
