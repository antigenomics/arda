"""Locus disambiguation for dual-use TRAV/DV V genes (αβ vs γδ).

Genes named ``TRAV*/DV*`` (e.g. TRAV14/DV4, TRAV29/DV5) can rearrange in either the
α locus (with a TRAJ) or the δ locus (with a TRDJ). The locus must therefore follow the
J segment, not the (ambiguous) V gene. These tests build a V domain from such a gene with
a TRDJ-style vs a TRAJ-style CDR3/FR4 and assert arda assigns TRD vs TRA accordingly.
"""

import pytest

from arda.annotate.mapper import annotate_records
from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

# Shared TRAV14/DV4 V-region (FR1..FR3) taken from the human TRD_5 reference scaffold.
_V = (
    "QKITQTQPGMFVQEKEAVTLDCTYD"  # FR1
    "TSDQSYG"  # CDR1
    "LFWYKQPSSGEMIFLIY"  # FR2
    "QGSYDEQN"  # CDR2
    "ATEGRYSLNFQKARKSANLVISASQLGDSAMYF"  # FR3
)
# CDR3 (without the leading C carried by FR3's terminal F→C junction) + FR4.
_DELTA = "ALGDPSGGYTDKLIF" + "FGKGTRVTVEP"  # TRDJ1 ending (…TDKLIF) + δ FR4
_ALPHA = "ALGDPSGGYNQGGKLIF" + "FGGGTKLIIKP"  # TRAJ-style ending + α FR4


@pytest.mark.parametrize(
    "name,query,exp_locus,exp_jprefix",
    [
        ("delta", _V + _DELTA, "TRD", "TRDJ"),
        ("alpha", _V + _ALPHA, "TRA", "TRAJ"),
    ],
)
def test_trav_dv_locus_follows_j(name, query, exp_locus, exp_jprefix):
    rec = annotate_records([(name, query)], organism="human", seqtype="aa")[0]
    assert rec["locus"] == exp_locus, f"{name}: locus {rec['locus']} (v={rec.get('v_call')})"
    assert "/DV" in (rec.get("v_call") or ""), f"{name}: expected a TRAV/DV V gene"
    assert (rec.get("j_call") or "").startswith(exp_jprefix), f"{name}: j={rec.get('j_call')}"
