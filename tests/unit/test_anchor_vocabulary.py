"""A functional segment's anchor residue, across every shipped organism.

The junction is Cys104..[FW]118 inclusive, so a functional V's templated region opens with Cys and
a functional J's closes with Phe or Trp. That holds for 3,129 of the 3,130 functional entries that
carry a templated region across human, mouse, rat, rhesus and rabbit.

⛔ The one exception is REAL GERMLINE, not a defect to repair: `TRAJ35*01` is an IMGT **functional**
gene whose anchor codon is `TGC` -- it decodes Cys, and `germline_nt` translates to `IGFGNVLHC` at
offset 1 and to nothing plausible at 0 or 2. Turning this assertion into a loader check that raises,
or into a `[FW]GXG` motif filter, silently deletes a functional gene from the vocabulary; every read
from it is then called as the nearest J that survived. That is why `cdr3fix` reads `anchor_nt` from
`cdr3_anchors.tsv` and never tests a motif, and why this lives in a test with the exception NAMED
rather than in `load_anchors`.

So the job here is to catch the SECOND one. A new non-canonical functional entry is either an IMGT
reclassification worth knowing about or a `refbuild` frame defect, and both want a human look.
"""

from __future__ import annotations

import pytest

from arda.cdr3fix import load_anchors

#: The shipped reference vocabulary, all five organisms.
ORGANISMS = ("human", "mouse", "rat", "rhesus_monkey", "rabbit")

#: Functional entries whose anchor residue is not the canonical one, each checked against IMGT.
#: Nothing may join this list without the three checks in CLAUDE.md, "Reference vocabulary".
KNOWN_NON_CANONICAL = {("human", "J", "TRAJ35*01")}


def test_a_functional_segment_carries_its_canonical_anchor():
    """Never: do not convert this into a loader assertion -- it would delete `TRAJ35*01`."""
    seen, outliers = 0, set()
    for org in ORGANISMS:
        for (seg, name), anchor in load_anchors(org).items():
            templated = anchor.templated_aa or ""
            if anchor.functionality != "F" or not templated:
                continue
            seen += 1
            canonical = templated[0] == "C" if seg == "V" else templated[-1] in "FW"
            if not canonical:
                outliers.add((org, seg, name))
    if seen == 0:
        pytest.skip("reference not built")
    assert seen > 3_000, f"only {seen} functional entries -- the reference looks truncated"
    assert outliers == KNOWN_NON_CANONICAL, (
        f"the set of non-canonical functional anchors moved: "
        f"new {sorted(outliers - KNOWN_NON_CANONICAL)}, gone {sorted(KNOWN_NON_CANONICAL - outliers)}")


def test_the_known_exception_is_readable_as_germline_and_stays_in_the_vocabulary():
    """`TRAJ35*01` must survive, anchor and all -- a motif check is what would drop it."""
    j = load_anchors("human").get(("J", "TRAJ35*01"))
    if j is None:
        pytest.skip("reference not built")
    assert j.functionality == "F" and j.templated_aa == "IGFGNVLHC"
    assert j.germline_nt[j.anchor_nt:j.anchor_nt + 3] == "TGC"    # Cys, and IMGT means it
