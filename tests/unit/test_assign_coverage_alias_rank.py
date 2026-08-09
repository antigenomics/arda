"""⛔ A coverage ALIAS must never outrank a real root in the k-mer posting cap.

`_assign_coverage` bounds each k-mer's postings at `cap` and fills them in descending abundance.
An ALIAS -- a junction the quality gate vacated, kept in the index only so partial reads that
covered it still reach the parent -- is ordered by its PARENT's count, which is high by
construction. So aliases sorted to the FRONT and evicted genuine low-abundance roots, and the
partial reads of those roots were either lost or credited to the alias's parent instead.

Measured at full depth on a TRA amplicon (SRR5233636, --ec-mode accurate, 23,360 aliases against
36,587 roots; every arm emits an IDENTICAL 36,587-row clonotype table, so this is purely read
ASSIGNMENT):

    cap 64,   aliases ordered by abundance   1,812,740   -25,473 vs --ec-mode fast
    cap 64,   aliases OFF                    1,838,181       -32
    cap 64,   ROOTS FIRST (the fix)          1,841,624    +3,411
    cap 1024, aliases ordered by abundance   1,869,556   +31,343

The mechanism was added to rescue 5 reads of 9,208 on Ramos and was costing 25,441 here.

⚠ `cap` is a keyword argument, so the ordering is pinned with cap=1 rather than a 64-target
fixture. The alias must also SHARE the contested k-mers -- an earlier version of this test used an
alias with no k-mer in common with the root, so nothing competed and it passed with the bug in.
"""

from __future__ import annotations

import polars as pl

from arda.rnaseq.correct import _assign_coverage

# Two 48 nt junctions. RARE is a surviving root with ONE spanning read; ABUNDANT is a big root.
RARE = "TGTGCCAGCAGTTTCTCGACCTGTTCGGCTAACTATGGCTACACCTTC"
ABUNDANT = "TGTGCCAGCAGTAAACTGGGGACAGGGCCCTTAGCAGTTTTCCCTTTC"


def _reads(seqs: dict[str, str]) -> pl.DataFrame:
    n = len(seqs)
    return pl.DataFrame({
        "sequence_id": list(seqs),
        "sequence": list(seqs.values()),
        "locus": ["TRB"] * n,
        "v_call": ["TRBV12-3*01"] * n,
        "j_call": ["TRBJ1-2*01"] * n,
        "c_call": ["TRBC1"] * n,
        "rev_comp": ["F"] * n,
        "junction": [""] * n,          # no junction of its own: only the alignment pass can place it
    })


def test_an_abundant_alias_must_not_steal_the_kmer_slots_of_a_rare_root():
    """The alias carries RARE's own junction (a call split the gate emptied) but points at the
    ABUNDANT root, so it is ordered by ABUNDANT's count and, at cap=1, wins every k-mer RARE needs.
    The read then lands on the wrong clonotype."""
    got = _assign_coverage(
        _reads({"r1": RARE}), [RARE, ABUNDANT], ["TRB", "TRB"], {},
        aliases=[(RARE, "TRB", 1)], root_counts=[1, 5000], cap=1)
    assert got[0] == ["r1"], "a read covering the RARE root was credited to the alias's parent"
    assert got[1] == []


def test_an_alias_still_places_a_read_that_only_covers_the_vacated_junction():
    """Roots-first must not disable aliases -- they are worth +3,411 reads at full depth."""
    got = _assign_coverage(
        _reads({"r1": RARE}), [ABUNDANT], ["TRB"], {},
        aliases=[(RARE, "TRB", 0)], root_counts=[5000], cap=64)
    assert got[0] == ["r1"], "a read covering only the vacated junction should reach the parent"
