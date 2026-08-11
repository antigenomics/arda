"""PCR/template-switch chimera flagging, and the germline trap that makes it hard.

⛔ THE TRAP, measured. A junction is `V 3' tail` + `N/P/D` + `J 5' head`, and both tails are
GERMLINE — every clonotype on a V starts with the same bases, every one on a J ends with them. On a
TRA amplicon the median V-templated prefix is 10 nt and the median J-templated suffix 25 nt against
a median clone-specific core of **5 nt**. So a prefix/suffix agreement test run on the raw junction
rediscovers the germline and calls **52.20 % of clonotypes chimeric** (35.50 % of reads). Excluding
the templated tails and requiring 6 non-templated nt each side takes that to **0.02 %**.

⛔ Cell Ranger's published rule (contigs sharing a V prefix ≥ 25 nt with different CDR3s) does not
port to bulk, and not because of the constant: it relies on the BARCODE PARTITION, where ~1 clone
per chain means a second V-sharing contig is an artifact. A polyclonal bulk repertoire has thousands
of real clones per V gene. What ports is UCHIME's shape — a query explained by two MORE ABUNDANT
parents — because the abundance ordering, not the partition, is what makes it a claim.
"""
from __future__ import annotations

import polars as pl
import pytest

from arda.rnaseq.correct import _flag_chimeras


def _frame(rows):
    """rows: (junction, count, v_call, j_call) -> the frame `_flag_chimeras` consumes."""
    return pl.DataFrame({
        "junction": [r[0] for r in rows],
        "duplicate_count": [r[1] for r in rows],
        "locus": ["TRA"] * len(rows),
        "v_call": [r[2] for r in rows],
        "j_call": [r[3] for r in rows],
    })


def _anchors(organism="human"):
    from arda.cdr3fix import load_anchors
    a = load_anchors(organism)
    if not a:
        pytest.skip("human reference not built")
    return a


def _real_junction(anchors, v, j, core):
    """A junction shaped like a real one: the V's germline tail + a clone-specific core + the J's
    germline head. Built from the SHIPPED anchors, so the germline spans are the real ones — a
    hand-written string would not exercise the exclusion at all."""
    vt = anchors[("V", v)].germline_nt
    jt = anchors[("J", j)].germline_nt
    return vt + core + jt


V, J1, J2 = "TRAV1-2*01", "TRAJ33*01", "TRAJ12*01"


def test_two_real_clones_sharing_ONLY_germline_are_not_chimeras():
    """⛔ THE 52 % TRAP, as a test that actually discriminates.

    ⚠ A first version of this test used three clones on one V and one J and **passed with the
    germline exclusion deleted** — the cores differed, so no suffix parent existed and the trap
    could not fire however wrong the code was. The construction below is the one that fires:

        Q = V1 germline + core + J germline
        A = V1 germline + otherCore + J germline    <- shares Q's V-TEMPLATED PREFIX, nothing else
        B = V2 germline + core + J germline         <- shares Q's core + J tail, differing only in
                                                       the V germline (a real call-split shape)

    Without the exclusion the breakpoint lands at the end of V1's germline: `Q[:vt]` is a prefix of
    A and `Q[vt:]` is a suffix of B, so Q is "explained" by two parents while sharing **no
    clone-specific sequence with either**. With the exclusion the breakpoint must be at least
    `min_specific` nt past the templated tail, where neither parent matches.
    """
    a = _anchors()
    V1, V2 = "TRAV11*01", "TRAV24*01"                 # same length, different germline
    core = "GGGCCCTTTAAACCC"
    rows = [
        (_real_junction(a, V1, J1, "TTTAAACCCGGGTTT"), 500, V1, J1),      # A
        (_real_junction(a, V2, J1, core), 300, V2, J1),                   # B
        (_real_junction(a, V1, J1, core), 5, V1, J1),                     # Q
    ]
    got = _flag_chimeras(_frame(rows), "human").to_list()
    assert got[2] == "", (
        f"flagged a real clone that shares only GERMLINE with each parent: {got[2]!r} — this is "
        f"the failure mode that calls 52.20 % of a repertoire chimeric"
    )


def test_a_constructed_chimera_IS_flagged():
    """The positive control: a query whose clone-specific core is parent A's first half joined to
    parent B's second half, both parents more abundant."""
    a = _anchors()
    core_a, core_b = "GGGCCCTTTAAA", "ACACACGTGTGT"
    rows = [
        (_real_junction(a, V, J1, core_a), 500, V, J1),
        (_real_junction(a, V, J1, core_b), 300, V, J1),
        (_real_junction(a, V, J1, core_a[:6] + core_b[6:]), 5, V, J1),   # the chimera
    ]
    got = _flag_chimeras(_frame(rows), "human").to_list()
    assert got[0] == "" and got[1] == "", "a parent was flagged as the chimera"
    assert got[2], "the constructed chimera was not flagged"
    assert "@" in got[2], f"flag should carry the breakpoint: {got[2]!r}"


@pytest.mark.parametrize("q_count,why", [
    (500, "the query is the MOST abundant of the three"),
    (300, "the query TIES its parents — a tie is not evidence either is the original"),
])
def test_parents_must_be_STRICTLY_more_abundant(q_count, why):
    """Abundance ordering is what turns 'shares a prefix and a suffix' into a chimera claim.

    ⚠ The strictly-greater test is not redundant with the descending sort, and getting a case that
    exercises it took two attempts. The prefix/suffix maps are filled most-abundant-first with
    `setdefault`, so a parent can never be *less* abundant than the query — but it can TIE it, and
    a tie is no reason to call either clone derived from the other. The tie only reaches the
    comparison when the query sorts AFTER its parents in the total order `(-count, junction)`, so
    the cores below are chosen to put Q last lexicographically; with Q sorting first it claims its
    own prefix, `a == qi`, and the check is never consulted however wrong it is.
    """
    a = _anchors()
    core_a, core_b = "GGGCCCAAAAAA", "ACACACGTGTGT"
    core_q = "GGGCCC" + core_b[6:]                    # sorts after core_a ('G' > 'A' at index 6)
    rows = [
        (_real_junction(a, V, J1, core_a), 300, V, J1),
        (_real_junction(a, V, J1, core_b), 300, V, J1),
        (_real_junction(a, V, J1, core_q), q_count, V, J1),
    ]
    assert _flag_chimeras(_frame(rows), "human").to_list()[2] == "", why


def test_a_point_mutant_is_not_a_chimera():
    """⛔ The SHM confound. On IG, hypermutation manufactures near-variants continuously; a query
    within 2 substitutions of a parent is already explained by the error model and must not be
    flagged, however the prefix/suffix search happens to land.

    ⚠ Both parents must genuinely match for this to test the guard: A supplies Q's prefix and B its
    suffix, so the search DOES find a breakpoint and only the Hamming check stops it. An earlier
    version had no suffix parent, so nothing was flagged whether the guard ran or not.
    """
    a = _anchors()
    core_a = "GGGCCCTTTAAACCC"
    core_q = "GGGCCCTTTAAACCG"                        # 1 substitution from core_a
    core_b = "ACACACTTTAAACCG"                        # shares Q's suffix, unrelated prefix
    rows = [
        (_real_junction(a, V, J1, core_a), 500, V, J1),
        (_real_junction(a, V, J1, core_b), 300, V, J1),
        (_real_junction(a, V, J1, core_q), 2, V, J1),
    ]
    assert _flag_chimeras(_frame(rows), "human").to_list()[2] == "", (
        "flagged a 1-substitution variant as a chimera — on IG that class is most of the table"
    )


def test_without_a_reference_it_flags_NOTHING():
    """⛔ Fail safe, not fail useful. Without anchors the germline cannot be excluded and the test
    would report ~52 %. `resolve_airr` already shipped the other behaviour once — degrading
    silently to an empty germline set and returning output that looked fine."""
    a = _anchors()
    rows = [(_real_junction(a, V, J1, c), n, V, J1) for c, n in
            (("GGGCCCTTTAAA", 500), ("ACACACGTGTGT", 300), ("GGGCCCGTGTGT", 5))]
    got = _flag_chimeras(_frame(rows), "no_such_organism_xyz").to_list()
    assert got == ["", "", ""], f"guessed without a reference: {got}"


def test_the_column_is_optional_and_named():
    s = _flag_chimeras(_frame([("TGTGCCTTT", 1, V, J1)]), "human")
    assert s.name == "chimera_parents"
    assert s.dtype == pl.Utf8


def test_a_call_split_is_not_a_chimera_of_itself():
    """⛔ Two clonotypes can share a junction BYTE-IDENTICALLY under different V/J calls — that is
    the call-split class (Jurkat: 130 of 14,531 reads, including an allele-level TRG split), and it
    is why `clonotype_key` exists. Such a twin claims every prefix AND every suffix of the query,
    so without the identical-parent guard the query is "explained" as a chimera of itself.
    """
    a = _anchors()
    jn = _real_junction(a, V, J1, "GGGCCCTTTAAACCC")
    rows = [
        (jn, 500, "TRAV11*01", J1),        # the abundant twin: same junction, different V call
        (_real_junction(a, V, J1, "ACACACGTGTGTACA"), 300, V, J1),
        (jn, 5, V, J1),                    # the query
    ]
    got = _flag_chimeras(_frame(rows), "human").to_list()
    assert got[2] == "", f"flagged a call split as a chimera of its own twin: {got[2]!r}"
