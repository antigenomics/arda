"""Tie lists: the germlines a read's alignment cannot rule out, and their library-wide ranking.

⛔ The premise, measured: on a Ramos library arda emitted **0 multi-gene tie lists in 504 calls**
against IgBLAST's **11.68 %**, and on the 104 reads where the two disagreed (`IGLV2-23` vs
`IGLV2-14`), aligning each read to BOTH germlines showed **59 of 60 fit identically** — typically
at identity 1.0000 over 63–70 nt. Neither tool was right; both were overconfident.
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.annotate.ties import TieResolver, _containing_py, rank_ties, resolve_airr

try:
    from arda import _denoise as _cpp
except ImportError:                                # pragma: no cover
    _cpp = None


# --- the C++ primitive against its Python reference ---------------------------------------------

@pytest.mark.skipif(_cpp is None, reason="extension not built")
def test_containing_matches_the_python_reference():
    cands = ["AAACCCGGGTTT", "TTTCCCGGGAAA", "AAACCCGGGAAA", "GGG", ""]
    for needle in ("CCCGGG", "AAACCC", "GGGTTT", "ZZZ", ""):
        assert list(_cpp.containing(needle, cands)) == _containing_py(needle, cands)


# --- membership ---------------------------------------------------------------------------------

_A = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC"       # 40 nt
_SHARED = _A[:35]


def _res(**extra):
    g = {"GENEA*01": _A, "GENEB*01": _SHARED + "GGGG", "GENEC*01": "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"}
    g.update(extra)
    return TieResolver(g)


def test_a_germline_containing_the_same_span_is_a_tie():
    r = _res()
    assert r.expand("GENEA*01", 1, 35) == "GENEA*01,GENEB*01"


def test_a_span_that_discriminates_stays_a_single_call():
    """The read reaches past the shared stretch, so the alignment DOES separate them."""
    r = _res()
    assert r.expand("GENEA*01", 1, 40) == "GENEA*01"


def test_a_short_span_gets_no_tie_list():
    """⛔ Below MIN_SPAN nearly every allele of a family contains its neighbours' stretch, so the
    'tie list' would be the whole family and say less than the single call it replaced."""
    r = _res()
    assert r.expand("GENEA*01", 1, 20) == "GENEA*01"


def test_missing_or_malformed_coordinates_leave_the_call_alone():
    r = _res()
    for gs, ge in ((None, 35), ("", ""), (0, 35), (35, 1), ("x", "y")):
        assert r.expand("GENEA*01", gs, ge) == "GENEA*01"


def test_a_runaway_tie_list_is_refused():
    """Above the cap the call is left as it was: an unusable call is worse than a wrong one."""
    many = {f"G{i}*01": _A for i in range(40)}
    r = TieResolver(many, max_ties=16)
    assert r.expand("G0*01", 1, 35) == "G0*01"


def test_output_order_does_not_depend_on_dict_order():
    a = TieResolver({"GENEA*01": _A, "GENEB*01": _SHARED + "GGGG"})
    b = TieResolver({"GENEB*01": _SHARED + "GGGG", "GENEA*01": _A})
    assert a.expand("GENEA*01", 1, 35) == b.expand("GENEA*01", 1, 35)


# --- the second pass ----------------------------------------------------------------------------

def test_the_library_consensus_leads_the_tie_list():
    calls = ["A,B", "A,B", "B", "B", "B"]
    assert rank_ties(calls) == ["B,A", "B,A", "B", "B", "B"]


def test_only_unambiguous_reads_vote():
    """⛔ A read whose own call is `A,B` cannot be evidence for A over B — counting it would let a
    common tie bootstrap itself, so the more confusable a pair is the more confidently it would
    elect one of them."""
    # A appears in 100 ties and never alone; B is named alone once. B must still win.
    calls = ["A,B"] * 100 + ["B"]
    assert rank_ties(calls)[0] == "B,A"


def test_ranking_uses_pre_expansion_evidence():
    """Ranking on the EXPANDED calls is self-defeating: expansion makes every read ambiguous, so
    the unambiguous-only rule has nothing to count and it degenerates to lexicographic order."""
    expanded = ["A,B", "A,B", "A,B"]
    evidence = ["B", "B", "A"]                      # what the tool actually called, pre-expansion
    assert rank_ties(expanded, evidence=evidence) == ["B,A", "B,A", "B,A"]
    assert rank_ties(expanded)[0] == "A,B", "without evidence it can only sort"


def test_ranking_never_changes_membership():
    calls = ["A,B,C", "C", "C"]
    out = rank_ties(calls)
    assert {x for x in out[0].split(",")} == {"A", "B", "C"}
    assert out[0].startswith("C")


def test_ranking_is_deterministic_under_row_permutation():
    a = rank_ties(["A,B", "B,A"])
    b = rank_ties(["B,A", "A,B"])
    assert a[0].split(",")[0] == b[0].split(",")[0]


def test_scores_must_match_the_calls():
    with pytest.raises(ValueError, match="scores"):
        rank_ties(["A,B"], scores=[1.0, 2.0])


# --- end to end ---------------------------------------------------------------------------------

def test_resolve_airr_expands_then_ranks(tmp_path):
    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["r1", "r2", "r3"],
        "v_call": ["IGLV2-23*01", "IGLV2-23*01", "IGLV2-14*01"],
        "v_germline_start": ["1", "1", "1"],
        "v_germline_end": ["70", "70", "70"],
        "mmseqs2_score": ["100", "100", "100"],
    }).write_csv(src, separator="\t", quote_style="never")
    out = tmp_path / "out.tsv"
    rep = resolve_airr(src, out, organism="human", segments=("v",))
    got = pl.read_csv(out, separator="\t", infer_schema_length=0).to_dicts()
    # The first 70 nt of IGLV2-14 and IGLV2-23 are identical, so every allele of both is a tie...
    assert all("IGLV2-14" in r["v_call"] and "IGLV2-23" in r["v_call"] for r in got)
    # ...and IGLV2-23 leads everywhere, because it has 2 unambiguous votes against IGLV2-14's 1.
    assert all(r["v_call"].startswith("IGLV2-23") for r in got)
    assert rep["expanded"]["v"] == 3 and rep["reranked"]["v"] == 1
