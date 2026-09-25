"""Tie lists: the germlines a read's alignment cannot rule out, and their library-wide ranking.

Never: The premise, measured: on a Ramos library arda emitted **0 multi-gene tie lists in 504 calls**
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
    """Never: Below MIN_SPAN nearly every allele of a family contains its neighbours' stretch, so the
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
    """Never: A read whose own call is `A,B` cannot be evidence for A over B — counting it would let a
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

def _reference_available() -> bool:
    """The COMMITTED reference, not the IMGT source tree.

    ⚠ This gate used to ask ``refbuild.imgt.load_functional_alleles`` whether the IMGT *source*
    directory was present -- which ``arda build-db`` downloads and neither CI nor a pip install
    has. So every test below it skipped in CI, silently, which is how ``resolve_airr`` shipped
    reading a germline source that did not exist on any plain install. ``markup.tsv`` and
    ``alleles.fasta`` are committed, so gating on them means these run everywhere.
    """
    from arda.paths import vdj_dir
    try:
        base = vdj_dir("human")
    except Exception:                              # no database at all
        return False
    return (base / "markup.tsv").exists() and (base / "alleles.fasta").exists()


def test_resolve_airr_raises_when_the_germlines_are_absent(tmp_path, monkeypatch):
    """Never: Raise, never degrade.

    Without germlines every call is left exactly as it was, so the output is byte-identical to the
    input and the report reads ``expanded: 0`` -- indistinguishable from a library that genuinely
    had no ties. That is the failure mode this project has shipped before, and here it kept a CI
    test red: the reference is not built in CI, so ``resolve-ties`` silently did nothing.
    """
    from arda import germline as germline_mod
    from arda.annotate.ties import resolve_airr as _resolve

    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["r1"], "v_call": ["IGLV2-23*01"],
        "v_germline_start": ["1"], "v_germline_end": ["70"],
    }).write_csv(src, separator="\t", quote_style="never")

    def _none(*a, **k):
        raise FileNotFoundError("no reference here")

    # `resolve_airr` imports the helper locally, so patch it at the source module.
    monkeypatch.setattr(germline_mod, "segment_germlines", _none)
    with pytest.raises(ValueError, match="germline"):
        _resolve(src, tmp_path / "out.tsv", organism="human", segments=("v",))


@pytest.mark.skipif(not _reference_available(), reason="no built reference")
def test_resolve_airr_does_not_touch_the_imgt_source_tree(tmp_path, monkeypatch):
    """Never: the germlines come from the COMMITTED reference, never from ``data/imgt``.

    ``resolve_airr`` used to read ``refbuild.imgt.load_functional_alleles``, whose input is the
    IMGT source directory that ``arda build-db`` downloads. Neither the wheel nor the reference
    tarball carries it, so ``arda resolve-ties`` **raised on every plain pip install** -- the same
    shape as the ``segments.fasta`` deploy trap, and invisible in a source checkout because one is
    sitting right there.

    So: break the IMGT path on purpose and require the command to work anyway.
    """
    from arda.annotate.ties import resolve_airr as _resolve
    from arda.refbuild import imgt as imgt_mod

    def _boom(*a, **k):
        raise AssertionError("resolve_airr must not read the IMGT source tree")

    monkeypatch.setattr(imgt_mod, "load_functional_alleles", _boom)
    monkeypatch.setattr(imgt_mod, "ungap_gene", _boom)

    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["r1"], "v_call": ["IGLV2-23*01"],
        "v_germline_start": ["1"], "v_germline_end": ["70"],
    }).write_csv(src, separator="\t", quote_style="never")

    report = _resolve(src, tmp_path / "out.tsv", organism="human", segments=("v",))
    assert report["rows"] == 1


def test_expansion_and_ranking_end_to_end_without_the_reference():
    """The whole pipeline -- expand, then rank -- on an explicit germline set.

    ⚠ Hermetic ON PURPOSE. The previous version of this test built its germlines from the installed
    IMGT reference, so it passed locally and failed in CI (which does not build one) for six
    releases' worth of commits before anyone looked. The biology it encodes is real and measured:
    the first 70 nt of ``IGLV2-14`` and ``IGLV2-23`` are identical, so a read aligned over that span
    cannot choose between them.
    """
    shared = "CAGTCTGCCCTGACTCAGCCTGCCTCCGTGTCTGGGTCTCCTGGACAGTCGATCACCATCTCCTGCACT"
    germ = {
        "IGLV2-14*01": shared + "GGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAG",
        "IGLV2-23*01": shared + "GGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAA",
        "IGLV3-1*01":  "TCCTATGAGCTGACACAGCCACCCTCGGTGTCAGTGTCCCCAGGACAAACGGCCAGGATCACCTGCTCT",
    }
    res = TieResolver(germ)
    # A read aligned over germline 1..70 cannot separate the two IGLV2 alleles...
    widened = [res.expand(c, 1, 70) for c in
               ("IGLV2-23*01", "IGLV2-23*01", "IGLV2-14*01")]
    for w in widened:
        assert "IGLV2-14*01" in w and "IGLV2-23*01" in w
    assert all("IGLV3-1*01" not in w for w in widened), "an unrelated germline must not join"

    # ...and the library then ranks them: IGLV2-23 has 2 unambiguous votes against IGLV2-14's 1.
    ranked = rank_ties(widened, None, evidence=["IGLV2-23*01", "IGLV2-23*01", "IGLV2-14*01"])
    assert all(r.startswith("IGLV2-23*01") for r in ranked)
    # Membership is unchanged -- only the order moves.
    assert all(sorted(a.split(",")) == sorted(b.split(",")) for a, b in zip(widened, ranked))


@pytest.mark.skipif(not _reference_available(), reason="no committed reference")
def test_resolve_airr_expands_then_ranks_against_the_real_reference(tmp_path):
    """The same claim against the SHIPPED germlines, skipped only when they are genuinely absent.

    ⚠ The skip is explicit and visible in the report, not a silent pass: the hermetic test above
    covers the logic unconditionally, so this one adds the reference and nothing else.
    """
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
    assert all("IGLV2-14" in r["v_call"] and "IGLV2-23" in r["v_call"] for r in got)
    assert all(r["v_call"].startswith("IGLV2-23") for r in got)
    assert rep["expanded"]["v"] == 3 and rep["reranked"]["v"] == 1


# --- candidates(): the refusal a restriction has to be able to see ------------------------------

def test_candidates_separates_a_refusal_from_a_unique_answer():
    """Never: ``None`` is not ``()``.

    ``expand`` returns the call unchanged for BOTH "the span was too short to say anything" and
    "the span is unique to this allele". Widening a call may conflate them; restricting one may
    not -- a genotype restriction that reads a refusal as a contradiction deletes the read's call.
    """
    r = _res()
    assert r.candidates("GENEA*01", 1, 40) == ("GENEA*01",)      # a real, unique answer
    assert r.candidates("GENEA*01", 1, 35) == ("GENEA*01", "GENEB*01")
    assert r.candidates("GENEA*01", 1, 20) is None               # below MIN_SPAN
    assert r.candidates("NOSUCH*01", 1, 40) is None              # not in this germline set
    assert r.candidates("", 1, 40) is None
    assert r.candidates("GENEA*01", None, 40) is None            # unusable coordinates


def test_a_runaway_tie_list_is_a_refusal_not_an_empty_answer():
    """Above ``max_ties`` there is no answer -- which is not the same as "compatible with nothing"."""
    many = {f"G{i}*01": _A for i in range(40)}
    assert TieResolver(many, max_ties=16).candidates("G0*01", 1, 35) is None


def test_expand_is_unchanged_by_the_refactor():
    """``expand`` is now a thin wrapper over ``candidates``; its contract must not have moved."""
    r = _res()
    assert r.expand("GENEA*01", 1, 35) == "GENEA*01,GENEB*01"
    assert r.expand("GENEA*01", 1, 40) == "GENEA*01"
    assert r.expand("GENEA*01", 1, 20) == "GENEA*01"
    assert r.expand("GENEA*01,GENEZ*99", 1, 35) == "GENEA*01,GENEZ*99,GENEB*01"
    for gs, ge in ((None, 35), ("", ""), (0, 35), (35, 1), ("x", "y")):
        assert r.expand("GENEA*01", gs, ge) == "GENEA*01"


# --- `--loci`: widening helps on some loci and hurts on others ---------------------------------
# Measured against an IgBLAST truth on two bulk libraries (benchmark round 32), exact `v_gene` SET
# agreement moves IGK .5707 -> .9373 and .5558 -> .9308, IGL .8586 -> .9137 and .8577 -> .9115,
# TRB .9299 -> .9694 -- and IGH .8255 -> .7238 and .8313 -> .7344, a clear loss on every IGH arm
# measured. So the useful control is per LOCUS, and a mixed library is the normal case.

@pytest.mark.skipif(not _reference_available(), reason="no committed reference")
def test_loci_widens_only_the_named_loci(tmp_path):
    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["igl1", "igl2", "igh1"],
        "locus": ["IGL", "IGL", "IGH"],
        "v_call": ["IGLV2-23*01", "IGLV2-14*01", "IGLV2-23*01"],
        "v_germline_start": ["1", "1", "1"],
        "v_germline_end": ["70", "70", "70"],
        "mmseqs2_score": ["100", "100", "100"],
    }).write_csv(src, separator="\t", quote_style="never")
    out = tmp_path / "out.tsv"
    rep = resolve_airr(src, out, organism="human", segments=("v",), loci=("IGL",))
    got = {r["sequence_id"]: r["v_call"] for r in
           pl.read_csv(out, separator="\t", infer_schema_length=0).to_dicts()}
    assert "IGLV2-14" in got["igl1"] and "IGLV2-23" in got["igl1"]
    # The IGH row carries the SAME call and the same span, so only the locus filter can be what
    # leaves it alone.
    assert got["igh1"] == "IGLV2-23*01"
    assert rep["expanded"]["v"] == 2 and rep["loci"] == ["IGL"]


@pytest.mark.skipif(not _reference_available(), reason="no committed reference")
def test_loci_is_case_and_space_insensitive(tmp_path):
    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["r1"], "locus": ["IGL"], "v_call": ["IGLV2-23*01"],
        "v_germline_start": ["1"], "v_germline_end": ["70"], "mmseqs2_score": ["100"],
    }).write_csv(src, separator="\t", quote_style="never")
    out = tmp_path / "out.tsv"
    rep = resolve_airr(src, out, organism="human", segments=("v",), loci=(" igl ",))
    assert rep["expanded"]["v"] == 1


def test_loci_without_a_locus_column_raises_rather_than_widening_everything(tmp_path):
    """Never: a filter that cannot be applied must not fall back to applying nothing.

    Silently widening every row would look exactly like success and would be the IGH regression
    the caller asked to avoid.
    """
    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["r1"], "v_call": ["IGLV2-23*01"],
        "v_germline_start": ["1"], "v_germline_end": ["70"],
    }).write_csv(src, separator="\t", quote_style="never")
    with pytest.raises(ValueError, match="no `locus` column"):
        resolve_airr(src, tmp_path / "out.tsv", organism="human", segments=("v",),
                     loci=("IGK", "IGL"))


@pytest.mark.skipif(not _reference_available(), reason="no committed reference")
def test_no_loci_is_every_locus_exactly_as_before(tmp_path):
    src = tmp_path / "in.tsv"
    pl.DataFrame({
        "sequence_id": ["igl1", "igh1"],
        "locus": ["IGL", "IGH"],
        "v_call": ["IGLV2-23*01", "IGLV2-23*01"],
        "v_germline_start": ["1", "1"], "v_germline_end": ["70", "70"],
        "mmseqs2_score": ["100", "100"],
    }).write_csv(src, separator="\t", quote_style="never")
    out = tmp_path / "out.tsv"
    rep = resolve_airr(src, out, organism="human", segments=("v",))
    assert rep["expanded"]["v"] == 2 and "loci" not in rep
