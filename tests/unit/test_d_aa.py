"""D mapping on amino-acid input, and D mapping at the clonotype level.

``_markup.d_local_align`` is a plain character comparison, so the same gapless local aligner
that finds a D in a nucleotide interior finds it in a translated one. What changes is the
alphabet: lambda = ln((1-p)/p) with p the chance two residues match, which is 1/4 for nt
(recovering ln 3) and a *measured* 0.0613 for aa -- N-region inserts and D germlines are both
G/S/Y-rich, so 1/20 would be wrong.

A trimmed D has no knowable reading frame, so the aa reference carries all three translations
of every allele. Two of them can tie on one span; that is one allele, not an ambiguity.
"""

from __future__ import annotations

import pytest

from arda.annotate.reference import _load_d_germlines, _load_d_germlines_aa
from arda.annotate.transfer import _best_d, _d_db_nt, _d_evalue, _d_gate, _d_min_score, _map_d
from arda.paths import vdj_dir
from arda.refbuild.translate import translate


def _gene(call: str) -> str:
    return call.split(",")[0].split("*")[0]


@pytest.fixture(scope="module")
def igh_aa():
    return _load_d_germlines_aa(vdj_dir("human"))["IGH"]


def test_aa_reference_holds_three_frames_per_allele(igh_aa):
    nt = _load_d_germlines(vdj_dir("human"))["IGH"]
    assert len(igh_aa) == 3 * len(nt), "one entry per (allele, reading frame)"
    by_allele = {}
    for allele, aa in igh_aa:
        by_allele.setdefault(allele, []).append(aa)
    allele, seq = nt[0]
    assert by_allele[allele] == [translate(seq[f:], 0) for f in (0, 1, 2)]


def test_the_aa_gate_uses_a_measured_lambda_not_a_uniform_one():
    import math

    lam_nt, e_nt, _ = _d_gate("nt")
    lam_aa, e_aa, _ = _d_gate("aa")
    assert lam_nt == pytest.approx(math.log(3.0))          # p = 1/4 => ln((1-p)/p) = ln 3
    assert lam_aa == pytest.approx(2.7285, abs=1e-4)       # p = 0.0613, measured
    assert lam_aa < math.log(19), "uniform 1/20 would overstate lambda; composition is biased"
    assert e_aa < e_nt, "small aa databases leave the E-value under-calibrated; gate tighter"


def test_aa_evalue_and_min_score_are_alphabet_aware():
    # Same score, same m and n: an aa hit is far less likely by chance than an nt one.
    assert _d_evalue(6, 30, 100, "aa") < _d_evalue(6, 30, 100, "nt")
    # ... so a shorter match suffices to clear the same E-value.
    assert _d_min_score(30, 100, 0.05, "aa") <= _d_min_score(30, 100, 0.05, "nt")


@pytest.mark.parametrize("seqtype", ["nt", "aa"])
def test_the_gate_is_monotone_in_both_directions(seqtype):
    """A better score is rarer; a tighter E-value demands a better score. The gate rests on it."""
    evalues = [_d_evalue(s, 30, 300, seqtype) for s in range(4, 11)]
    assert evalues == sorted(evalues, reverse=True) and len(set(evalues)) == len(evalues)
    thresholds = [_d_min_score(30, 300, e, seqtype) for e in (0.5, 0.2, 0.1, 0.05, 0.01)]
    assert thresholds == sorted(thresholds)


def test_the_gate_refuses_a_degenerate_search_rather_than_dividing_by_zero():
    huge = 1 << 30
    assert _d_min_score(0, 300) == huge          # no interior to search
    assert _d_min_score(30, 0) == huge           # no database
    assert _d_min_score(30, 300, 0.0) == huge    # an E-value of zero is unreachable
    assert _d_evalue(0, 30, 300) == float("inf")
    assert _d_evalue(6, 0, 300) == float("inf")


@pytest.mark.parametrize("seqtype", ["nt", "aa"])
def test_d_support_lets_a_consumer_reproduce_the_gate(seqtype):
    """``d_support`` is the E-value at the called score, so re-thresholding on it is exact.

    A consumer who keeps rows with ``d_support <= _D_MAX_EVALUE`` must keep exactly the rows
    arda called. That only holds if the emitted number is the same quantity the gate used.
    """
    m, n = 30, 300
    _, max_e, _ = _d_gate(seqtype)
    min_score = _d_min_score(m, n, seqtype=seqtype)
    assert _d_evalue(min_score, m, n, seqtype) <= max_e, "the accepted score clears the gate"
    assert _d_evalue(min_score - 1, m, n, seqtype) > max_e, "one point less would not have"


def test_a_planted_d_frame_is_found_in_an_aa_interior(igh_aa):
    """Plant one translated D germline verbatim; the caller must recover its gene and span."""
    nt = dict(_load_d_germlines(vdj_dir("human"))["IGH"])
    allele = "IGHD3-10*01"
    frame_aa = translate(nt[allele], 0)
    np1, np2 = "TP", "YFE"
    vpref, jsuf = "CARDS", "WGQGT"
    query = vpref + np1 + frame_aa + np2 + jsuf
    v_end, j_start = len(vpref), len(vpref) + len(np1 + frame_aa + np2) + 1

    rec: dict = {}
    _map_d(rec, query, v_end, j_start, igh_aa, "IGHJ4*02", seqtype="aa")
    assert _gene(rec["d_call"]) == "IGHD3-10", rec["d_call"]
    assert rec["np1"] == np1 and rec["np2"] == np2
    span = query[int(rec["d_sequence_start"]) - 1 : int(rec["d_sequence_end"])]
    assert span == frame_aa
    # Offsets index a reading frame, not the D germline, so germline coords must be withheld.
    assert not rec.get("d_germline_start") and not rec.get("d_cigar")


def test_two_frames_of_one_allele_are_not_reported_as_an_ambiguity(igh_aa):
    """A palindromic-ish D can align in two frames at the same span. That is one allele."""
    nt = dict(_load_d_germlines(vdj_dir("human"))["IGH"])
    frame_aa = translate(nt["IGHD3-10*01"], 0)
    interior = "TP" + frame_aa + "YFE"
    hit = _best_d(interior, igh_aa, _d_min_score(len(interior), _d_db_nt(igh_aa), seqtype="aa"))
    assert hit is not None
    alleles = hit[2]
    assert len(alleles) == len(set(alleles)), f"duplicate allele in ambiguity list: {alleles}"


def test_vj_locus_gets_no_d_on_aa_input():
    aa = _load_d_germlines_aa(vdj_dir("human"))
    assert "TRA" not in aa and "IGK" not in aa, "VJ loci have no D germlines, in any alphabet"


def test_clonotype_d_is_called_on_the_corrected_junction():
    """`correct` maps D once per clonotype, from (junction, v_call, j_call) alone."""
    import polars as pl

    from arda.rnaseq.correct import _clonotype_d

    # A real human TRB junction: TRBV20-1 germline .. TRBD1 .. TRBJ2-1 germline.
    junction = "TGCAGTGCTAGAGA" + "TAC" + "GGGACAGGGGGC" + "AAC" + "CTCCTACAATGAGCAGTTCTTC"
    df = pl.DataFrame({"junction": [junction, "TGTGCTTTTTTT"],
                       "v_call": ["TRBV20-1*01", "TRAV1-1*01"],
                       "j_call": ["TRBJ2-1*01", "TRAJ33*01"]})
    cols = {s.name: s for s in _clonotype_d(df, "human")}
    assert cols["d_call"][0].startswith("TRBD1"), cols["d_call"][0]
    assert float(cols["d_support"][0]) < 0.2
    assert cols["d_call"][1] == "" and cols["d2_call"][1] == "", "TRA is a VJ locus: no D"
    # The markup columns come with the call: np1 + D + np2 must rebuild the interior, so a
    # consumer of the CLONOTYPE table (which has no read) can cut the junction up unaided.
    lo, hi = cols["d_sequence_start"][0], cols["d_sequence_end"][0]
    interior = junction[cols["v_sequence_end"][0] : cols["j_sequence_start"][0] - 1]
    assert cols["np1"][0] + junction[lo - 1 : hi] + cols["np2"][0] == interior
    # A VJ-locus row carries no coordinates at all rather than a bogus zero-length span.
    assert cols["d_sequence_start"][1] == -1 and cols["np1"][1] == ""


def test_the_d_evalue_cli_default_does_not_loosen_the_aa_gate():
    """⛔ `--d-max-evalue` must default to None, never to the literal nt value.

    The calibrated operating point is alphabet-dependent -- 0.2 for nt, 0.05 for aa (the aa D
    database is 22-38 residues per TR locus and under-calibrated, so it gates tighter). A CLI
    default of 0.2 looks like a no-op and silently loosens every `--seqtype aa` run by 4x.
    """
    import inspect

    from arda import cli

    commands = [cli.annotate, cli.rnaseq_map, cli.rnaseq_correct, cli.rnaseq_assemble,
                cli.rnaseq_mode, cli.amplicon_mode, cli.rnaseq_reduce]
    for fn in commands:
        param = inspect.signature(fn).parameters["d_max_evalue"]
        assert param.default.default is None, f"{fn.__name__} defaults to {param.default.default}"

    # ...and None really does resolve per alphabet, rather than to one shared constant.
    from arda.annotate.transfer import _D_AA_MAX_EVALUE, _D_MAX_EVALUE, _d_min_score

    assert _D_AA_MAX_EVALUE < _D_MAX_EVALUE
    assert _d_min_score(30, 300, None, "aa") == _d_min_score(30, 300, _D_AA_MAX_EVALUE, "aa")
    assert _d_min_score(30, 300, None, "nt") == _d_min_score(30, 300, _D_MAX_EVALUE, "nt")
