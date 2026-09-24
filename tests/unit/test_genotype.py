"""Genotype inference and call restriction.

⚠ **Hermetic ON PURPOSE**, following ``test_ties.py``: the germlines are synthetic and injected, so
these run with no reference, no mmseqs and no network. The same file records why -- an earlier
version of the tie tests built its germlines from an installed reference, passed locally, and
failed in CI for six releases before anyone looked.
"""

from __future__ import annotations

import random

import polars as pl
import pytest

from arda.genotype import (GENOTYPE_COLUMNS, gene_of, infer_genotype, read_genotype, restrict,
                           write_genotype)

# --- synthetic germlines ------------------------------------------------------------------------
# One gene, three alleles, arranged so that a read spanning germline 1..60 can separate `*02` from
# the others but CANNOT separate `*01` from `*03`. That is the partial-separability shape, and it
# is the reason this module exists in the form it does.

_rng = random.Random(7)
_BASE = "".join(_rng.choice("ACGT") for _ in range(120))
_FLIP = {"A": "C", "C": "A", "G": "T", "T": "G"}


def _mut(s: str, i: int) -> str:
    return s[:i] + _FLIP[s[i]] + s[i + 1:]


GERMLINES = {
    "G*01": _BASE,
    "G*02": _mut(_BASE, 30),        # differs INSIDE a 1..60 span -- separable
    "G*03": _mut(_BASE, 90),        # differs only OUTSIDE it -- NOT separable from *01
    "H*01": "".join(_rng.choice("ACGT") for _ in range(120)),   # a one-allele gene
}


@pytest.fixture(autouse=True)
def germlines(monkeypatch):
    import arda.genotype as gt
    monkeypatch.setattr(gt, "segment_germlines", lambda organism, segment="v": dict(GERMLINES))


def _airr(tmp_path, reads, name="in.tsv"):
    """``reads`` is ``(v_call, germline_end, junction)`` per row; span always starts at 1."""
    path = tmp_path / name
    pl.DataFrame({
        "sequence_id": [f"r{i}" for i in range(len(reads))],
        "locus": ["TRB"] * len(reads),
        "v_call": [c for c, _e, _j in reads],
        "j_call": ["TRBJ1-1*01"] * len(reads),
        "junction": [j for _c, _e, j in reads],
        "v_germline_start": ["1"] * len(reads),
        "v_germline_end": [str(e) for _c, e, _j in reads],
        "v_mutations": [""] * len(reads),
        "v_anchor_nt": ["120"] * len(reads),
    }).write_csv(path, separator="\t", quote_style="never")
    return path


def _calls(rows):
    """``{gene: (alleles, note)}`` from ``infer_genotype`` rows."""
    out: dict[str, tuple[list[str], str]] = {}
    for _locus, gene, allele, _c, _r, _gc, _bf, note in rows:
        alleles, _ = out.setdefault(gene, ([], note))
        if allele:
            alleles.append(allele)
    return {g: (sorted(a), n) for g, (a, n) in out.items()}


# --- the regression test this whole design exists for -------------------------------------------

def test_an_ambiguous_clonotype_constrains_without_choosing(tmp_path):
    """Never: dropping ambiguous votes is NOT the conservative choice -- it is the wrong answer.

    The donor is ``G*01/G*02`` and the reads span 1..60, where ``*02`` is separable and ``*01`` is
    not (it expands to ``{*01,*03}``). Counting only UNAMBIGUOUS votes sees 100 % ``*02``, calls
    the gene homozygous ``*02``, and then restriction empties every ``*01`` read -- systematically,
    reproducibly, for every donor with this genotype.

    Requiring the chosen set to also HIT the ``{*01,*03}`` votes is what rules ``{*02}`` out. Here
    ``*01`` and ``*03`` are themselves inseparable at this span, so the honest answer is that the
    gene cannot be resolved -- **not** a confident homozygous call.
    """
    reads = [("G*02", 60, f"c{i}") for i in range(10)] + [("G*01", 60, f"d{i}") for i in range(10)]
    rows, _ = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    alleles, note = _calls(rows)["G"]
    assert alleles != ["G*02"], "the homozygous-*02 answer is the bug"
    assert note in ("undetermined", "unexplained"), note


def test_the_same_shape_resolves_once_the_reads_are_long_enough(tmp_path):
    """The positive control for the test above: at the full span every allele separates."""
    reads = [("G*02", 120, f"c{i}") for i in range(10)] + [("G*01", 120, f"d{i}") for i in range(10)]
    rows, report = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    assert _calls(rows)["G"] == (["G*01", "G*02"], "ok")
    assert report["het"] == 1 and report["hom"] == 0


# --- the ordinary calls -------------------------------------------------------------------------

def test_a_homozygous_gene_is_called_with_one_allele(tmp_path):
    reads = [("G*02", 120, f"c{i}") for i in range(12)]
    rows, report = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    assert _calls(rows)["G"] == (["G*02"], "ok")
    assert report["hom"] == 1


def test_a_gene_with_one_catalogued_allele_is_carried_without_being_inferred(tmp_path):
    reads = [("H*01", 120, f"c{i}") for i in range(3)]
    rows, _ = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    assert _calls(rows)["H"] == (["H*01"], "single_allele")


def test_a_gene_below_min_clonotypes_is_reported_not_guessed(tmp_path):
    """Never: reported WITH ITS REASON, not dropped -- a gene missing from the table is
    indistinguishable from one that was never in the reference."""
    reads = [("G*02", 120, f"c{i}") for i in range(3)]
    rows, _ = infer_genotype(_airr(tmp_path, reads), min_clonotypes=10)
    assert _calls(rows)["G"] == ([], "low_support")


def test_one_expanded_clone_cannot_decide_the_gene(tmp_path):
    """Never: the vote unit is the CLONOTYPE, not the read.

    1,000 reads of ONE clone on ``*02`` against 10 separate clones on ``*01``. Counting reads makes
    ``*02`` 0.99 of the library and the gene is called homozygous ``*02`` -- off a single molecule.
    Counting clonotypes sees 1 against 10. This is the role TIgGER's ``j_max`` filter plays;
    ``stats.allele_candidate`` has no equivalent and does not weight by ``duplicate_count`` either.
    """
    reads = [("G*02", 120, "one-big-clone")] * 1000 + [("G*01", 120, f"d{i}") for i in range(10)]
    rows, report = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    alleles, _note = _calls(rows)["G"]
    assert "G*02" not in alleles
    assert report["clonotypes"] == 11, "1,010 reads collapse to 11 clonotypes"


def test_a_mutated_read_casts_no_vote(tmp_path):
    """Never: ``candidates`` slices the GERMLINE of the called allele -- the read's own bases never
    enter it -- so a read with a sequencing error at a discriminating base votes exactly as
    confidently as one that matches. The unmutated filter is what ties the vote to the read."""
    path = _airr(tmp_path, [("G*02", 120, f"c{i}") for i in range(12)])
    df = pl.read_csv(path, separator="\t", quote_char=None, infer_schema_length=0)
    df.with_columns(v_mutations=pl.lit("G45A")).write_csv(path, separator="\t", quote_style="never")
    rows, report = infer_genotype(path, min_clonotypes=5)
    assert report["clonotypes"] == 0 and not rows

    rows, report = infer_genotype(path, min_clonotypes=5, unmutated_only=False)
    assert report["clonotypes"] == 12


def test_the_unmutated_guard_raises_rather_than_silently_skipping(tmp_path):
    path = _airr(tmp_path, [("G*02", 120, "c0")])
    df = pl.read_csv(path, separator="\t", quote_char=None, infer_schema_length=0)
    df.drop("v_mutations").write_csv(path, separator="\t", quote_style="never")
    with pytest.raises(ValueError, match="v_mutations"):
        infer_genotype(path)


def test_output_does_not_depend_on_row_order(tmp_path):
    """Determinism is a requirement: polars `group_by` is a multithreaded hash aggregation."""
    reads = [("G*02", 120, f"c{i}") for i in range(10)] + [("G*01", 120, f"d{i}") for i in range(10)]
    a, _ = infer_genotype(_airr(tmp_path, reads, "a.tsv"), min_clonotypes=5)
    b, _ = infer_genotype(_airr(tmp_path, list(reversed(reads)), "b.tsv"), min_clonotypes=5)
    assert a == b


# --- restriction --------------------------------------------------------------------------------

_GENOTYPE = {"G": ("G*01", "G*02")}


def test_restriction_narrows_to_the_carried_allele():
    assert restrict(("G*01", "G*03"), "G*01", _GENOTYPE) == "G*01"


def test_a_refusal_to_expand_is_not_a_contradiction():
    """Never: ``None`` means the tie machinery declined -- short span, unknown allele, runaway list.
    Reading that as "contradicts the genotype" would delete the call of every short read."""
    assert restrict(None, "G*03", _GENOTYPE) == "G*03"


def test_a_contradicted_read_gets_an_empty_call_not_a_guess():
    assert restrict(("G*03",), "G*03", _GENOTYPE) == ""


def test_a_gene_that_was_not_genotyped_passes_through():
    assert restrict(("H*01",), "H*01", _GENOTYPE) == "H*01"


def test_restrict_airr_adds_a_column_and_leaves_v_call_alone(tmp_path):
    src = _airr(tmp_path, [("G*01", 60, "a"), ("G*02", 120, "b"), ("G*03", 120, "c")])
    from arda.genotype import restrict_airr
    out = tmp_path / "out.tsv"
    report = restrict_airr(src, out, _GENOTYPE, echo=lambda _m: None)
    before = pl.read_csv(src, separator="\t", quote_char=None, infer_schema_length=0)
    after = pl.read_csv(out, separator="\t", quote_char=None, infer_schema_length=0)
    assert after["v_call"].to_list() == before["v_call"].to_list()
    assert [x or "" for x in after["v_call_genotyped"].to_list()] == ["G*01", "G*02", ""]
    assert report["contradicted"] == 1


# --- the file -------------------------------------------------------------------------------------

def test_the_simplest_user_genotype_is_a_column_of_allele_names(tmp_path):
    p = tmp_path / "user.tsv"
    p.write_text("allele\nG*01\nG*02\n")
    assert read_genotype(p, organism="x") == {"G": ("G*01", "G*02")}


def test_an_allele_the_reference_cannot_call_raises(tmp_path):
    """Never: 884 human V alleles are functional in IMGT and 801 reach a scaffold, so naming one of
    the others is an easy mistake that would otherwise restrict nothing, silently."""
    p = tmp_path / "user.tsv"
    p.write_text("allele\nG*01\nG*99\n")
    with pytest.raises(ValueError, match="G\\*99"):
        read_genotype(p, organism="x")


def test_uncalled_genes_round_trip_as_absent_not_as_carrying_nothing(tmp_path):
    rows = [("TRB", "G", "G*01", "8", "19", "10", "12.30", "ok"),
            ("TRB", "K", "", "", "", "4", "0.00", "low_support")]
    out = tmp_path / "g.tsv"
    assert write_genotype(rows, out, params={"organism": "x", "alleles": len(GERMLINES)}) == 2
    assert out.read_text().startswith("# organism: x\n")
    got = read_genotype(out, organism=None)
    assert got == {"G": ("G*01",)}, "an uncalled gene must be ABSENT, so restriction leaves it be"
    assert restrict(("K*07",), "K*07", got) == "K*07"


def test_the_written_columns_are_the_declared_ones(tmp_path):
    out = tmp_path / "g.tsv"
    write_genotype([("TRB", "G", "G*01", "8", "19", "10", "12.30", "ok")], out)
    header = out.read_text().splitlines()[0].split("\t")
    assert header == list(GENOTYPE_COLUMNS)


def test_gene_of():
    assert gene_of("TRBV19*01") == "TRBV19"
    assert gene_of("TRBV19") == "TRBV19"


def test_a_tie_list_spanning_an_ungenotyped_gene_keeps_that_gene_s_allele():
    """Never: an allele is ruled out only by ITS OWN gene's genotype.

    Tie lists routinely span genes -- ``resolve-ties`` exists because IGLV2-14 and IGLV2-23 are
    indistinguishable over 70 nt. Asking "is this allele carried" gene-blind empties every read
    whose tie list merely brushes a genotyped gene. Measured on the 453-read fixture with a
    one-gene genotype: 79 rows came back empty, every one of them a read of some other gene.
    """
    # G is genotyped (*01,*02); H is not. A tie list holding the uncarried G*03 and H*01 must
    # keep H*01 rather than collapsing to nothing.
    assert restrict(("G*03", "H*01"), "G*03", _GENOTYPE) == "H*01"
    # ...and a tie list entirely inside the genotyped gene still contradicts.
    assert restrict(("G*03",), "G*03", _GENOTYPE) == ""


def test_a_row_with_no_v_call_is_not_counted_as_contradicted(tmp_path):
    """Never: a metric with no input is OMITTED, never 0.

    A J-only read has no V call, so a restriction has nothing to assess -- scoring it as a
    contradiction made a correct restriction look like it had rejected a sixth of the fixture
    (74 of 453 rows).
    """
    from arda.genotype import restrict_airr

    src = _airr(tmp_path, [("G*01", 120, "a"), ("", 120, "b"), ("G*03", 120, "c")])
    report = restrict_airr(src, tmp_path / "out.tsv", _GENOTYPE, echo=lambda _m: None)
    assert report == {"rows": 3, "no_call": 1, "narrowed": 0, "contradicted": 1,
                      "unchanged": 1, "genes": 1}


# --- the likelihood test ------------------------------------------------------------------------

def test_a_single_observation_is_undetermined_not_a_confident_call(tmp_path):
    """Never: a ratio of counts is not a confidence, and this is the shape that proved it.

    Real ``TRBV20-1`` on a TRB amplicon: **2,544 clonotypes, exactly 1** of which could name an
    allele. The parsimony rule this replaced reported ``explained = 1.0000, note = ok`` and called
    the gene homozygous off that one observation -- identically to ``TRBV11-2``, which it called
    off **754 of 757**. A likelihood ratio separates them; a coverage fraction cannot.
    """
    # 1 clonotype can discriminate; 40 cannot (their span covers no difference).
    reads = [("G*01", 120, "informative")] + [("G*01", 60, f"blind{i}") for i in range(40)]
    rows, _ = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    alleles, note = _calls(rows)["G"]
    assert alleles == [] and note in ("undetermined", "low_support"), note


def test_overwhelming_evidence_gets_a_large_bayes_factor(tmp_path):
    """The other side of the same claim: real depth must produce a real number."""
    reads = [("G*01", 120, f"c{i}") for i in range(200)]
    rows, _ = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    row = next(r for r in rows if r[1] == "G")
    assert row[2] == "G*01" and row[7] == "ok"
    assert float(row[6]) > 50, f"log10 BF {row[6]} is too small for 200 concordant clonotypes"


def test_read_depth_is_reported_beside_the_clonotype_count(tmp_path):
    """Never: read coverage and clonotype coverage are different numbers and both are wanted.

    The clonotype is the independent observation -- a junction is de facto a UMI, so one junction
    is one rearrangement however many reads carry it -- but depth is what says whether that
    clonotype's own allele assignment can be trusted. Reporting only one of the two hides half the
    question.
    """
    reads = ([("G*01", 120, "a")] * 5 + [("G*01", 120, "b")] * 3
             + [("G*01", 120, f"c{i}") for i in range(10)])
    rows, report = infer_genotype(_airr(tmp_path, reads), min_clonotypes=5)
    row = next(r for r in rows if r[1] == "G")
    assert row[3] == "12", "12 clonotypes"
    assert row[4] == "18", "18 reads behind them"
    assert report["clonotypes"] == 12 and report["reads"] == 18


def test_the_error_rate_is_measured_from_within_clonotype_disagreement(tmp_path):
    """Never: reads of ONE junction are one rearrangement and so one allele -- a read naming
    another is hypermutation or sequencing error by construction, never allelic variation. That
    makes the miscall rate measurable on the library instead of a constant somebody picked.

    ⚠ And it must be measured over ALL reads, not the germline-exact subset the assignment uses:
    those were selected for carrying no mismatch, so they never disagree and the estimate collapses
    silently onto its own floor. On a real TRB amplicon that was 0 discordant reads of 77,345.
    """
    from arda.genotype import MIN_ERROR_RATE, error_rate

    clean = _clonotype_records(tmp_path, [("G*01", 120, f"c{i}") for i in range(10)])
    assert error_rate(clean) == MIN_ERROR_RATE, "no disagreement -> the floor, not zero"

    # One junction whose reads disagree about the allele: 9 say *01, 1 says *02.
    mixed = _clonotype_records(
        tmp_path, [("G*01", 120, "same")] * 9 + [("G*02", 120, "same")])
    assert error_rate(mixed) == pytest.approx(0.1)


def _clonotype_records(tmp_path, reads):
    from arda.annotate.airr_out import read_airr
    from arda.annotate.ties import TieResolver
    from arda.genotype import _clonotypes
    df = read_airr(_airr(tmp_path, reads, f"e{len(reads)}.tsv"))
    return _clonotypes(df, TieResolver(dict(GERMLINES)), unmutated_only=True, scope="framework")
