"""The parts of :mod:`arda.cdr3fix` the unit suite never touched.

Four public symbols (``markup_batch``, ``to_frame``, ``format_report``, ``MARKUP_COLUMNS``)
had no test at all, and three of the seven VDJdb fix types were only reachable through the
realworld fixture. Neither ``Cdr3Markup.explain`` nor ``Cdr3Error.__str__`` was ever called,
though both are what ``arda markup --report`` prints.
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.cdr3fix import (MARKUP_COLUMNS, format_report, load_anchors, markup_batch,
                          markup_cdr3, markup_records, resolve_allele, to_frame)

V, J = "TRBV9*01", "TRBJ2-2*01"
CLEAN = "CASSARSGELFF"


# --------------------------------------------------------------------------- batch / frame


def _frame(rows):
    return pl.DataFrame({"cdr3": [r[0] for r in rows], "v": [r[1] for r in rows],
                         "j": [r[2] for r in rows], "species": ["human"] * len(rows)})


def test_to_frame_of_nothing_still_has_the_schema():
    empty = to_frame([])
    assert empty.height == 0
    assert empty.columns == list(MARKUP_COLUMNS)


def test_markup_batch_returns_one_row_per_record_in_order():
    src = _frame([(CLEAN, V, J), ("CASSARSGELF", V, J)])
    out = markup_batch(src)
    assert out.columns == list(MARKUP_COLUMNS)
    assert out.height == 2
    assert out["cdr3"].to_list() == [CLEAN, "CASSARSGELF"]
    assert out["cdr3_repaired"].to_list() == [CLEAN, CLEAN]     # the terminal F restored
    assert out["good"].to_list() == [True, True]


def test_to_frame_serialises_errors_as_the_report_renders_them():
    recs = markup_records(_frame([("CASSARSGELF", V, J)]))
    row = to_frame(recs).to_dicts()[0]
    assert row["j_fix"] == "FixAdd" and row["n_errors"] == 1
    assert "missing 'F'" in row["errors"]
    assert row["v_end"] == 4 and row["j_start"] >= 0


def test_format_report_counts_and_lists():
    recs = markup_records(_frame([(CLEAN, V, J), ("CASSARSGELF", V, J), ("F", V, J)]))
    txt = format_report(recs, show_ok=True)
    assert "cdr3fix report: 3 records" in txt
    assert "correct (no fix needed) : 1" in txt
    assert "repaired                : 1" in txt
    assert "failed                  : 1" in txt
    assert "[OK]" in txt and "[FIXED]" in txt and "[FAILED]" in txt
    # show_ok=False hides the clean record's line but keeps it in the count
    assert "[OK]" not in format_report(recs, show_ok=False)


def test_explain_and_error_str_are_what_the_report_prints():
    mk = markup_cdr3("CASSARSGELF", V, J, "human")
    assert str(mk.errors[0]) == "J del@10 missing 'F' d=0"
    assert "CASSARSGELF" in mk.explain() and "CASSARSGELFF" in mk.explain()


# --------------------------------------------------------------------------- fix types


def test_v_side_missing_cys_is_added():
    """The mirror of the J-side missing-F case, which was the only one covered."""
    mk = markup_cdr3("ASSARSGELFF", V, J, "human")
    assert mk.cdr3_repaired == CLEAN
    assert mk.v_fix == "FixAdd" and mk.v_canonical is False
    assert [(e.side, e.kind, e.to) for e in mk.errors if e.applied] == [("V", "del", "C")]


@pytest.mark.parametrize("cdr3,side,fix", [
    ("XCASSARSGELFF", "V", "v_fix"),      # a residue before Cys104
    ("CASSARSGELFFW", "J", "j_fix"),      # a residue past Phe118
])
def test_extra_residue_at_an_anchor_is_trimmed(cdr3, side, fix):
    mk = markup_cdr3(cdr3, V, J, "human")
    assert mk.cdr3_repaired == CLEAN
    assert getattr(mk, fix) == "FixTrim"
    assert [(e.side, e.kind) for e in mk.errors if e.applied] == [(side, "ins")]


@pytest.mark.parametrize("cdr3,which", [("C", "j_fix"), ("F", "v_fix")])
def test_a_junction_with_nothing_to_align_fails_rather_than_guesses(cdr3, which):
    mk = markup_cdr3(cdr3, V, J, "human")
    assert getattr(mk, which) == "FailedNoAlignment"
    assert not mk.good
    assert mk.cdr3_repaired == cdr3, "a failed side must not rewrite the junction"


def test_max_replace_bounds_how_deep_a_repair_may_reach():
    """``max_replace`` is a *distance from the anchor*, so 0 still repairs the anchor itself.

    A substitution one residue in (dist=1) is refused at 0 and applied at 1. This is the knob
    that keeps arda from rewriting the V/N boundary as if it were a curation error.
    """
    off_by_one = "CXSSARSGELFF"                       # sub at index 1, dist 1 from Cys104
    assert markup_cdr3(off_by_one, V, J, "human", max_replace=0).cdr3_repaired == off_by_one
    assert markup_cdr3(off_by_one, V, J, "human", max_replace=1).cdr3_repaired == CLEAN

    reported = markup_cdr3(off_by_one, V, J, "human", max_replace=0)
    assert reported.v_fix == "NoFixNeeded", "a reported-only error leaves the fix type clean"
    assert reported.errors and not any(e.applied for e in reported.errors)

    # dist=0 is *at* the anchor and is repaired even at max_replace=0.
    assert markup_cdr3("CASSARSGELF", V, J, "human", max_replace=0).cdr3_repaired == CLEAN


# --------------------------------------------------------------------------- resolution


def test_allele_ladder_falls_back_to_the_first_allele_of_the_gene():
    """IGLV3-4 ships no *01, so the ladder lands on its lowest-numbered allele."""
    anchors = load_anchors("human")
    assert ("V", "IGLV3-4*01") not in anchors
    assert resolve_allele("IGLV3-4", "V", anchors) == "IGLV3-4*02"


def test_an_unknown_gene_resolves_to_nothing():
    assert resolve_allele("NOSUCHV1-1*01", "V", load_anchors("human")) == ""


def test_an_ambiguous_call_uses_its_first_allele():
    """`v_call` may be a comma list, as arda's own AIRR output emits. Take the first."""
    both = markup_cdr3(CLEAN, f"{V},TRBV20-1*01", f"{J},TRBJ1-4*01", "human")
    first = markup_cdr3(CLEAN, V, J, "human")
    assert (both.v_end, both.j_start) == (first.v_end, first.j_start)
    assert both.cdr3_repaired == first.cdr3_repaired
