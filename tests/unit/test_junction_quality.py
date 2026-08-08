"""Quality-aware error correction: ``map --junction-quality`` and ``correct --min-junction-q``.

The abundance error model in :func:`arda.rnaseq.correct._parents` cannot tell a sequencing miscall
from a real low-frequency variant -- both are "a rare neighbour of an abundant clonotype" -- and
``--error-rate`` only trades them off globally. Phred *can*, because it is a different measurement:
measured at the mismatching base over 302,172 real MIGEC windows, the two published spike-in
variants sit at median Q 34-35 and the 1-substitution error cloud around them at median Q 24
(54 % below Q30, against 5 % for the parent clone's own bases).

Two halves, tested here:

* Stage 1 is the only place the FASTQ quality is still in hand, so ``map --junction-quality``
  carries it forward as a ``junction_quality`` column aligned base-for-base with ``junction``.
  ⛔ The failure mode that no downstream check can catch is a quality string of the RIGHT LENGTH
  taken from the wrong strand or offset, so the tests below assert the *content*, on both strands.
* Stage 2 gates each read on the Phred of the bases that differ from its putative parent.
  Matching bases are not evidence and are never looked at.
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.annotate.airr_out import airr_header, format_rows, _format_rows_py, _format_rows_cpp
from arda.annotate.transfer import AIRR_COLUMNS
from arda.rnaseq.correct import EC_MODES, correct_airr
from arda.rnaseq.map import JUNCTION_QUALITY, junction_quality, map_rnaseq, read_pairs
from arda.refbuild.translate import reverse_complement

# ---------------------------------------------------------------------------------------------
# A junction placed inside a synthetic read. `junction` = CDR3 flanked by the two anchor codons,
# so cdr3_start = junction_start + 3 and cdr3_end = junction_end - 3 (1-based, closed).
_LEAD = "ACGTACGTAC"                                   # 10 nt before the junction
_JUNC = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"               # 30 nt, C...F, in frame
_TAIL = "GGGACCAGG"
_SEQ = _LEAD + _JUNC + _TAIL
_QUAL = "I" * len(_LEAD) + "".join("I5I"[i % 3] for i in range(len(_JUNC))) + "I" * len(_TAIL)


def _rec(**kw) -> dict:
    rec = {"sequence": _SEQ, "junction": _JUNC, "rev_comp": "F",
           "cdr3_start": len(_LEAD) + 4, "cdr3_end": len(_LEAD) + len(_JUNC) - 3}
    rec.update(kw)
    return rec


# ------------------------------------------------------------------ junction_quality extraction
def test_forward_read_quality_is_the_junction_slice():
    q = junction_quality(_rec(), _QUAL)
    assert q == _QUAL[len(_LEAD):len(_LEAD) + len(_JUNC)]
    assert len(q) == len(_JUNC)


def test_reverse_complement_read_reverses_the_quality_it_does_not_complement_it():
    """⛔ The one bug a length check cannot catch.

    For ``rev_comp == "T"`` the stored ``sequence`` is the read as submitted and everything else
    in the record is on the coding strand, so the quality runs the other way. Reversing it is
    load-bearing; not reversing yields a same-length string of the WRONG bases' qualities.
    """
    rec = _rec(sequence=reverse_complement(_SEQ), rev_comp="T")
    submitted_qual = _QUAL[::-1]                       # the quality as it sits in the FASTQ
    q = junction_quality(rec, submitted_qual)
    assert q == _QUAL[len(_LEAD):len(_LEAD) + len(_JUNC)]
    # And the un-reversed answer is a different string of the same length -- i.e. the assertion
    # above is not satisfied by accident.
    naive = submitted_qual[len(_LEAD):len(_LEAD) + len(_JUNC)]
    assert len(naive) == len(q) and naive != q


def test_no_quality_available_returns_empty():
    assert junction_quality(_rec(), "") == ""


def test_no_junction_returns_empty():
    assert junction_quality(_rec(junction=""), _QUAL) == ""


def test_missing_coordinates_fall_back_to_a_search():
    rec = _rec()
    rec["cdr3_start"] = rec["cdr3_end"] = ""
    assert junction_quality(rec, _QUAL) == _QUAL[len(_LEAD):len(_LEAD) + len(_JUNC)]


def test_wrong_coordinates_do_not_produce_a_misaligned_slice():
    """Coordinates that do not reproduce the junction are discarded, not trusted."""
    rec = _rec(cdr3_start=1, cdr3_end=len(_JUNC) - 6)
    assert junction_quality(rec, _QUAL) == _QUAL[len(_LEAD):len(_LEAD) + len(_JUNC)]


def test_junction_absent_from_the_sequence_returns_empty():
    assert junction_quality(_rec(sequence="AAAAAAAAAA"), "IIIIIIIIII") == ""


def test_quality_of_a_different_length_than_the_sequence_returns_empty():
    assert junction_quality(_rec(), _QUAL[:-5]) == ""


# ------------------------------------------------------------------------- the extra AIRR column
def test_header_and_rows_are_unchanged_without_the_flag():
    assert airr_header() == "\t".join(AIRR_COLUMNS)
    rec = dict.fromkeys(AIRR_COLUMNS, "")
    assert len(format_rows([rec]).rstrip("\n").split("\t")) == len(AIRR_COLUMNS)


def test_extra_column_is_appended_last_and_empty_when_the_key_is_absent():
    extra = (JUNCTION_QUALITY,)
    assert airr_header(extra).split("\t")[-1] == JUNCTION_QUALITY
    rec = dict.fromkeys(AIRR_COLUMNS, "")
    fields = format_rows([rec], extra).rstrip("\n").split("\t")
    assert len(fields) == len(AIRR_COLUMNS) + 1 and fields[-1] == ""
    rec[JUNCTION_QUALITY] = "IIII"
    assert format_rows([rec], extra).rstrip("\n").split("\t")[-1] == "IIII"


@pytest.mark.skipif(_format_rows_cpp is None, reason="C++ extension not built")
def test_python_and_cpp_agree_on_the_extra_column():
    cols = tuple(AIRR_COLUMNS) + (JUNCTION_QUALITY,)
    recs = [dict.fromkeys(AIRR_COLUMNS, ""), {**dict.fromkeys(AIRR_COLUMNS, ""),
                                              JUNCTION_QUALITY: "IIII"}]
    assert _format_rows_cpp(recs, cols) == _format_rows_py(recs, cols)


# ------------------------------------------------------------------------------ reading quality
def _fastq(path, records) -> None:
    path.write_text("".join(f"@{i}\n{s}\n+\n{q}\n" for i, s, q in records))


def test_read_pairs_with_qual_yields_the_quality_single_end(tmp_path):
    f = tmp_path / "r.fq"
    _fastq(f, [("a", "ACGT", "IIII"), ("b", "TTTT", "!!!!")])
    assert list(read_pairs(f, with_qual=True)) == [("a", "ACGT", "IIII"), ("b", "TTTT", "!!!!")]
    assert list(read_pairs(f)) == [("a", "ACGT"), ("b", "TTTT")]


def test_read_pairs_with_qual_yields_the_quality_paired(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _fastq(r1, [("a", "ACGT", "IIII")])
    _fastq(r2, [("a", "TTTT", "!!!!")])
    assert list(read_pairs(r1, r2, with_qual=True)) == [
        ("a/1", "ACGT", "IIII"), ("a/2", "TTTT", "!!!!")]


def test_junction_quality_refuses_reconstruct(tmp_path):
    """A merged fragment's bases come from two reads, so no input quality describes it."""
    f = tmp_path / "r.fq"
    _fastq(f, [("a", "ACGT", "IIII")])
    with pytest.raises(ValueError, match="reconstruct"):
        map_rnaseq(f, tmp_path / "o.tsv", r2=f, reconstruct=True, with_junction_quality=True)


# ------------------------------------------------------------------------------- the Stage-2 gate
# One abundant parent and two 1-substitution children of equal abundance. They differ from the
# parent at the SAME position; the only thing separating them is the Phred at that base.
_PARENT = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"
_CHILD = "TGTGCCAGCAGCTTAGACGGGACAGGGTTT"       # differs at the last base (position 29)
_MISPOS = 29


def _q(at_mispos: str, rest: str = "I", n: int = len(_PARENT)) -> str:
    return "".join(at_mispos if i == _MISPOS else rest for i in range(n))


def _airr(tmp_path, rows, name="in.tsv"):
    p = tmp_path / name
    pl.DataFrame(rows).write_csv(p, separator="\t")
    return p


def _reads(junction, n, qual, *, prefix):
    return [{"sequence_id": f"{prefix}{k}", "sequence": junction, "junction": junction,
             "junction_aa": "CASSLDGTGF",
             "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01", "locus": "TRB",
             "junction_quality": qual} for k in range(n)]


def _clonotypes(path, out, **kw):
    correct_airr(path, out, map_d=False, **kw)
    return {r["junction"]: r["duplicate_count"]
            for r in pl.read_csv(out, separator="\t").to_dicts()}


def test_gate_raises_when_the_column_is_absent(tmp_path):
    """⛔ Raise, never degrade: a silently unapplied gate looks exactly like a gate that passed."""
    rows = _reads(_PARENT, 5, "I" * len(_PARENT), prefix="p")
    for r in rows:
        del r["junction_quality"]
    with pytest.raises(ValueError, match="junction_quality"):
        correct_airr(_airr(tmp_path, rows), tmp_path / "o.tsv", map_d=False, min_junction_q=20)


def test_low_quality_child_is_dropped_and_high_quality_child_survives(tmp_path):
    # `error_rate` low enough that the abundance model keeps both children: the gate is then the
    # ONLY thing that can separate them, which is the point of the test.
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, _q("I"), prefix="hi"))
    both = _clonotypes(_airr(tmp_path, rows, "hi.tsv"), tmp_path / "hi.out",
                       error_rate=1e-6, min_junction_q=20)
    assert both.get(_CHILD) == 20                       # high-Q discriminating base: kept

    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, _q("5"), prefix="lo"))  # '5' = Q20-... see below
    low = _clonotypes(_airr(tmp_path, rows, "lo.tsv"), tmp_path / "lo.out",
                      error_rate=1e-6, min_junction_q=25)
    assert _CHILD not in low                            # low-Q discriminating base: gone
    assert low[_PARENT] == 420                          # ...and its reads went to the parent


def test_the_threshold_boundary_is_inclusive(tmp_path):
    """Exactly at the threshold is KEPT; one below is dropped."""
    at = chr(33 + 20)                                   # Phred 20
    below = chr(33 + 19)
    for qchar, expect in ((at, 20), (below, None)):
        rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
                + _reads(_CHILD, 20, _q(qchar), prefix="c"))
        cl = _clonotypes(_airr(tmp_path, rows, f"b{ord(qchar)}.tsv"),
                         tmp_path / f"b{ord(qchar)}.out",
                         error_rate=1e-6, min_junction_q=20)
        assert cl.get(_CHILD) == expect


def test_only_the_mismatching_bases_are_looked_at(tmp_path):
    """Every base except the discriminating one is Q0; the read survives on that one base."""
    qual = "".join(chr(33 + 40) if i == _MISPOS else "!" for i in range(len(_PARENT)))
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, qual, prefix="c"))
    cl = _clonotypes(_airr(tmp_path, rows, "m.tsv"), tmp_path / "m.out",
                     error_rate=1e-6, min_junction_q=30)
    assert cl.get(_CHILD) == 20


def test_all_zero_quality_drops_every_gated_read_but_never_the_parent(tmp_path):
    """A clonotype with no more-abundant neighbour has no hypothesis to test, so it is not gated
    -- even at Q0 across the board. Otherwise an all-``!`` library would report nothing at all."""
    rows = (_reads(_PARENT, 400, "!" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, "!" * len(_PARENT), prefix="c"))
    cl = _clonotypes(_airr(tmp_path, rows, "z.tsv"), tmp_path / "z.out",
                     error_rate=1e-6, min_junction_q=20)
    assert cl == {_PARENT: 420}


def test_a_clonotype_with_no_more_abundant_neighbour_is_untouched(tmp_path):
    """Zero mismatching positions to test -> the gate must be a no-op, at any threshold."""
    lone = "TGTGCCAGCAGTAAACCCGGGTTTAAACCCTTC"
    rows = _reads(lone, 30, "!" * len(lone), prefix="x")
    cl = _clonotypes(_airr(tmp_path, rows, "l.tsv"), tmp_path / "l.out",
                     error_rate=1e-6, min_junction_q=40)
    assert cl == {lone: 30}


def test_a_read_without_a_quality_string_is_kept(tmp_path):
    """Absent evidence is not evidence of error."""
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, "", prefix="c"))
    cl = _clonotypes(_airr(tmp_path, rows, "n.tsv"), tmp_path / "n.out",
                     error_rate=1e-6, min_junction_q=40)
    assert cl.get(_CHILD) == 20


def test_the_report_counts_what_the_gate_removed(tmp_path):
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, _q("!"), prefix="c"))
    rep = correct_airr(_airr(tmp_path, rows), tmp_path / "o.tsv", map_d=False,
                       error_rate=1e-6, min_junction_q=20)
    assert rep.reads_low_quality == 20
    assert rep.clonotypes_low_quality == 1


# ------------------------------------------------------------------------------ the mode selector
def test_fast_mode_is_the_shipped_default(tmp_path):
    """``--ec-mode fast`` must be a no-op: the default output cannot move."""
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, _q("!"), prefix="c"))
    src = _airr(tmp_path, rows)
    correct_airr(src, tmp_path / "a.tsv", map_d=False, error_rate=1e-6)
    correct_airr(src, tmp_path / "b.tsv", map_d=False, error_rate=1e-6, ec_mode="fast")
    assert (tmp_path / "a.tsv").read_text() == (tmp_path / "b.tsv").read_text()
    assert EC_MODES["fast"]["min_junction_q"] == 0


def test_accurate_mode_applies_the_quality_gate(tmp_path):
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, _q("!"), prefix="c"))
    cl = _clonotypes(_airr(tmp_path, rows), tmp_path / "o.tsv",
                     error_rate=1e-6, ec_mode="accurate")
    assert _CHILD not in cl


def test_an_explicit_knob_overrides_the_mode(tmp_path):
    rows = (_reads(_PARENT, 400, "I" * len(_PARENT), prefix="p")
            + _reads(_CHILD, 20, _q("!"), prefix="c"))
    cl = _clonotypes(_airr(tmp_path, rows), tmp_path / "o.tsv",
                     error_rate=1e-6, ec_mode="accurate", min_junction_q=0)
    assert cl.get(_CHILD) == 20


def test_unknown_mode_raises(tmp_path):
    rows = _reads(_PARENT, 5, "I" * len(_PARENT), prefix="p")
    with pytest.raises(ValueError, match="ec_mode"):
        correct_airr(_airr(tmp_path, rows), tmp_path / "o.tsv", map_d=False, ec_mode="turbo")


def test_ec_mode_on_rnaseq_run_implies_the_stage1_column_it_needs():
    """⛔ `run` does Stage 1 and Stage 2 in ONE call, so a user cannot turn on
    `map --junction-quality` by hand. If asking for the gate did not imply producing its input,
    `--ec-mode accurate` would silently do nothing there -- the exact class of failure this
    codebase keeps hitting (a flag that is accepted, changes nothing, and reports success).
    """
    import inspect

    from arda.rnaseq import pipeline

    src = inspect.getsource(pipeline.run)
    assert "with_junction_quality=" in src, "run() must drive the Stage-1 column"
    assert 'ec_mode != "fast"' in src and "min_junction_q is not None" in src, (
        "asking for the gate, by either route, must imply producing its input")


def test_ec_mode_and_min_junction_q_reach_correct_through_the_pipeline():
    """Both stages of `run` must forward the knobs, or the preset stops at the CLI boundary."""
    import inspect

    from arda.rnaseq import pipeline

    for fn in (pipeline.run, pipeline.finish):
        sig = inspect.signature(fn)
        assert "ec_mode" in sig.parameters, fn.__name__
        assert "min_junction_q" in sig.parameters, fn.__name__
    assert "ec_mode=ec_mode" in inspect.getsource(pipeline.finish)
