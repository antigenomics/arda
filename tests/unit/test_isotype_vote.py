"""A clonotype's isotype is decided by ONE VOTE PER FRAGMENT — per molecule, not per row.

⛔ Report the isotype CLASS (IGHG/IGHM/IGHA), never the subclass: IGHG1-4 are ~95 % identical, so
the top *gene* ties 26.7 % of the time while the top *class* never does.

The vote is a two-level dedup and both levels are load-bearing:

* ``_dominant_ccall`` deduplicates the clonotype's reads down to fragments, so a fragment whose two
  mates were both ASSIGNED does not vote twice;
* ...and each fragment's own calls collapse to one vote, because ``frag_iso`` holds one entry per
  ROW carrying a ``c_class``. Without that second level a fragment whose two mates BOTH carry a
  constant call voted twice, and an assembly-rescued fragment voted a third time, since
  ``assemble`` copies the member read's ``c_class`` onto the rescued row while the Stage-1 row
  still carries it. A one-fragment minority could then outvote a two-fragment majority.
"""

from __future__ import annotations

import polars as pl

from arda.rnaseq.correct import correct_airr

_JN = "TGTGCGAAAGGGGCCCTTCAGAAAACATTACGTTTGGGGGAGTCTATACCCCTAAATCCTTTTGATGTCTGG"


def _row(sid, **kw):
    r = {"sequence_id": sid, "sequence": _JN, "junction": _JN, "junction_aa": "CAKGALQKW",
         "v_call": "IGHV3-23*01", "j_call": "IGHJ3*01", "locus": "IGH", "c_class": ""}
    r.update(kw)
    return r


def _isotype(tmp_path, rows, name):
    src = tmp_path / f"{name}.tsv"
    pl.DataFrame(rows).write_csv(src, separator="\t")
    out = tmp_path / f"{name}.out"
    correct_airr(src, out, map_d=False)
    got = pl.read_csv(out, separator="\t", infer_schema_length=0).to_dicts()
    assert len(got) == 1, got
    return got[0]["c_call"]


def test_a_fragment_with_two_constant_mates_votes_once(tmp_path):
    """TWO fragments say IGHM; ONE fragment says IGHG but carries the call on BOTH its mates.

    Per fragment that is 2-1 for IGHM. Counting rows it is 2-2, and the tie then falls to whichever
    class the row order happened to put first.
    """
    rows = [
        _row("fragM1/1", c_class="IGHM"),
        _row("fragM2/1", c_class="IGHM"),
        _row("fragG/1", c_class="IGHG"),
        _row("fragG/2", c_class="IGHG"),      # same MOLECULE, second mate
    ]
    assert _isotype(tmp_path, rows, "mates") == "IGHM"


def test_a_duplicated_constant_call_on_one_fragment_cannot_outvote_two_fragments(tmp_path):
    """The assembly shape: the rescued row copies the member read's own ``c_class``, so the same
    fragment's call appears on two rows of the concatenated frame."""
    rows = [
        _row("fragA/1", c_class="IGHA"),
        _row("fragA2/1", c_class="IGHA"),
        _row("fragG/1", c_class="IGHG"),
        _row("fragG/1", c_class="IGHG"),      # duplicate row for the SAME read
        _row("fragG/1", c_class="IGHG"),
    ]
    assert _isotype(tmp_path, rows, "dup") == "IGHA"


def test_a_resolved_class_beats_the_generic_one(tmp_path):
    """``isotype_class`` emits the generic ``IGHC`` only on cross-class ambiguity, so a handful of
    ambiguous reads must never outvote a resolved call."""
    rows = [_row(f"amb{i}/1", c_class="IGHC") for i in range(5)] + [_row("res/1", c_class="IGHG")]
    assert _isotype(tmp_path, rows, "generic") == "IGHG"


def test_the_vote_does_not_depend_on_row_order(tmp_path):
    """⛔ A tie must break lexicographically, never on encounter order -- row order here comes from
    a threaded mmseqs search, so an order-dependent tie is not reproducible."""
    rows = [_row("f1/1", c_class="IGHG"), _row("f2/1", c_class="IGHM")]
    a = _isotype(tmp_path, rows, "ord_a")
    b = _isotype(tmp_path, list(reversed(rows)), "ord_b")
    assert a == b == "IGHG"          # tie -> lexicographic winner
