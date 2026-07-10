"""The committed junction-markup example must still reproduce, byte for byte.

`examples/README.md` quotes this output. `examples/example.airr.tsv` went stale for four
release rounds — written with 49 AIRR columns while the schema grew to 83 — because nothing
checked it. These artifacts need no mmseqs, so the guard is a plain unit test.

If this fails after an intentional change, run ``python examples/regenerate.py`` and read the
diff before committing it: the diff *is* the behaviour change.
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

from arda.cdr3fix import format_report, markup_records, to_frame

EXAMPLES = Path(__file__).resolve().parent.parent.parent / "examples"


@pytest.fixture(scope="module")
def marked():
    src = pl.read_csv(EXAMPLES / "junctions.tsv", separator="\t", infer_schema_length=0)
    return markup_records(src, cdr3="cdr3", v="v", j="j", species="species", sequence_id="id")


def test_committed_markup_tsv_still_reproduces(marked):
    got = to_frame(marked)
    want = pl.read_csv(EXAMPLES / "junctions.markup.tsv", separator="\t", infer_schema_length=0)
    # The committed file carries the --d-posterior columns too; compare the markup columns.
    shared = [c for c in got.columns if c in want.columns]
    assert shared, "committed example lost every markup column"
    for col in shared:
        assert got[col].cast(pl.Utf8).fill_null("").to_list() == \
            want[col].cast(pl.Utf8).fill_null("").to_list(), f"column {col!r} drifted"


def test_committed_report_still_reproduces(marked):
    want = (EXAMPLES / "junctions.report.txt").read_text()
    assert format_report(marked, show_ok=True).rstrip("\n") == want.rstrip("\n")


def test_every_documented_repair_outcome_is_present(marked):
    """examples/README.md tabulates one record per outcome. Keep them all reachable."""
    by_id = {m.sequence_id: m for m in marked}
    assert set(by_id) == {"clean-trb", "v-anchor-sub", "j-anchor-missing-F",
                          "j-anchor-extra-G", "reported-not-repaired",
                          "flanking-fr3-trimmed", "bad-segment-refused"}

    assert (by_id["clean-trb"].v_fix, by_id["clean-trb"].j_fix) == ("NoFixNeeded", "NoFixNeeded")
    assert by_id["clean-trb"].good and not by_id["clean-trb"].errors

    m = by_id["v-anchor-sub"]                       # F -> C at Cys104
    assert m.cdr3_repaired == "CLVGPQGSSASKIIF" and m.v_fix == "FixReplace"

    m = by_id["j-anchor-missing-F"]                 # the Phe118 anchor restored
    assert m.cdr3_repaired == "CAIRDDKIIF" and m.j_fix == "FixAdd"

    m = by_id["j-anchor-extra-G"]                   # a residue past Phe118 trimmed
    assert m.cdr3_repaired == "CATSSPGLASDEQFF" and m.j_fix == "FixTrim"

    m = by_id["flanking-fr3-trimmed"]               # framework on BOTH flanks, both trimmed
    assert m.cdr3_repaired == "CASPGGIQYF"          # exactly what VDJdb's own fixer emits
    assert (m.v_fix, m.j_fix) == ("FixTrim", "FixTrim") and m.good
    assert {(e.side, e.frm) for e in m.errors if e.applied} == {("V", "YF"), ("J", "GAG")}


def test_every_repair_lands_on_a_canonical_junction(marked):
    """A repair exists to restore the anchors. `good` implies both are there."""
    for m in marked:
        if m.good:
            assert m.cdr3_repaired.startswith("C"), m.sequence_id
            assert m.cdr3_repaired.endswith(("F", "W")), m.sequence_id
        assert m.v_canonical == m.cdr3_repaired.startswith("C")
        assert m.j_canonical == m.cdr3_repaired.endswith(("F", "W"))


def test_the_deep_error_is_reported_and_not_repaired(marked):
    """The example that exists to show detection and repair are separate decisions."""
    m = next(x for x in marked if x.sequence_id == "reported-not-repaired")
    assert m.cdr3_repaired == "CASSSPLLSSDTQYF"
    deep = [e for e in m.errors if not e.applied]
    near = [e for e in m.errors if e.applied]
    assert len(deep) == 1 and deep[0].kind == "sub" and deep[0].dist == 6
    assert "reported, not repaired" in str(deep[0])
    assert len(near) == 1 and near[0].kind == "ins" and near[0].dist == 0
    # only the anchor-adjacent edit reached the output
    assert m.cdr3_repaired == "CASSSPLLSSDTQYFG"[:-1]


def test_a_bad_segment_is_flagged_never_guessed(marked):
    """TRAV1-2*02 is a real IMGT allele arda ships with no usable anchor."""
    m = next(x for x in marked if x.sequence_id == "bad-segment-refused")
    assert m.cdr3_repaired == "CAVRSMDSNYQLIW", "a refused record must come back untouched"
    assert not m.good and m.v_fix == "FailedBadSegment" and m.v_end == -1
    assert m.j_fix == "NoFixNeeded", "the J side is fine; only the V allele has no anchor"
    assert not any(e.applied for e in m.errors)


def test_applied_edits_exist_exactly_when_the_junction_changed(marked):
    """The invariant the ``applied`` flag is supposed to carry, on every example."""
    for m in marked:
        assert any(e.applied for e in m.errors) == (m.cdr3_repaired != m.cdr3), m.sequence_id


def test_cdr3fix_json_is_vdjdb_shaped(marked):
    blob = next(x for x in marked if x.sequence_id == "j-anchor-missing-F").to_cdr3fix()
    assert blob["cdr3"] == "CAIRDDKIIF" and blob["cdr3_old"] == "CAIRDDKII"
    assert blob["jFixType"] == "FixAdd" and blob["good"] is True
    # round-trips through JSON the way VDJdb's own column does
    assert json.loads(json.dumps(blob))["vFixType"] == "NoFixNeeded"
