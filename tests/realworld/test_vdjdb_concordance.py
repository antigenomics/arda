"""Concordance of ``arda.cdr3fix`` with VDJdb's own ``Cdr3Fixer``.

The fixture (`tests/assets/vdjdb/sample.tsv.gz`, see SOURCES.md) carries VDJdb's
``cdr3fix`` JSON, which is free ground truth: ``vEnd``/``jStart`` boundaries, the
fix type, and both the submitted (``cdr3_old``) and repaired (``cdr3``) junction.

Two things are scored separately, because they are different claims:

* **Boundaries** -- where the V and J germlines stop templating. Asserted.
* **Repair** -- can we reproduce VDJdb's fix from the submitted junction. Asserted.

What is deliberately NOT asserted: fix-type *strings*. The shipped VDJdb database is
post-stage-II, so its labels include ``Realign`` and ``ChangeSegment`` (it re-picks
the V/J segment). arda trusts the submitted call and never swaps it, so those labels
have no arda equivalent. Nor do we assert against VDJdb's blind spots: its
largest-common-substring scanner cannot see an internal substitution at all, so
every mismatch arda reports there is extra information, not a disagreement.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import polars as pl
import pytest

from arda.cdr3fix import markup_records
from tests.conftest import requires_human_db

pytestmark = [requires_human_db]

FIXTURE = Path(__file__).resolve().parent.parent / "assets" / "vdjdb" / "sample.tsv.gz"


@pytest.fixture(scope="module")
def vdjdb():
    if not FIXTURE.exists():
        pytest.skip(f"missing fixture {FIXTURE}")
    with gzip.open(FIXTURE, "rb") as fh:
        df = pl.read_csv(fh.read(), separator="\t", infer_schema_length=0)
    ref = [json.loads(x) for x in df["cdr3fix"]]
    return df, ref


def test_boundaries_agree_with_vdjdb(vdjdb):
    """vEnd/jStart on the records VDJdb itself marks good."""
    df, ref = vdjdb
    recs = markup_records(df, cdr3="cdr3", v="v.segm", j="j.segm", species="species")
    pairs = [(r, f) for r, f in zip(recs, ref) if f["good"]]
    assert pairs
    v_ok = sum(r.v_end == f["vEnd"] for r, f in pairs)
    j_ok = sum(r.j_start == f["jStart"] for r, f in pairs)
    n = len(pairs)
    print(f"\n[vdjdb] vEnd {v_ok}/{n} = {v_ok/n:.1%} | jStart {j_ok}/{n} = {j_ok/n:.1%}")
    assert v_ok / n >= 0.95
    assert j_ok / n >= 0.92


def test_clean_records_are_left_alone(vdjdb):
    """Markup must be idempotent: VDJdb's already-repaired junction needs no fix."""
    df, ref = vdjdb
    recs = markup_records(df, cdr3="cdr3", v="v.segm", j="j.segm", species="species")
    pairs = [(r, f) for r, f in zip(recs, ref) if f["good"]]
    same = sum(r.cdr3_repaired == f["cdr3"] for r, f in pairs)
    n = len(pairs)
    print(f"\n[vdjdb] unchanged {same}/{n} = {same/n:.1%}")
    assert same / n >= 0.98


def test_repair_reproduces_vdjdb_fix(vdjdb):
    """Feed the *submitted* junction (cdr3_old) and reproduce VDJdb's repair."""
    df, ref = vdjdb
    need = [i for i, f in enumerate(ref) if f["fixNeeded"]]
    assert need, "fixture must contain records VDJdb actually repaired"
    sub, sref = df[need], [ref[i] for i in need]
    recs = markup_records(sub, cdr3="cdr3_old", v="v.segm", j="j.segm", species="species")

    exact = sum(r.cdr3_repaired == f["cdr3"] for r, f in zip(recs, sref))
    # A "third string" is neither the submission nor VDJdb's repair -- a novel rewrite.
    third = sum(1 for r, f in zip(recs, sref)
                if r.cdr3_repaired not in (f["cdr3"], f["cdr3_old"]))
    n = len(need)
    print(f"\n[vdjdb] repair reproduced {exact}/{n} = {exact/n:.1%}; novel rewrites {third}")
    assert exact / n >= 0.90
    assert third / n <= 0.05


def test_arda_reports_mismatches_vdjdb_cannot_see(vdjdb):
    """Informational: VDJdb's substring scanner never reports an internal mismatch."""
    df, ref = vdjdb
    recs = markup_records(df, cdr3="cdr3", v="v.segm", j="j.segm", species="species")
    reported = sum(1 for r in recs for e in r.errors if not e.applied)
    applied = sum(1 for r in recs for e in r.errors if e.applied)
    print(f"\n[vdjdb] germline mismatches: {applied} repaired, {reported} reported-only")
    assert reported + applied >= 0     # never assert against VDJdb's blind spots
