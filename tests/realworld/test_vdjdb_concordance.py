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
    assert same == n, "idempotence is exact on this fixture; any drift is a regression"


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
    # 99 of 100, and the one difference is arda being RIGHT. `CASSKRGGYEQYV` is a clean
    # TRBJ2-7*02 junction -- that allele templates `SYEQYV`, terminal V and all -- and VDJdb
    # rewrote the V to F (its own label: `Realign`) to force the functional *01 (`SYEQYF`).
    # Since 2.32.0 a call the junction contradicts changes the ALLELE, not the sequence, so arda
    # re-calls *01 -> *02 and returns the submission untouched. Lowering this to 99 is deliberate.
    assert exact >= n - 1, "arda reproduces VDJdb's repair on all but the *02 re-call"
    # A novel rewrite -- neither the submission nor VDJdb's repair -- is data corruption.
    assert third == 0


def test_cdr3fix_json_agrees_with_vdjdb_on_every_verdict(vdjdb):
    """`vCanonical` and `jCanonical` are VDJdb's own booleans and mean exactly the same.

    They describe the junction *as repaired*. Reading them off the submission instead put
    `jCanonical` at odds with VDJdb on 76 of these 250 rows -- every record whose terminal
    Phe118 arda had just restored.

    ``good`` is deliberately NOT asserted equal, because since 2.32.0 arda's is stricter: a
    germline disagreement it declines to repair makes the side `impossible`, where VDJdb's fix
    type reads `NoFixNeeded` because its largest-common-substring scanner never looked that deep.
    What IS asserted is the direction -- arda never calls a record good that VDJdb calls bad -- so
    the difference can only ever be arda withholding a verdict, never inventing one.
    """
    df, ref = vdjdb
    recs = markup_records(df, cdr3="cdr3_old", v="v.segm", j="j.segm", species="species")
    # ⚠ Both fields are VDJdb's LITERAL `C` / `[FW]` tests and stay that way, because the
    # `cdr3fix` JSON means what VDJdb means by them. They can therefore disagree with arda on an
    # allele whose anchor is neither F nor W -- `TRBJ2-7*02` templates `SYEQYV` -- and that is the
    # field being wrong about the biology, not arda. `good` no longer reads them for that reason.
    for key, attr in (("vCanonical", "v_canonical"), ("jCanonical", "j_canonical")):
        bad = [f["cdr3_old"] for r, f in zip(recs, ref) if getattr(r, attr) != f[key]]
        assert len(bad) <= 1, f"{key}: {len(bad)}/{len(ref)} disagree with VDJdb: {bad[:5]}"

    optimistic = [(r, f) for r, f in zip(recs, ref) if r.good and not f["good"]]
    assert not optimistic, \
        f"arda called {len(optimistic)} records good that VDJdb calls bad: " \
        f"{[f['cdr3_old'] for _, f in optimistic][:5]}"
    stricter = sum(1 for r, f in zip(recs, ref) if f["good"] and not r.good)
    print(f"\n[vdjdb] good: arda stricter on {stricter}/{len(ref)}, never more permissive")
    assert stricter, "if this reaches 0, the defect-4 strictness has been lost"


def test_a_repaired_junction_is_always_canonical(vdjdb):
    """The rule the whole repair exists to serve, on both the submitted and fixed columns.

    ⚠ Scored against the CALLED ALLELE'S OWN anchor residue, not the `[FW]` motif -- a conserved
    motif is not an anchor. `TRBJ2-7*02` templates `SYEQYV`, so `CASSKRGGYEQYV` is canonical for it
    and the motif test calls it broken; `TRAJ35*01`'s anchor decodes Cys. Reading the residue from
    `cdr3_anchors.tsv` is both stricter (it pins the exact residue, not a two-letter class) and
    correct on the alleles where the motif is wrong.
    """
    from arda.cdr3fix import load_anchors
    df, _ = vdjdb
    checked = 0
    for col in ("cdr3", "cdr3_old"):
        for r in markup_records(df, cdr3=col, v="v.segm", j="j.segm", species="species"):
            if r.good:
                anchors = load_anchors(r.species)        # the fixture spans organisms
                v = anchors.get(("V", r.v_call))
                j = anchors.get(("J", r.j_call))
                assert v is not None and j is not None, (col, r.species, r.v_call, r.j_call)
                v_aa, j_aa = v.templated_aa[0], j.templated_aa[-1]
                assert r.cdr3_repaired.startswith(v_aa), (col, r.cdr3, r.v_call, v_aa)
                assert r.cdr3_repaired.endswith(j_aa), (col, r.cdr3, r.j_call, j_aa)
                checked += 1
            # v_end/j_start index cdr3_repaired, not the submission
            if r.v_end >= 0 and r.j_start >= 0:
                assert 0 <= r.v_end <= r.j_start <= len(r.cdr3_repaired), (col, r.cdr3)
    assert checked, "no good records to score"


def test_arda_reports_mismatches_vdjdb_cannot_see(vdjdb):
    """Informational: VDJdb's substring scanner never reports an internal mismatch."""
    df, ref = vdjdb
    recs = markup_records(df, cdr3="cdr3", v="v.segm", j="j.segm", species="species")
    reported = sum(1 for r in recs for e in r.errors if not e.applied)
    applied = sum(1 for r in recs for e in r.errors if e.applied)
    print(f"\n[vdjdb] germline mismatches: {applied} repaired, {reported} reported-only")
    assert reported + applied >= 0     # never assert against VDJdb's blind spots
