"""``arda.cdr3fix`` against generative ground truth.

VDJdb tells us what its own fixer decided; OLGA tells us what actually happened.
Sampling a rearrangement gives the exact number of germline nucleotides the V and
J contributed, so we can check arda's error-vs-boundary heuristic against truth
rather than against another heuristic.

Two properties matter, and they pull in opposite directions:

* **Never under-report.** Residues wholly templated by germline match exactly, so
  the alignment must always credit them. ``v_end >= truth.v_end`` is an invariant,
  not a statistic. (This test caught a real bug: a tie-break on ``i + j`` made the
  aligner prepend ``CA`` to a junction that already began with the conserved Cys.)
* **Never repair a clean sequence.** OLGA junctions contain no typos, so *any*
  applied edit is a false positive. This is the honest price of repairing the
  anchor-adjacent substitution (``CASS`` -> ``CCSS``) that VDJdb's fixer cannot do
  at all, and it is measured here rather than assumed.

Over-reporting by a residue or two is expected and harmless: the codon straddling
the V/N boundary is part-germline, and a chance N-region residue matches germline
about 1 time in 20.

The false-repair rate is small (~0.3%), so a clean run proves little on its own.
``test_injected_typo_is_detected_and_repaired`` is the positive control that keeps
the suite from passing vacuously if markup silently stopped working.

Coverage: OLGA ships human TRA/TRB/IGH/IGK/IGL and mouse TRA/TRB. Setting
``$ARDA_VDJREARM`` adds human TRD and TRG (see SOURCES.md), giving all 7 human loci.
"""

from __future__ import annotations

import os

import pytest

from arda.cdr3fix import load_anchors, markup_cdr3
from tests.conftest import requires_human_db, requires_olga

pytestmark = [requires_olga, requires_human_db]

# Fixed so the reported rates are reproducible; 400 keeps the whole module under a second.
N_SEQS, SEED = 400, 2


@pytest.fixture(scope="module")
def truths():
    """``{(organism, locus): [Truth, ...]}`` for every loadable model."""
    from tests.synthetic.olga_truth import generate, olga_model_dirs

    out = {}
    for model_dir, (org, locus, is_vdj) in sorted(olga_model_dirs().items(),
                                                  key=lambda kv: kv[1]):
        anchors = load_anchors(org)
        if not anchors:
            continue
        ts = generate(model_dir, is_vdj, N_SEQS, anchors, seed=SEED)
        # OLGA's germlines are not IMGT's (its TRBV3-1*01 starts 3 nt later), so a
        # ground-truth comparison is only meaningful where the two agree.
        out[(org, locus)] = [t for t in ts if t.germline_matches_imgt]
    assert out, "no OLGA models loadable"
    return out


def _markups(truths, **kw):
    for (org, _locus), ts in truths.items():
        for t in ts:
            yield t, markup_cdr3(t.cdr3_aa, t.v_call, t.j_call, org, **kw)


def test_models_cover_all_seven_human_loci_when_vdjrearm_present(truths):
    """OLGA alone lacks TRG/TRD; $ARDA_VDJREARM supplies them."""
    human = {locus for (org, locus) in truths if org == "human"}
    assert {"TRA", "TRB", "IGH", "IGK", "IGL"} <= human
    if os.environ.get("ARDA_VDJREARM"):
        assert {"TRD", "TRG"} <= human


def test_every_generated_record_resolves(truths):
    """Guard: if markup silently no-ops, every other assertion here passes vacuously."""
    pairs = list(_markups(truths))
    assert len(pairs) >= 1000
    unresolved = [t.cdr3_aa for t, m in pairs if not (m.v_call and m.j_call and m.v_end >= 0)]
    assert unresolved == [], f"{len(unresolved)} records failed to resolve"


def test_never_under_reports_the_templated_span(truths):
    """Invariant: a residue wholly templated by germline always matches, so the
    alignment must credit it. Under-reporting means the aligner dropped real signal."""
    bad = [(t.cdr3_aa, t.v_call, m.v_end, t.v_end, m.j_start, t.j_start)
           for t, m in _markups(truths)
           if m.v_end < t.v_end or m.j_start > t.j_start]
    assert bad == [], f"{len(bad)} under-reported spans, e.g. {bad[:3]}"


def test_boundaries_land_within_one_residue_of_truth(truths):
    """The straddling codon at the V/N boundary makes +1 common and legitimate."""
    close = total = 0
    for t, m in _markups(truths):
        total += 2
        close += (m.v_end - t.v_end) <= 1
        close += (t.j_start - m.j_start) <= 1
    print(f"\n[olga] boundaries within 1 residue: {close}/{total} = {close/total:.1%}")
    assert close / total >= 0.95


def test_clean_sequences_are_rarely_repaired(truths):
    """OLGA junctions have no typos, so every applied edit is a false positive.

    Measured (seed=2, n=400/model): 11/3180 = 0.35% overall, worst locus ~1%. That is
    the cost of repairing the anchor-adjacent substitution; `max_replace=0` cuts it to
    0.09% and gives up the `CASS` -> `CCSS` fix entirely.
    """
    total_bad = total_n = 0
    for (org, locus), ts in truths.items():
        bad = sum(markup_cdr3(t.cdr3_aa, t.v_call, t.j_call, org).fix_needed for t in ts)
        n = len(ts)
        if not n:
            continue
        total_bad, total_n = total_bad + bad, total_n + n
        print(f"\n[olga] {org} {locus}: false repairs {bad}/{n} = {bad/n:.2%}")
        assert bad / n <= 0.04, f"{org} {locus} repairs clean sequences too often"
    print(f"\n[olga] TOTAL false repairs {total_bad}/{total_n} = {total_bad/total_n:.2%}")
    assert total_bad / total_n <= 0.015


def test_injected_typo_is_detected_and_repaired(truths):
    """Positive control. Mutate the residue beside the Cys -- a real, unambiguous
    error -- and require arda to find it and restore the original."""
    detected = repaired = eligible = 0
    for (org, _locus), ts in truths.items():
        for t in ts:
            m0 = markup_cdr3(t.cdr3_aa, t.v_call, t.j_call, org)
            if m0.v_end < 2 or m0.fix_needed:
                continue                      # need a clean record templating index 1
            eligible += 1
            aa = t.cdr3_aa
            mutant = aa[0] + ("A" if aa[1] != "A" else "G") + aa[2:]
            m = markup_cdr3(mutant, t.v_call, t.j_call, org)
            detected += any(e.side == "V" and e.kind == "sub" and e.pos == 1
                            for e in m.errors)
            repaired += m.cdr3_repaired == aa
    print(f"\n[olga] injected typo @1: detected {detected}/{eligible}, "
          f"repaired to original {repaired}/{eligible}")
    assert eligible >= 200
    assert detected / eligible >= 0.85
    assert repaired / eligible >= 0.85


def test_deep_typo_is_reported_but_not_repaired(truths):
    """Mirror of the above: a typo far from the anchor is localized, never rewritten."""
    reported = changed = eligible = 0
    for (org, _locus), ts in truths.items():
        for t in ts:
            m0 = markup_cdr3(t.cdr3_aa, t.v_call, t.j_call, org)
            aa = t.cdr3_aa
            # a J-templated residue at least 3 from the [FW]118 anchor
            pos = len(aa) - 4
            if m0.fix_needed or m0.j_start < 0 or pos <= m0.j_start or pos <= m0.v_end:
                continue
            eligible += 1
            mutant = aa[:pos] + ("A" if aa[pos] != "A" else "G") + aa[pos + 1:]
            m = markup_cdr3(mutant, t.v_call, t.j_call, org)
            reported += any(e.side == "J" and e.kind == "sub" and not e.applied
                            for e in m.errors)
            changed += m.fix_needed
    print(f"\n[olga] deep typo: reported {reported}/{eligible}, rewritten {changed}")
    assert eligible >= 100
    assert reported / eligible >= 0.80
    assert changed == 0, "a mismatch far from the anchor must never be repaired"


def test_max_replace_zero_is_strictly_more_conservative(truths):
    """The knob does what it says: fewer repairs of clean sequences."""
    strict = sum(m.fix_needed for _, m in _markups(truths, max_replace=0))
    lenient = sum(m.fix_needed for _, m in _markups(truths, max_replace=1))
    print(f"\n[olga] false repairs: max_replace=0 -> {strict}, max_replace=1 -> {lenient}")
    assert strict < lenient


def test_harness_exposes_d_ground_truth():
    """VDJ models must yield the D allele and its exact nt span -- branch 2 needs it."""
    from tests.synthetic.olga_truth import generate, olga_model_dirs

    vdj = [(d, v) for d, v in olga_model_dirs().items() if v[2]]
    assert vdj, "expected at least one VDJ model (human_T_beta)"
    model_dir, (org, _locus, _) = vdj[0]
    ts = generate(model_dir, True, 20, load_anchors(org), seed=3)
    assert all(t.d_call for t in ts)
    for t in ts:
        # A fully trimmed-away D yields an empty span; otherwise it sits inside the junction.
        if t.d_end >= t.d_start:
            assert 0 <= t.d_start <= t.d_end < len(t.cdr3_nt)
