"""The junction boundary in GERMLINE coordinates must be in the AIRR, not only in the reference.

⛔ Why this exists. `v_mutations` positions are 1-based in the called V allele and `j_mutations`
positions are 1-based in the called J allele, and **both lists span the junction**: the V germline's
3' tail and the J germline's 5' head lie inside it, so exonuclease chew-back and non-templated N/P
bases are emitted as substitutions against a germline that does not template them.

Separating framework from junction needs exactly two numbers -- the Cys104 offset in the V allele
and the [FW]118 offset in the J allele -- and until these columns existed they lived only in
`cdr3_anchors.tsv` inside the reference, so no downstream consumer could do the split from the TSV
alone.

Measured on SRR5233636 (a TRA amplicon, 500,000 reads; T-cell receptors do not hypermutate, so
every entry there is spurious by construction): 1.046 V and 1.658 J entries per read, with 86.2 % of
J entries at J position <= 10. Against the anchors, the recurrent ones are junction-internal:
`TRAV8-6*01` positions 281/282 at frequency 0.88 with its anchor at 270, and `TRAJ8*01` position 1
at 0.67 with its anchor at 26 -- an allele difference in the templated V/J tail, which is neither
somatic mutation nor N/P diversity. ⚠ Frequency alone does NOT separate those classes; frequency
AND position against the anchor does. That is why the columns are emitted raw and unclassified.
"""

from __future__ import annotations

import pytest

from arda.annotate.transfer import AIRR_COLUMNS


def test_the_anchor_columns_are_appended_last():
    """⛔ New columns go LAST. Adding one mid-list silently shifts every later column for a
    consumer that reads the shipped set by position -- which has happened here before."""
    assert AIRR_COLUMNS[-2:] == ("v_anchor_nt", "j_anchor_nt") or \
        list(AIRR_COLUMNS[-2:]) == ["v_anchor_nt", "j_anchor_nt"]


def test_the_anchors_classify_the_measured_recurrent_variants():
    """The two columns must reproduce the split measured on SRR5233636.

    These are the real recurrent variants from that library, with the real anchor offsets: they are
    junction-internal, which is why frequency alone read them as 'alleles' and position is needed.

    ⛔ Deliberately NOT marked `requires_human_db`. In this project that skip is how a reference
    defect regressed unnoticed; the guard below skips only when the reference genuinely is not
    built, which is visible in the report.
    """
    from arda.cdr3fix import load_anchors

    anchors = load_anchors("human")
    if not anchors:
        pytest.skip("human reference not built")

    v = anchors.get(("V", "TRAV8-6*01"))
    j = anchors.get(("J", "TRAJ8*01"))
    if v is None or j is None:
        pytest.skip("alleles absent from this reference build")

    assert v.anchor_nt == 270, "TRAV8-6*01 Cys104 offset moved; the measurement below is keyed to it"
    assert j.anchor_nt == 26, "TRAJ8*01 [FW]118 offset moved"

    # v_mutations position p (1-based in the V allele) is junction-internal iff p > v_anchor_nt.
    for p in (281, 282):                      # measured at frequency 0.88
        assert p > v.anchor_nt, f"V position {p} should be inside the junction"
    for p in (226, 254):                      # measured at low frequency, framework
        assert p <= v.anchor_nt, f"V position {p} should be framework"

    # j_mutations position p is junction-internal iff p <= j_anchor_nt + 3 (the anchor codon is 3 nt
    # and the junction INCLUDES it -- junction != CDR3).
    for p in (1, 5):                          # measured at frequency 0.67 / 0.66
        assert p <= j.anchor_nt + 3, f"J position {p} should be inside the junction"


def test_framework_only_identity_is_what_arda_now_emits():
    """Since 2.16.0 arda applies the framework scope itself — this checks it on its own output.

    ``v_identity`` used to run to ``t_vend``, the V germline's END, which is PAST Cys104, so it was
    depressed by junction diversity rather than by mutation load. It now stops at ``v_anchor_nt``,
    and the anchors still ship so the split stays checkable rather than merely trusted: recomputing
    it from ``v_anchor_nt`` + ``v_germline_start/end`` + ``v_mutations`` must reproduce the emitted
    number.

    ⛔ The two TR records are the check that matters: T-cell receptors do not hypermutate, so a
    framework-only identity below 1.0 there would mean the scope is wrong.
    """
    import re
    from pathlib import Path

    import polars as pl

    tsv = Path(__file__).resolve().parents[2] / "examples" / "example.airr.tsv"
    if not tsv.exists():                      # pragma: no cover - source checkout without examples
        pytest.skip("examples/example.airr.tsv not present")

    mut_re = re.compile(r"^([ACGTN])(\d+)([ACGTN])$")
    df = pl.read_csv(tsv, separator="\t", infer_schema_length=0)
    seen = {}
    for r in df.to_dicts():
        if not r.get("v_call") or not r.get("v_anchor_nt"):
            continue
        gs, ge = int(r["v_germline_start"]), int(r["v_germline_end"])
        anchor = int(r["v_anchor_nt"])
        pos = [int(m.group(2)) for m in
               (mut_re.match(t) for t in (r.get("v_mutations") or "").split(",")) if m]
        span = max(0, min(ge, anchor) - gs + 1)
        assert span > 0, f"{r['v_call']}: no framework left to measure"
        fw = sum(1 for p in pos if p <= anchor)
        seen[r["v_call"].split(",")[0]] = round(1 - fw / span, 4)

    # TR cannot hypermutate: framework-only identity must be exactly 1.
    assert seen.get("TRBV28*02") == 1.0, seen
    assert seen.get("TRAV12-2*02") == 1.0, seen
    # ...and that is now what arda REPORTS, not merely what a downstream recipe can recover.
    # Before 2.16.0 this record's v_identity was 0.8723, depressed purely by junction diversity.
    trb = next(r for r in df.to_dicts() if (r.get("v_call") or "").startswith("TRBV28*02"))
    assert float(trb["v_identity"]) == 1.0, "v_identity is still scoped past Cys104"
    # No v_mutations entry may sit inside the junction on ANY record of the committed example.
    for r in df.to_dicts():
        if not r.get("v_anchor_nt"):
            continue
        anchor = int(r["v_anchor_nt"])
        bad = [m.group(0) for m in
               (mut_re.match(t) for t in (r.get("v_mutations") or "").split(",")) if m
               and int(m.group(2)) > anchor]
        assert not bad, f"{r['sequence_id']}: junction-internal v_mutations {bad}"
