"""Real-world concordance vs IgBLAST, from committed GenBank fixtures.

Uses versioned fixtures under ``tests/data/realworld/`` (GenBank mRNA + IgBLAST
AIRR reference, built by ``scripts/build_test_fixtures.py``), so this runs offline
and reproducibly — it needs only mmseqs + the human reference DB, not network or
IgBLAST. Covers IGH (BCR, with D + long CDR3) and TRB (TCR).

Asserts:
  * region concordance with IgBLAST (per locus),
  * the AIRR CDR3/junction invariants — ``junction_aa`` starts with C and ends with
    F/W, and ``cdr3_aa == junction_aa[1:-1]`` exactly (no off-by-one),
  * trimmed inputs (V-only and J-side fragments) annotate sensibly.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import polars as pl

from arda.annotate.mapper import annotate_records
from arda.annotate.io import read_sequences
from tests.conftest import requires_mmseqs, requires_human_db

pytestmark = [requires_mmseqs, requires_human_db]

FIXTURES = Path(__file__).resolve().parent.parent / "data" / "realworld"
LOCI = ["igh", "trb"]
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3")


def _load(locus: str):
    fa = FIXTURES / f"{locus}_mrna.fasta"
    airr = FIXTURES / f"{locus}_mrna.igblast.airr.tsv"
    if not fa.exists() or not airr.exists():
        pytest.skip(f"missing fixture for {locus} (run scripts/build_test_fixtures.py)")
    queries = list(read_sequences(fa))
    igb = {r["sequence_id"]: r for r in pl.read_csv(
        airr, separator="\t", infer_schema_length=0).iter_rows(named=True)}
    arda = {r["sequence_id"]: r for r in annotate_records(queries, "human", "nt", threads=4)}
    return queries, igb, arda


@pytest.mark.parametrize("locus", LOCI)
def test_region_concordance(locus):
    _, igb, arda = _load(locus)
    compared = agree = 0
    for sid, ig in igb.items():
        ar = arda.get(sid)
        if not ar:
            continue
        for r in REGIONS:
            i = (ig.get(f"{r}_aa") or "").strip()
            a = (ar.get(f"{r}_aa") or "").strip()
            if not i or not a:
                continue
            compared += 1
            agree += (i in a or a in i) if r in ("fwr1", "cdr3") else (i == a)
    assert compared > 0
    frac = agree / compared
    print(f"\n[{locus}] region concordance: {agree}/{compared} = {frac:.1%}")
    assert frac >= 0.95


@pytest.mark.parametrize("locus", LOCI)
def test_cdr3_junction_airr_invariants(locus):
    """junction_aa = C...[FW] and cdr3_aa == junction_aa[1:-1], matching IgBLAST.

    The structural invariant ``cdr3_aa == junction_aa[1:-1]`` must hold for EVERY
    record (it is guaranteed by construction). The biological anchors (C-start,
    F/W-end, exact IgBLAST CDR3) must hold for the vast majority; a rare oddball
    record (unusual/non-productive rearrangement) may differ as it does for IgBLAST.
    """
    _, igb, arda = _load(locus)
    n = c_start = fw_end = cdr3_ok = 0
    for sid, ar in arda.items():
        ja, c3 = (ar.get("junction_aa") or ""), (ar.get("cdr3_aa") or "")
        if not ja:
            continue
        n += 1
        # Structural invariant — strict, must always hold.
        assert c3 == ja[1:-1], f"{sid}: cdr3 {c3!r} != junction[1:-1] {ja[1:-1]!r}"
        c_start += ja.startswith("C")
        fw_end += ja.endswith(("F", "W"))
        cdr3_ok += c3 == (igb.get(sid, {}).get("cdr3_aa") or "")
    assert n > 0
    print(f"\n[{locus}] junction C-start {c_start}/{n}, F/W-end {fw_end}/{n}, "
          f"cdr3==igblast {cdr3_ok}/{n}")
    assert c_start / n >= 0.97 and fw_end / n >= 0.97 and cdr3_ok / n >= 0.95


def test_trimmed_inputs():
    """V-only (FR1-FR3) and J-side (CDR3-FR4) fragments annotate without error."""
    queries, _, arda = _load("igh")
    qmap = dict(queries)
    sid = next(s for s in arda if arda[s].get("fwr3_end") and arda[s].get("fwr4_end"))
    full, ar = qmap[sid], arda[sid]
    v_only = full[int(ar["fwr1_start"]) - 1 : int(ar["fwr3_end"])]
    j_side = full[int(ar["cdr3_start"]) - 1 : int(ar["fwr4_end"])]

    rv = annotate_records([("v", v_only)], "human", "nt", threads=4)[0]
    assert rv["fwr1_aa"] and rv["cdr1_aa"] and rv["fwr3_aa"]      # V regions present
    assert not rv["fwr4_aa"]                                       # no J in a V fragment

    rj = annotate_records([("j", j_side)], "human", "nt", threads=4)[0]
    assert rj["cdr3_aa"] and rj["fwr4_aa"]                         # J-side regions present
