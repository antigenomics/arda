"""Real-world concordance vs IgBLAST, from committed gzipped GenBank fixtures.

Fixtures under ``tests/data/realworld/<organism>.fasta.gz`` (+ IgBLAST AIRR
reference ``<organism>.igblast.airr.tsv.gz``) are built by
``scripts/build_test_fixtures.py`` — a balanced ~10k mRNA set across all five
organisms and their loci (IG for all; TR for human/mouse). Tests run offline,
needing only mmseqs + the per-organism reference DB.

Asserts, per organism:
  * region concordance with IgBLAST,
  * AIRR CDR3/junction invariants (``junction_aa`` = C...[FW];
    ``cdr3_aa == junction_aa[1:-1]`` exactly; CDR3 matches IgBLAST),
  * trimmed inputs (V-only / J-side fragments) annotate sensibly.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
import polars as pl

from arda.annotate.mapper import annotate_records
from arda.annotate.io import read_sequences
from tests.conftest import requires_mmseqs, requires_human_db

pytestmark = [requires_mmseqs, requires_human_db]

FIXTURES = Path(__file__).resolve().parent.parent / "data" / "realworld"
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3")
ORGANISMS = sorted(p.name[: -len(".fasta.gz")] for p in FIXTURES.glob("*.fasta.gz"))


def _read_airr_gz(path: Path) -> dict[str, dict]:
    with gzip.open(path, "rb") as fh:
        df = pl.read_csv(fh.read(), separator="\t", infer_schema_length=0)
    return {r["sequence_id"]: r for r in df.iter_rows(named=True)}


def _load(org: str):
    fa = FIXTURES / f"{org}.fasta.gz"
    airr = FIXTURES / f"{org}.igblast.airr.tsv.gz"
    if not fa.exists() or not airr.exists():
        pytest.skip(f"missing fixture for {org} (run scripts/build_test_fixtures.py)")
    queries = list(read_sequences(fa))
    igb = _read_airr_gz(airr)
    arda = {r["sequence_id"]: r for r in annotate_records(queries, org, "nt", threads=4)}
    return queries, igb, arda


if not ORGANISMS:
    pytestmark.append(pytest.mark.skip(reason="no fixtures (run scripts/build_test_fixtures.py)"))


@pytest.mark.parametrize("org", ORGANISMS)
def test_region_concordance(org):
    """Compare arda vs IgBLAST on records IgBLAST itself annotates cleanly.

    GenBank contains genomic / partial / non-productive entries that confuse both
    tools; we compare on productive records and skip IgBLAST regions with stops, so
    this measures annotation agreement on real rearrangements, not input junk.
    """
    _, igb, arda = _load(org)
    compared = agree = 0
    for sid, ig in igb.items():
        if ig.get("productive") != "T":
            continue
        ar = arda.get(sid)
        if not ar:
            continue
        for r in REGIONS:
            i = (ig.get(f"{r}_aa") or "").strip()
            a = (ar.get(f"{r}_aa") or "").strip()
            if not i or not a or "*" in i:
                continue
            compared += 1
            agree += (i in a or a in i) if r in ("fwr1", "cdr3") else (i == a)
    assert compared > 0
    frac = agree / compared
    print(f"\n[{org}] region concordance (productive): {agree}/{compared} = {frac:.1%}")
    assert frac >= 0.97


@pytest.mark.parametrize("org", ORGANISMS)
def test_cdr3_junction_airr_invariants(org):
    _, igb, arda = _load(org)
    n = c_start = fw_end = cdr3_ok = 0
    for sid, ar in arda.items():
        ja, c3 = (ar.get("junction_aa") or ""), (ar.get("cdr3_aa") or "")
        # The structural invariant must hold for EVERY emitted junction (by
        # construction): a junction is only emitted with both conserved flanks.
        if ja:
            assert c3 == ja[1:-1], f"{sid}: cdr3 {c3!r} != junction[1:-1] {ja[1:-1]!r}"
        # Biological checks on records IgBLAST calls productive.
        if igb.get(sid, {}).get("productive") != "T" or not ja:
            continue
        n += 1
        c_start += ja.startswith("C")
        fw_end += ja.endswith(("F", "W"))
        cdr3_ok += c3 == (igb.get(sid, {}).get("cdr3_aa") or "")
    assert n > 0
    print(f"\n[{org}] junction C-start {c_start}/{n}, F/W-end {fw_end}/{n}, "
          f"cdr3==igblast {cdr3_ok}/{n}")
    assert c_start / n >= 0.97 and fw_end / n >= 0.97 and cdr3_ok / n >= 0.97


def test_trimmed_inputs():
    """V-only (FR1-FR3) and J-side (CDR3-FR4) fragments annotate without error."""
    queries, _, arda = _load("human")
    qmap = dict(queries)
    sid = next(s for s in arda if arda[s].get("fwr3_end") and arda[s].get("fwr4_end")
               and arda[s].get("fwr1_start"))
    full, ar = qmap[sid], arda[sid]
    v_only = full[int(ar["fwr1_start"]) - 1 : int(ar["fwr3_end"])]
    j_side = full[int(ar["cdr3_start"]) - 1 : int(ar["fwr4_end"])]
    rv = annotate_records([("v", v_only)], "human", "nt", threads=4)[0]
    assert rv["fwr1_aa"] and rv["cdr1_aa"] and rv["fwr3_aa"] and not rv["fwr4_aa"]
    rj = annotate_records([("j", j_side)], "human", "nt", threads=4)[0]
    assert rj["cdr3_aa"] and rj["fwr4_aa"]
