"""Real full-length receptor mRNAs from GenBank (human + mouse).

Frozen fixtures in ``tests/data/genbank_receptors{,_mouse}.fa`` (see the README there for provenance
and re-fetch commands). Unlike the synthetic scaffolds these carry real SHM, real allele diversity,
and -- as it turned out -- real vector cloning sites, so they exercise the annotation path the way
user data does. The tests skip cleanly without mmseqs / the reference DB.

The load-bearing test is FR4 agreement between the two scaffold kinds. A `J + C` scaffold now
carries FR4 (derived from IgBLAST's aux file); a V-J scaffold always has. FR4 lies wholly inside J,
so for the *same molecule* the two must report the same FR4 -- otherwise a J->C read and a V-J read
of one clone would disagree on their only shared region.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from arda import paths
from arda.annotate.mapper import annotate_records
from arda.annotate.transfer import AIRR_COLUMNS

from tests.conftest import requires_mmseqs, requires_human_db

_DATA = Path(__file__).parent.parent / "data"
_RC = str.maketrans("ACGTacgtN", "TGCATGCAN")


def _load(fname: str) -> list[tuple[str, str]]:
    recs, name, seq = [], None, []
    for line in (_DATA / fname).read_text().splitlines():
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name:
        recs.append((name, "".join(seq)))
    return recs


def _fr4_agreement(recs, organism):
    """Shared invariant: a V-less J->C fragment cut from each molecule reports an FR4 that contains
    the full-length V-J call. Returns (n_checked, n_truncation_recovered)."""
    full = {r["sequence_id"]: r for r in annotate_records(recs, organism, "nt", threads=8)}
    seqs = dict(recs)
    frags, expect = [], {}
    for acc, r in full.items():
        if not r["fwr4_start"] or not r["fwr4_aa"]:
            continue
        frag = seqs[acc][max(0, int(r["fwr4_start"]) - 1 - 20):][:100]   # spans FR4, runs into C
        if len(frag) >= 80:
            frags.append((acc, frag))
            expect[acc] = r["fwr4_aa"]
    jc = {x["sequence_id"]: x for x in annotate_records(frags, organism, "nt", threads=8)}
    checked = truncation_recovered = 0
    for acc, _ in frags:
        b = jc.get(acc)
        assert b is not None and b["locus"], f"{acc}: J->C fragment did not map"
        assert not b["v_call"], f"{acc}: fragment unexpectedly retained a V"
        assert b["c_call"], f"{acc}: J->C fragment got no isotype (c_call)"
        frag_aa = b["fwr4_aa"] or ""
        assert expect[acc] in frag_aa, f"{acc}: FR4 {frag_aa!r} lost the V-J call {expect[acc]!r}"
        checked += 1
        truncation_recovered += (frag_aa != expect[acc])
    return checked, truncation_recovered


# ============================== human =============================================================

pytestmark = [requires_mmseqs, requires_human_db]


@pytest.fixture(scope="module")
def human():
    recs = _load("genbank_receptors.fa")
    assert len(recs) == 29
    return recs


@pytest.fixture(scope="module")
def human_annot(human):
    return {r["sequence_id"]: r for r in annotate_records(human, "human", "nt", threads=8)}


def test_every_human_record_maps_to_its_locus(human_annot):
    for acc, r in human_annot.items():
        assert r["locus"] in ("IGH", "IGK", "IGL", "TRA", "TRB"), f"{acc}: {r['locus']!r}"
        assert r["j_call"], f"{acc}: no j_call"


def test_full_length_records_are_productive_with_a_canonical_junction(human_annot):
    """A full-length receptor yields a canonical C..[FW] junction and is productive -- except
    PQ879427.1, whose *own GenBank CDS translation* is frameshifted after CDR3 (never reaches
    FGGGTKLTVL). arda calls it non-productive with an out-of-frame ``_`` marker, which is correct."""
    for acc, r in human_annot.items():
        if acc == "PQ879427.1":
            assert r["productive"] == "F" and "_" in (r["junction_aa"] or "")
            continue
        if not r["v_call"]:
            continue                            # BC100294.1: a real V-less 5' truncation
        jaa = r["junction_aa"] or ""
        assert jaa.startswith("C") and jaa[-1] in "FW", f"{acc}: non-canonical junction {jaa!r}"


def test_v_less_truncation_still_gets_isotype(human_annot):
    """BC100294.1 is a real 5'-truncated TRA cDNA: the alignment starts inside J (no V), runs into
    C. No junction, but j_call + c_call -- isotype from a V-less read, which is what C is for."""
    r = human_annot["BC100294.1"]
    assert not r["v_call"] and r["j_call"] and r["c_call"]


def test_human_fr4_agrees_between_vj_and_jc_scaffolds(human):
    """Cut a V-less J->C fragment from each molecule; its FR4 must contain the full-length call.
    ``contains`` not ``equals``: 3 IGH records carry an XhoI (CTCGAG) cloning site a few nt into J,
    so the V-J alignment stops early and truncates FR4 while the J+C scaffold, anchored in C beyond
    the site, recovers it. That recovery is the whole point of extending the reference into C."""
    checked, recovered = _fr4_agreement(human, "human")
    assert checked >= 25
    assert recovered >= 3, "expected the XhoI-bearing IGH records to show FR4 recovery"


def test_airr_alignment_fields_are_populated_and_consistent(human_annot):
    """The AIRR alignment fields added in Phase B: aligned strings are equal length; V-germline
    coords start at 1 and V identity is a fraction on a V-covered read; vj_in_frame agrees with
    productive; a real V-less truncation (BC100294.1) has J-germline coords but no V-side fields."""
    for acc, r in human_annot.items():
        sa, ga = r.get("sequence_alignment") or "", r.get("germline_alignment") or ""
        assert sa and ga and len(sa) == len(ga), f"{acc}: alignment strings absent/unequal"
        if r["v_call"]:
            assert int(r["v_germline_start"]) == 1, f"{acc}: V germline should start at 1"
            vid = float(r["v_identity"])
            assert 0.5 < vid <= 1.0, f"{acc}: implausible v_identity {vid}"
            # productive requires an in-frame V/J; the frameshifted PQ879427.1 must read F/F.
            assert (r["vj_in_frame"] == "T") == (r["productive"] == "T")
    v_less = human_annot["BC100294.1"]
    assert not v_less["v_germline_start"] and not v_less["v_identity"]
    assert int(v_less["j_germline_start"]) >= 1     # J coverage recorded even without V


def test_records_pass_the_airr_rearrangement_schema(human_annot):
    """The payoff of the whole phase: every emitted record validates against the official AIRR
    Rearrangement schema -- all 14 required fields present and correctly typed -- including the
    V-less truncation and the out-of-frame record. arda's TSV is a real AIRR file, not a subset."""
    pytest.importorskip("airr")
    from airr.schema import RearrangementSchema

    missing = [f for f in RearrangementSchema.required if f not in AIRR_COLUMNS]
    assert not missing, f"AIRR_COLUMNS is missing required fields: {missing}"
    for acc, r in human_annot.items():
        RearrangementSchema.validate_row(r)          # raises on any schema violation


def test_reverse_complemented_input_yields_the_same_junction(human, human_annot):
    """arda searches both strands, so an antisense molecule must give the same junction as its
    sense form -- the exact property a stranded library's R2 reads depend on.

    Per AIRR, ``sequence`` keeps the read AS SUBMITTED and ``rev_comp=T`` signals that the output
    data are on its reverse complement. So for an antisense input, ``sequence`` is the (antisense)
    submitted read, while ``sequence_alignment`` and the junction are on the coding strand."""
    rc_in = {acc: seq.translate(_RC)[::-1] for acc, seq in human}
    rc_annot = {r["sequence_id"]: r for r in annotate_records(list(rc_in.items()), "human", "nt", threads=8)}
    compared = 0
    for acc, fwd in human_annot.items():
        if not fwd["junction_aa"]:
            continue
        rev = rc_annot[acc]
        assert rev["junction_aa"] == fwd["junction_aa"], f"{acc}: RC junction differs"
        assert rev["rev_comp"] == "T" and fwd["rev_comp"] == "F"
        assert rev["sequence"] == rc_in[acc], f"{acc}: sequence is not the submitted read"
        # the coding-strand alignment is the reverse complement of the submitted sequence
        assert rev["sequence_alignment"] == fwd["sequence_alignment"]
        compared += 1
    assert compared >= 20


# ============================== mouse =============================================================

_HAS_MOUSE_DB = (paths.vdj_dir("mouse") / "alleles.fasta").exists()
requires_mouse_db = pytest.mark.skipif(not _HAS_MOUSE_DB, reason="mouse reference DB not built")


@pytest.fixture(scope="module")
def mouse():
    recs = _load("genbank_receptors_mouse.fa")
    assert len(recs) == 20
    return recs


@requires_mouse_db
def test_mouse_records_map_and_are_productive(mouse):
    """A whole organism the real-data suite otherwise never touches -- and a different aux file,
    which mixes legacy ``JH1``/``JK1`` names with ``IGHJ1``-style ones."""
    for r in annotate_records(mouse, "mouse", "nt", threads=8):
        assert r["locus"] in ("IGH", "IGK", "TRA", "TRB"), f"{r['sequence_id']}: {r['locus']!r}"
        jaa = r["junction_aa"] or ""
        assert jaa.startswith("C") and jaa[-1] in "FW", f"{r['sequence_id']}: {jaa!r}"


@requires_mouse_db
def test_mouse_fr4_agrees_between_vj_and_jc_scaffolds(mouse):
    """The FR4 derivation is organism-agnostic code over an organism-specific aux file; mouse must
    hold the same invariant as human."""
    checked, _ = _fr4_agreement(mouse, "mouse")
    assert checked >= 18
