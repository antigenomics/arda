"""Synthetic rearrangement tests.

Take real human reference scaffolds, apply controlled mutations (substitutions,
in-frame indels, N runs), annotate them through the full mmseqs+C++ pipeline,
and assert:

* the **round-trip invariant** — ``query[start-1:end] == region_seq`` — which must
  hold for every covered region regardless of how mmseqs places the alignment;
* analytic coordinate behaviour — substitutions don't move coordinates; an
  insertion of L nt inside CDR3 shifts FR4 downstream by L and preserves the
  ``[FW]GXG`` motif.
"""

import random
import re

import pytest
import polars as pl

from arda.paths import vdj_dir
from arda.annotate.mapper import annotate_records
from tests.conftest import requires_mmseqs, requires_human_db

REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
pytestmark = [requires_mmseqs, requires_human_db]


@pytest.fixture(scope="module")
def human_ref():
    from pathlib import Path
    from arda.refbuild.imgt import read_fasta

    fa = dict(read_fasta(Path(vdj_dir("human") / "alleles.fasta")))
    nt = pl.read_csv(vdj_dir("human") / "markup.tsv", separator="\t", infer_schema_length=0)
    ref = {r["scaffold_id"]: r for r in nt.iter_rows(named=True)}
    # Prefer productive scaffolds with no N in regions for clean checks.
    ids = [i for i in ref if ref[i]["productive"] == "T"]
    random.seed(42)
    random.shuffle(ids)
    return fa, ref, ids


def _round_trip_ok(query, rec):
    for r in REGIONS:
        if rec.get(r):
            s, e = int(rec[f"{r}_start"]), int(rec[f"{r}_end"])
            if query[s - 1 : e] != rec[r]:
                return False
    return True


def test_round_trip_invariant_with_substitutions(human_ref):
    fa, ref, ids = human_ref
    queries = []
    chosen = ids[:25]
    rng = random.Random(1)
    for sid in chosen:
        s = list(fa[sid])
        # 3 substitutions inside FR3 (won't move boundaries)
        f3s, f3e = int(ref[sid]["fwr3_start"]), int(ref[sid]["fwr3_end"])
        for _ in range(3):
            p = rng.randint(f3s, f3e) - 1
            s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
        queries.append((sid, "".join(s)))
    out = annotate_records(queries, "human", "nt", threads=4)
    qmap = dict(queries)
    for rec in out:
        assert _round_trip_ok(qmap[rec["sequence_id"]], rec)


def test_insertion_in_cdr3_shifts_fr4(human_ref):
    fa, ref, ids = human_ref
    # Pick scaffolds whose CDR3 region is annotated and FR4 present.
    sid = next(i for i in ids if ref[i]["cdr3_start"] and ref[i]["fwr4_start"])
    s = fa[sid]
    c3s = int(ref[sid]["cdr3_start"])
    fr4s_ref = int(ref[sid]["fwr4_start"])
    mut = s[: c3s + 2] + "AAA" + s[c3s + 2 :]  # +3 nt inside CDR3
    rec = annotate_records([("mut", mut)], "human", "nt", threads=4)[0]
    assert _round_trip_ok(mut, rec)
    assert rec["cdr3"] and rec["fwr4"]
    # FR4 should shift downstream by exactly 3 nt.
    assert int(rec["fwr4_start"]) == fr4s_ref + 3
    # FR4 protein still starts with the conserved [FW]GXG motif.
    assert re.match(r"[FW]G.G", rec["fwr4_aa"])


def test_mmseqs2_score_exposed(human_ref):
    # The mmseqs2 alignment quality of the chosen scaffold hit is surfaced on the record:
    # a strong, exact hit scores well above zero, near-perfect identity, tiny E-value.
    fa, ref, ids = human_ref
    sid = ids[0]
    rec = annotate_records([(sid, fa[sid])], "human", "nt", threads=4)[0]
    assert isinstance(rec["mmseqs2_score"], float) and rec["mmseqs2_score"] > 0
    assert rec["mmseqs2_identity"] > 90.0
    assert rec["mmseqs2_evalue"] < 1e-3
    # A query with no germline similarity yields no hit and a blank (filterable) score.
    miss = annotate_records([("junk", "ACGT" * 30)], "human", "nt", threads=4)[0]
    assert miss["mmseqs2_score"] == "" and miss["locus"] == ""


def test_reverse_complement_strand(human_ref):
    from arda.refbuild.translate import reverse_complement

    fa, ref, ids = human_ref
    sid = ids[0]
    fwd = annotate_records([(sid, fa[sid])], "human", "nt", threads=4)[0]
    rc = annotate_records([("rc", reverse_complement(fa[sid]))], "human", "nt", threads=4)[0]
    assert fwd["rev_comp"] == "F" and rc["rev_comp"] == "T"
    # Region calls must be identical regardless of input strand.
    for k in ("v_call", "j_call", "cdr1_aa", "cdr3_aa", "fwr4_aa"):
        assert fwd[k] == rc[k]


def _aa_ref():
    from pathlib import Path
    from arda.refbuild.imgt import read_fasta
    fa = dict(read_fasta(Path(vdj_dir("human") / "alleles.aa.fasta")))
    aa = pl.read_csv(vdj_dir("human") / "markup.aa.tsv", separator="\t", infer_schema_length=0)
    ref = {r["scaffold_id"]: r for r in aa.iter_rows(named=True)}
    return fa, ref


def test_region_deletion_nt(human_ref):
    """Deleting an internal region (CDR2) collapses it; neighbours stay correct."""
    fa, ref, ids = human_ref
    sid = next(i for i in ids if i.startswith("IGH_") and ref[i]["cdr2_start"])
    s, r = fa[sid], ref[sid]
    c2s, c2e = int(r["cdr2_start"]), int(r["cdr2_end"])
    deleted = s[: c2s - 1] + s[c2e:]                       # remove whole CDR2
    out = annotate_records([("del", deleted)], "human", "nt", threads=4)[0]
    # Flanking and distal regions remain present and round-trip exactly.
    for reg in ("fwr1", "cdr1", "fwr2", "fwr3", "cdr3", "fwr4"):
        if out[reg]:
            s_, e_ = int(out[f"{reg}_start"]), int(out[f"{reg}_end"])
            assert deleted[s_ - 1 : e_] == out[reg]
    # The deleted region is gone (empty or collapsed to <= a boundary residue).
    assert len(out["cdr2_aa"]) <= 1
    # CDR3/FR4 still annotated downstream of the deletion, with valid junction.
    assert out["cdr3_aa"] and out["fwr4_aa"]
    assert out["junction_aa"].startswith("C") and out["junction_aa"].endswith(("F", "W"))


def test_region_deletion_aa():
    """Same, for protein input against the aa reference."""
    fa, ref = _aa_ref()
    sid = next(i for i in fa if i.startswith("IGH_") and ref[i]["cdr2_start"]
               and ref[i].get("fwr4"))
    p, r = fa[sid], ref[sid]
    c2s, c2e = int(r["cdr2_start"]), int(r["cdr2_end"])
    deleted = p[: c2s - 1] + p[c2e:]
    out = annotate_records([("del", deleted)], "human", "aa", threads=4)[0]
    for reg in ("fwr1", "cdr1", "fwr2", "fwr3", "cdr3", "fwr4"):
        if out[reg]:
            s_, e_ = int(out[f"{reg}_start"]), int(out[f"{reg}_end"])
            assert deleted[s_ - 1 : e_] == out[reg]
    assert len(out["cdr2"]) <= 1


def test_out_of_frame_junction(human_ref):
    """A 1-nt CDR3 insertion makes V/J out of frame: junction is still reported,
    translated with an N-bridge ('_' codon), FR4 stays readable, and the V/J split
    positions are populated."""
    fa, ref, ids = human_ref
    sid = next(i for i in ids if i.startswith("IGH_") and ref[i]["cdr3_start"])
    s = fa[sid]
    c3s = int(ref[sid]["cdr3_start"])
    frameshift = s[: c3s + 4] + "A" + s[c3s + 4:]   # +1 nt inside CDR3
    out = annotate_records([("oof", frameshift)], "human", "nt", threads=4)[0]
    assert out["productive"] == "F"                  # out of frame
    assert out["junction_aa"]                         # still reported
    assert "_" in out["junction_aa"]                  # N-bridge marker
    assert out["cdr3_aa"] == out["junction_aa"][1:-1]  # invariant holds with '_'
    assert re.match(r"[FW]G.G", out["fwr4_aa"])        # FR4 still reads correctly
    assert int(out["v_sequence_end"]) > 0 and int(out["j_sequence_start"]) > 0


def test_vj_positions_reported(human_ref):
    """v_sequence_end and j_sequence_start are transferred and reported for the
    vast majority of hits (a few scaffolds lack the extended igblast markup)."""
    fa, ref, ids = human_ref
    out = annotate_records([(i, fa[i]) for i in ids[:40]], "human", "nt", threads=4)
    have = [r for r in out if r["v_call"]]
    assert have
    both = sum(1 for r in have if str(r["v_sequence_end"]).isdigit()
               and str(r["j_sequence_start"]).isdigit())
    assert both / len(have) >= 0.9


@pytest.fixture(scope="module")
def human_d_germlines():
    from arda.annotate.reference import load_reference
    return load_reference("human", "nt").d_germlines


def _gene(call):
    return call.split(",")[0].split("*")[0] if call else ""


def _igh_scaffold_with_pad(fa, ref, ids):
    """An IGH scaffold whose V/J germline-end markup is present (so we can build a
    synthetic junction in its N-pad)."""
    sid = next(i for i in ids if i.startswith("IGH_")
               and str(ref[i].get("v_sequence_end") or "").isdigit()
               and str(ref[i].get("j_sequence_start") or "").isdigit())
    r = ref[sid]
    return sid, fa[sid], int(r["v_sequence_end"]), int(r["j_sequence_start"])


def test_d_segment_recovered(human_ref, human_d_germlines):
    """Build a synthetic junction (np + full D + np) in a scaffold's V..J interior
    and recover d_call + d_sequence_start/end (AIRR, query coords)."""
    fa, ref, ids = human_ref
    # Longest IGH D -> an unambiguous, full-length match.
    d_allele, d_seq = max(human_d_germlines["IGH"], key=lambda x: len(x[1]))
    _, s, v_end, j_start = _igh_scaffold_with_pad(fa, ref, ids)
    # Replace the N-pad between the V and J germline ends with np + D + np.
    query = s[:v_end] + "CAGAT" + d_seq + "ACTGG" + s[j_start - 1 :]
    rec = annotate_records([("dins", query)], "human", "nt", threads=4)[0]
    assert rec["d_call"], "no D segment called"
    assert _gene(rec["d_call"]) == _gene(d_allele)
    ds, de = int(rec["d_sequence_start"]), int(rec["d_sequence_end"])
    # The recovered span is essentially the inserted D (gapless local alignment may
    # extend by a base if a flanking np base coincidentally matches a D terminus).
    assert abs((de - ds + 1) - len(d_seq)) <= 3


def test_double_d_junction(human_d_germlines):
    """D-D fusion mapping: a V..J interior holding two tandem D segments yields
    d_call + d2_call ordered 5'->3' with np1/np2/np3 between V, the Ds, and J.

    This exercises the mapping logic directly: an ~85 nt synthetic junction
    (two full germline D + bridges) is longer than mmseqs reliably aligns *through*
    when projecting the J germline start onto the query, so we feed the interior
    coordinates to ``_map_d`` rather than round-tripping the whole read."""
    from arda.annotate.transfer import _map_d

    igh = sorted(human_d_germlines["IGH"], key=lambda x: len(x[1]), reverse=True)
    (a1, d1) = igh[0]
    # Second D from a *different* gene, to keep the two calls distinguishable.
    (a2, d2) = next((al, sq) for al, sq in igh if _gene(al) != _gene(a1))
    np1, np2, np3 = "CAGAT", "TTAACGGTTAAC", "ACTGG"
    vpref, jsuf = "ACGT" * 6, "TGCA" * 6                  # flanking V / J germline
    junction = np1 + d1 + np2 + d2 + np3
    query = vpref + junction + jsuf
    v_end = len(vpref)                                    # last V base (1-based)
    j_start = len(vpref) + len(junction) + 1             # first J base (1-based)

    rec = {}
    _map_d(rec, query, "IGH", v_end, j_start, human_d_germlines["IGH"])
    assert rec.get("d_call") and rec.get("d2_call"), "expected two D segments"
    # d_call is the 5' (first) D, d2_call the 3' one.
    assert _gene(rec["d_call"]) == _gene(a1)
    assert _gene(rec["d2_call"]) == _gene(a2)
    # Spans are ordered, non-overlapping, and recover the inserted germline D.
    assert int(rec["d_sequence_end"]) < int(rec["d2_sequence_start"])
    assert query[int(rec["d_sequence_start"]) - 1 : int(rec["d_sequence_end"])] == d1
    assert query[int(rec["d2_sequence_start"]) - 1 : int(rec["d2_sequence_end"])] == d2
    # np regions partition the junction between V, the two D, and J.
    assert rec["np1"] == np1 and rec["np2"] == np2 and rec["np3"] == np3


def test_d_mapping_option_toggles_output(human_ref, human_d_germlines):
    """`map_d=False` suppresses all D output; `map_d=True` (default) restores it,
    leaving the V/J markup identical either way."""
    fa, ref, ids = human_ref
    d_allele, d_seq = max(human_d_germlines["IGH"], key=lambda x: len(x[1]))
    _, s, v_end, j_start = _igh_scaffold_with_pad(fa, ref, ids)
    query = s[:v_end] + "CAGAT" + d_seq + "ACTGG" + s[j_start - 1 :]

    on = annotate_records([("q", query)], "human", "nt", threads=4, map_d=True)[0]
    off = annotate_records([("q", query)], "human", "nt", threads=4, map_d=False)[0]

    assert on["d_call"] and on["d_sequence_start"]
    # Disabled: every D/np field is empty.
    for c in ("d_call", "d2_call", "d_sequence_start", "d_sequence_end",
              "d2_sequence_start", "d2_sequence_end", "np1", "np2", "np3"):
        assert (off.get(c) or "") == ""
    # V/J annotation is unaffected by the D option.
    for c in ("v_call", "j_call", "cdr3_aa", "fwr4_aa", "v_sequence_end",
              "j_sequence_start"):
        assert on[c] == off[c]


def test_no_hit_yields_empty_record():
    rec = annotate_records([("random", "ACGT" * 40)], "human", "nt", threads=4)[0]
    assert rec["sequence_id"] == "random"
    # A random sequence should not produce region calls.
    assert rec["cdr3"] in ("", None)
