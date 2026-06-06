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


def test_no_hit_yields_empty_record():
    rec = annotate_records([("random", "ACGT" * 40)], "human", "nt", threads=4)[0]
    assert rec["sequence_id"] == "random"
    # A random sequence should not produce region calls.
    assert rec["cdr3"] in ("", None)
