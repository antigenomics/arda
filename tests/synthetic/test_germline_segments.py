"""Bare germline segment annotation.

A V-only or J-only query carries no V-J junction, but arda applies no coverage
filter, so each maps to its scaffold and only the regions inside the query's
coverage are emitted. This is the path mirpy uses to annotate isolated germline
V and J alleles:

* a bare germline **V** -> fwr1, cdr1, fwr2, cdr2, fwr3 (no fwr4);
* a bare germline **J** -> fwr4 (no V-side regions).
"""

import random

import pytest
import polars as pl

from arda.paths import vdj_dir
from arda.annotate.mapper import annotate_records
from tests.conftest import requires_mmseqs, requires_human_db

pytestmark = [requires_mmseqs, requires_human_db]

_V_REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")


@pytest.fixture(scope="module")
def human_ref():
    from pathlib import Path
    from arda.refbuild.imgt import read_fasta

    fa = dict(read_fasta(Path(vdj_dir("human") / "alleles.fasta")))
    nt = pl.read_csv(vdj_dir("human") / "markup.tsv", separator="\t", infer_schema_length=0)
    ref = {r["scaffold_id"]: r for r in nt.iter_rows(named=True)}
    ids = [i for i in ref if ref[i]["productive"] == "T"
           and int(ref[i]["v_sequence_end"] or 0) > 0
           and int(ref[i]["j_sequence_start"] or 0) > 0]
    random.seed(7)
    random.shuffle(ids)
    return fa, ref, ids[:20]


def test_bare_v_segment_yields_v_side_regions(human_ref):
    fa, ref, ids = human_ref
    queries = []
    meta = {}
    for sid in ids:
        v_end = int(ref[sid]["v_sequence_end"])
        v_seq = fa[sid][:v_end]
        qid = f"{sid}__V"
        queries.append((qid, v_seq))
        meta[qid] = v_seq
    recs = {r["sequence_id"]: r for r in
            annotate_records(queries, organism="human", seqtype="nt",
                             strand="forward", map_d=False)}
    covered = 0
    for qid, v_seq in meta.items():
        rec = recs[qid]
        if not rec.get("fwr1_start"):
            continue  # allele that did not map; allowed for a minority
        covered += 1
        for region in _V_REGIONS:
            assert rec.get(f"{region}_start"), f"{qid} missing {region}"
            s, e = int(rec[f"{region}_start"]), int(rec[f"{region}_end"])
            assert v_seq[s - 1:e] == rec[region]  # round-trip invariant
        # A bare V has no J side.
        assert not rec.get("fwr4_start")
    assert covered >= len(meta) // 2  # the bulk should map


def test_bare_j_segment_yields_fwr4(human_ref):
    fa, ref, ids = human_ref
    queries = []
    meta = {}
    for sid in ids:
        j_start = int(ref[sid]["j_sequence_start"])
        j_seq = fa[sid][j_start - 1:]
        qid = f"{sid}__J"
        queries.append((qid, j_seq))
        meta[qid] = j_seq
    recs = {r["sequence_id"]: r for r in
            annotate_records(queries, organism="human", seqtype="nt",
                             strand="forward", map_d=False)}
    covered = 0
    for qid, j_seq in meta.items():
        rec = recs[qid]
        if not rec.get("fwr4_start"):
            continue
        covered += 1
        s, e = int(rec["fwr4_start"]), int(rec["fwr4_end"])
        assert j_seq[s - 1:e] == rec["fwr4"]  # round-trip invariant
        assert rec.get("fwr4_aa")
    assert covered >= len(meta) // 2
