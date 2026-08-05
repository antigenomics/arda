"""The segment reference must be smaller AND coordinate-correct.

Smaller is easy and worthless on its own — a reference whose region coordinates do not land on
the right bases produces confidently wrong markup. So every test here is about the coordinates
surviving the V/J slice, not about the target count.

Needs a built human reference (it is derived from `markup.tsv` + `alleles.fasta`), so it skips
on a bare checkout. The counting logic that needs no DB is exercised separately below.
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.annotate.reference import REGIONS
from arda.refbuild.segments import SegmentStats, build_segment_reference
from tests.conftest import requires_human_db


def _read_fasta(path):
    seqs, sid, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if sid:
        seqs[sid] = "".join(buf)
    return seqs


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Build into a scratch copy so the shipped reference is never touched."""
    import shutil

    from arda.paths import vdj_dir

    src = vdj_dir("human")
    d = tmp_path_factory.mktemp("segref")
    for name in ("markup.tsv", "alleles.fasta"):
        shutil.copy(src / name, d / name)
    stats = build_segment_reference("human", out_dir=d)
    return d, stats


def test_stats_arithmetic_without_a_reference():
    """SegmentStats is pure arithmetic — no DB needed, so this never skips."""
    s = SegmentStats(v_targets=775, j_targets=124, jc_targets=345, source_scaffolds=15414)
    assert s.total == 1244
    assert round(s.reduction, 2) == 12.39
    assert s.as_dict()["total"] == 1244
    assert SegmentStats().reduction == 0.0          # no divide-by-zero on an empty build


@requires_human_db
def test_collapses_the_vxj_product(built):
    d, stats = built
    assert stats.total < stats.source_scaffolds / 5, (
        f"expected a large reduction, got {stats.source_scaffolds} -> {stats.total}")
    assert stats.v_targets > 100 and stats.j_targets > 10
    assert stats.jc_targets == 345, "the shipped J+C scaffolds must be carried through"
    fasta = _read_fasta(d / "segments.fasta")
    assert len(fasta) == stats.total


@requires_human_db
def test_every_region_coordinate_lands_inside_its_target(built):
    """The whole point: a sliced target must not keep coordinates from before the slice."""
    d, _ = built
    fasta = _read_fasta(d / "segments.fasta")
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    bad = []
    for r in m.iter_rows(named=True):
        L = len(fasta[r["scaffold_id"]])
        for reg in REGIONS:
            s, e = int(r[f"{reg}_start"]), int(r[f"{reg}_end"])
            if s < 0:
                continue                              # -1 = region absent, legitimate
            if not (1 <= s <= e <= L):
                bad.append((r["scaffold_id"], reg, s, e, L))
    assert not bad, f"{len(bad)} coordinates outside their target, e.g. {bad[:3]}"


@requires_human_db
def test_v_targets_carry_the_v_regions_and_no_j(built):
    d, _ = built
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    v = m.filter(pl.col("segment") == "V")
    assert v.height > 0
    for r in v.iter_rows(named=True):
        assert r["v_call"] and not r["j_call"]
        assert int(r["fwr1_start"]) > 0 and int(r["fwr3_end"]) > 0, r["scaffold_id"]
        # FR4 is J-side: a V target must not claim it
        assert int(r["fwr4_start"]) == -1, f"{r['scaffold_id']} claims fwr4"


@requires_human_db
def test_j_targets_carry_fr4_and_it_matches_the_sequence(built):
    """Shifting J coordinates by `j_sequence_start - 1` must land on the real FR4 bases."""
    d, _ = built
    fasta = _read_fasta(d / "segments.fasta")
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    j = m.filter(pl.col("segment") == "J")
    assert j.height > 0
    checked = 0
    for r in j.iter_rows(named=True):
        assert r["j_call"] and not r["v_call"]
        s, e = int(r["fwr4_start"]), int(r["fwr4_end"])
        assert s > 0, f"{r['scaffold_id']} lost fwr4"
        seq = fasta[r["scaffold_id"]][s - 1:e]
        assert len(seq) == e - s + 1
        assert set(seq) <= set("ACGTN"), seq
        checked += 1
    assert checked == j.height


@requires_human_db
def test_ighj1_fr4_is_the_textbook_sequence(built):
    """One hard-coded biological anchor, so a silent off-by-one cannot pass."""
    d, _ = built
    fasta = _read_fasta(d / "segments.fasta")
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    row = next((r for r in m.iter_rows(named=True) if r["j_call"] == "IGHJ1*01"), None)
    if row is None:
        pytest.skip("IGHJ1*01 not in this reference build")
    s, e = int(row["fwr4_start"]), int(row["fwr4_end"])
    assert fasta[row["scaffold_id"]][s - 1:e] == "TGGGGCCAGGGCACCCTGGTCACCGTCTCCTCA"


@requires_human_db
def test_jc_targets_are_carried_through_v_less(built):
    """J+C targets already have no V; `RefEntry.is_jc` and transfer's `t0 > 0` guard rely on it."""
    d, _ = built
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    jc = m.filter(pl.col("segment") == "JC")
    assert jc.height == 345
    for r in jc.iter_rows(named=True):
        assert not r["v_call"], f"{r['scaffold_id']} has a v_call"
        assert r["c_call"], f"{r['scaffold_id']} lost its c_call"


@requires_human_db
def test_targets_are_unique_and_non_empty(built):
    d, stats = built
    fasta = _read_fasta(d / "segments.fasta")
    assert len(set(fasta)) == stats.total, "duplicate target ids"
    assert all(len(s) > 0 for s in fasta.values()), "empty target sequence"
    # A V target should be a substantial chunk of a V gene, not a stub.
    v = [s for i, s in fasta.items() if i.startswith("V|")]
    assert min(len(s) for s in v) > 50, "a V target is implausibly short"
