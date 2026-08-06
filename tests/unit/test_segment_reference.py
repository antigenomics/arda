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
    s = SegmentStats(v_targets=775, j_targets=124, c_targets=25, source_scaffolds=15414)
    assert s.total == 924
    assert round(s.reduction, 2) == 16.68
    assert s.as_dict()["total"] == 924
    assert SegmentStats().reduction == 0.0          # no divide-by-zero on an empty build


@requires_human_db
def test_collapses_the_vxj_product(built):
    d, stats = built
    assert stats.total < stats.source_scaffolds / 5, (
        f"expected a large reduction, got {stats.source_scaffolds} -> {stats.total}")
    assert stats.v_targets > 100 and stats.j_targets > 10
    # The 345 J+C scaffolds are a J x C product over 25 distinct C alleles, in which every scaffold
    # of a locus ends in the SAME constant sequence -- 27.7 % of the targets drawing 76.4 % of the
    # alignments. One target per C allele is what is left once that product is collapsed too.
    assert stats.c_targets == 25, "the J+C product must collapse to one target per C allele"
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
def test_c_targets_are_the_constant_region_alone(built):
    """A C target carries `c_call` and nothing else — no V, no J, no regions.

    It replaces the 345 J+C scaffolds, which were a J x C product: every scaffold of a locus ends
    in the SAME constant sequence, so a read reaching C was aligned against all of them to learn
    one `c_call`. The J half is already covered by the `J|` targets, and a read spanning J into C
    hits both.
    """
    d, stats = built
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    c = m.filter(pl.col("segment") == "C")
    assert c.height == stats.c_targets == 25
    assert not m.filter(pl.col("segment") == "JC").height, "J+C scaffolds must not be copied through"
    fasta = _read_fasta(d / "segments.fasta")
    for r in c.iter_rows(named=True):
        assert not r["v_call"], f"{r['scaffold_id']} has a v_call"
        assert not r["j_call"], f"{r['scaffold_id']} has a j_call — the J is its own target"
        assert r["c_call"], f"{r['scaffold_id']} lost its c_call"
        assert r["scaffold_id"] == f"C|{r['c_call']}"
        assert fasta[r["scaffold_id"]], f"{r['scaffold_id']} is empty"
        for reg in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"):
            assert int(r[f"{reg}_start"]) == -1, f"{r['scaffold_id']} claims a {reg}"


@requires_human_db
def test_every_c_allele_of_every_locus_gets_a_target(built):
    """Not only the loci where the C call is informative — coverage is the other job.

    Only IGH's constant genes separate anything worth reporting (11 alleles, 7 classes = the
    isotype); TRA/TRD/IGK have one C allele each, TRB/TRG two, and IGL's seven are all one class.
    Dropping those 14 targets is measurably faster. It also deletes the only segment target a read
    lying wholly inside the constant region can hit, so such a read hits nothing, never enters
    `seen`, and is never rescued: **14 of 453 reads vanish** on the real-read fixture, all of them
    V-less J->C reads.
    """
    d, _ = built
    m = pl.read_csv(d / "segments.markup.tsv", separator="\t", infer_schema_length=0)
    src = pl.read_csv(d / "markup.tsv", separator="\t", infer_schema_length=0)
    want = {r["c_call"] for r in src.iter_rows(named=True)
            if (r["c_call"] or "") and not (r["v_call"] or "")}
    got = {r["c_call"] for r in m.filter(pl.col("segment") == "C").iter_rows(named=True)}
    assert got == want, f"missing C targets: {sorted(want - got)}"
    assert len({a[:3] for a in got}) > 1, "every locus with a constant region must be represented"


@requires_human_db
def test_targets_are_unique_and_non_empty(built):
    d, stats = built
    fasta = _read_fasta(d / "segments.fasta")
    assert len(set(fasta)) == stats.total, "duplicate target ids"
    assert all(len(s) > 0 for s in fasta.values()), "empty target sequence"
    # A V target should be a substantial chunk of a V gene, not a stub.
    v = [s for i, s in fasta.items() if i.startswith("V|")]
    assert min(len(s) for s in v) > 50, "a V target is implausibly short"


def test_segment_targets_resolve_through_the_reference():
    """A `V|`/`J|` target must return a RefEntry, or every segment hit is dropped as unmapped.

    `Reference` is built from `markup.tsv`, which describes scaffolds. Until `segments.markup.tsv`
    was loaded alongside it, `ref.get("V|IGHV3-7*02")` returned None and `_annotate_chunk` dropped
    the read — measured: searching the segment reference produced **0** annotated reads against
    the scaffold reference's 278 on the same input, even though the search found them. The search
    being fast is useless if nothing downstream can resolve what it found.
    """
    from arda.annotate.reference import load_reference

    ref = load_reference("human", "nt")
    scaffold = ref.get("IGH_1053")
    assert scaffold is not None and scaffold.v_call and scaffold.j_call

    v = ref.get("V|IGHV3-7*02")
    assert v is not None, "V segment targets do not resolve — segment hits would be dropped"
    assert v.v_call == "IGHV3-7*02" and v.locus == "IGH"
    assert not v.j_call, "a V segment has no J half"

    j = ref.get("J|IGHJ4*02")
    assert j is not None, "J segment targets do not resolve"
    assert j.j_call == "IGHJ4*02" and not j.v_call


def test_v_segment_region_coordinates_match_their_scaffolds():
    """The claim that makes segment-only annotation possible, asserted rather than assumed.

    A scaffold is `V + pad + J` with the V left-aligned, so a V allele occupies position 1 of both
    its segment and every scaffold built from it. Measured: all 775 V segments agree with their
    scaffolds on `fwr1/cdr1/fwr2/cdr2/fwr3` start and end — 0 differ. A V-only read annotated from
    the segment therefore gets the same region coordinates it would get from a scaffold, and never
    needed the scaffold at all (a scaffold exists to carry a junction such a read does not have).
    """
    import polars as pl

    from arda.annotate.reference import load_reference

    ref = load_reference("human", "nt")
    d = ref.target_fasta.parent
    seg = {r["scaffold_id"].split("|", 1)[1]: r
           for r in pl.read_csv(d / "segments.markup.tsv", separator="\t",
                                infer_schema_length=0).iter_rows(named=True)
           if str(r["scaffold_id"]).startswith("V|")}
    scaf_by_v = {}
    for r in pl.read_csv(d / "markup.tsv", separator="\t",
                         infer_schema_length=0).iter_rows(named=True):
        scaf_by_v.setdefault(r.get("v_call"), r)

    cols = [f"{r}_{e}" for r in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3") for e in ("start", "end")]
    checked = 0
    for allele, srow in seg.items():
        ref_row = scaf_by_v.get(allele)
        if ref_row is None:
            continue
        checked += 1
        for c in cols:
            assert (srow.get(c) or "") == (ref_row.get(c) or ""), (
                f"{allele}: {c} differs between segment ({srow.get(c)}) and scaffold "
                f"({ref_row.get(c)}) — segment-only annotation would shift the regions")
    assert checked >= 700, f"only {checked} V segments compared; the fixture looks wrong"
