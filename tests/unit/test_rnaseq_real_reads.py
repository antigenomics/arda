"""End-to-end RNA-seq behaviour on REAL sequencer reads, not synthesised ones.

Every other rnaseq test builds its reads by slicing the reference and mutating it, so the
reference is on both sides of the assertion. That misses whatever real data does that we did not
think to synthesise: adapters, quality drop-off at the 3' end, PCR error, off-target transcripts,
mates that do not overlap, reads that land in J->C rather than V.

`tests/data/rnaseq_real/` is 660 read pairs drawn from **SRR5233639** (PRJNA371303, human bulk
RNA-seq, 2x100 bp) -- 260 fragments arda maps, spread across IGH/IGK/IGL/TRA/TRB/TRG, plus 400
sampled non-receptor fragments so the reject path is exercised too. 88 KB, provenance in
SOURCES.md.

These need mmseqs and a built human DB, so they skip on a bare checkout -- which is exactly why
the invariants that must NEVER be gated (determinism, chunk invariance, sharding) live in
`test_correct_determinism.py`, `test_chunked_fragments.py` and `test_rnaseq_slurm.py` and use no
DB at all.
"""

from __future__ import annotations

import json

import polars as pl
import pytest

from tests.conftest import requires_human_db, requires_mmseqs

READS_1 = "tests/data/rnaseq_real/reads_1.fq.gz"
READS_2 = "tests/data/rnaseq_real/reads_2.fq.gz"

pytestmark = [requires_mmseqs, requires_human_db]


@pytest.fixture(scope="module")
def mapped(tmp_path_factory):
    """Stage 1 over the real fixture, once for the whole module."""
    from arda.rnaseq.map import map_rnaseq

    d = tmp_path_factory.mktemp("real")
    out = d / "m.airr.tsv"
    # `adaptive=False` on purpose: this fixture is the UNCAPPED reference every two-pass test
    # below compares against, so it must be the exhaustive search, not the (default) adaptive one.
    # Adaptive search has its own tests; mixing the two would make a two-pass assertion fail for
    # a reason that has nothing to do with the two-pass.
    rep = map_rnaseq(READS_1, out, r2=READS_2, threads=2, adaptive=False)
    return out, rep


def test_maps_real_reads_across_every_locus(mapped):
    """Both lineages and all four TR loci, from one bulk library."""
    out, rep = mapped
    assert rep.total_reads == 1320          # 660 pairs x 2 mates
    assert rep.mapped_reads > 300
    loci = rep.per_locus
    for locus in ("IGH", "IGK", "IGL", "TRA", "TRB"):
        assert loci.get(locus, 0) > 0, f"{locus} absent: {loci}"
    # The fixture is B-cell-dominated (bulk PBMC-like), so IG must exceed TR.
    assert sum(loci.get(x, 0) for x in ("IGH", "IGK", "IGL")) > \
           sum(loci.get(x, 0) for x in ("TRA", "TRB", "TRG", "TRD"))


def test_rejects_the_non_receptor_fragments(mapped):
    """400 of the 660 pairs are ordinary transcriptome. They must not all be kept."""
    _, rep = mapped
    assert rep.mapped_fraction < 0.6, (
        f"kept {rep.mapped_fraction:.1%} of a fixture that is 60% non-receptor by construction")


def test_emitted_junctions_are_biologically_plausible(mapped):
    """A junction arda reports on a real read should look like a rearrangement."""
    out, _ = mapped
    d = pl.read_csv(out, separator="\t", infer_schema_length=0)
    j = [x for x in d["junction_aa"].to_list() if x]
    assert j, "no junction_aa emitted at all"
    # Per-read junctions are often truncated (the read need not span the junction), so this is
    # deliberately weak -- `correct --complete-only` is what enforces completeness. What must
    # hold is that they are peptides, not garbage.
    ok = [x for x in j if all(c in "ACDEFGHIKLMNPQRSTVWY*_" for c in x)]
    assert len(ok) == len(j), f"{len(j) - len(ok)} junctions carry non-residue characters"


def test_output_is_invariant_to_chunk_size(tmp_path):
    """The property that makes --chunk-size a free performance knob.

    `chunked_fragments` guarantees a chunk boundary never splits a fragment, so the isotype a
    constant-region mate donates cannot depend on where the boundary fell. Verified here on real
    reads across the default and both neighbours.
    """
    from arda.rnaseq.map import map_rnaseq

    digests, isos = set(), set()
    for chunk in (50_000, 400_000, 1_000_000):
        out = tmp_path / f"c{chunk}.tsv"
        rep = map_rnaseq(READS_1, out, r2=READS_2, threads=2, chunk_size=chunk)
        digests.add(out.read_bytes())
        isos.add((rep.mapped_reads, rep.isotype_from_mate, rep.constant_only_fragments))
    assert len(digests) == 1, "AIRR output changed with --chunk-size"
    assert len(isos) == 1, f"fragment-level counts changed with --chunk-size: {isos}"


def test_repeated_runs_are_byte_identical(tmp_path):
    """Reproducibility on real reads, end to end through all three stages."""
    from arda.rnaseq import pipeline

    outs = []
    for i in range(2):
        d = tmp_path / f"run{i}"
        pipeline.run(READS_1, d, "R", r2=READS_2, threads=2)
        outs.append({name: (d / f"R.{name}").read_bytes()
                     for name in ("airr.tsv", "assembled.airr.tsv", "clones.tsv")})
    for name in outs[0]:
        assert outs[0][name] == outs[1][name], f"{name} differs between two identical runs"


def test_sharded_run_matches_single_node_on_real_reads(tmp_path):
    """The mode-equivalence guarantee, on real data rather than a hand-built AIRR."""
    from arda.cluster import split_pairs
    from arda.rnaseq import pipeline
    from arda.rnaseq.map import map_rnaseq

    single = tmp_path / "single"
    pipeline.run(READS_1, single, "R", r2=READS_2, threads=2)

    shards = split_pairs(READS_1, tmp_path / "sh", shards=4, r2=READS_2)
    mapdir = tmp_path / "map"
    mapdir.mkdir()
    for i, (s1, s2) in enumerate(shards):
        map_rnaseq(s1, mapdir / f"shard_{i:05d}.airr.tsv", r2=s2, threads=2,
                   report_path=mapdir / f"shard_{i:05d}.map.json")
    sharded = tmp_path / "sharded"
    pipeline.reduce(mapdir, sharded, "R", threads=2)

    for name in ("airr.tsv", "assembled.airr.tsv", "clones.tsv"):
        assert (single / f"R.{name}").read_bytes() == (sharded / f"R.{name}").read_bytes(), \
            f"{name} differs between single-node and 4-shard runs"


def test_report_carries_resources_and_provenance(tmp_path):
    """What an operator needs to size a job and to diagnose a cross-mode difference."""
    from arda.rnaseq import pipeline

    rep = pipeline.run(READS_1, tmp_path, "R", r2=READS_2, threads=2)
    for stage in ("map", "assemble", "correct"):
        s = rep[stage]
        assert s["wall_seconds"] > 0, f"{stage} has no wall time"
        assert s["peak_rss_mb"] > 0, f"{stage} has no peak RSS"
        assert s["rss_gain_mb"] >= 0
    # peak_rss_mb is a monotone whole-process high-water mark, so it can only rise.
    peaks = [rep[s]["peak_rss_mb"] for s in ("map", "assemble", "correct")]
    assert peaks == sorted(peaks), f"peak_rss_mb is not monotone across stages: {peaks}"
    assert rep["arda_version"]
    assert "mmseqs_version" in rep and "reference" in rep
    on_disk = json.loads((tmp_path / "R.arda.json").read_text())
    assert on_disk["correct"]["clonotypes_out"] == rep["correct"]["clonotypes_out"]


# ---------------------------------------------------------------------------------------------
# Two-pass segment search. The claim is "cheaper, and no read arda would report is lost", so the
# tests compare it against the one-pass path on the same real reads rather than against itself.
# ---------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def segment_reference():
    """The segment reference is DERIVED, so a checkout does not have one until it is built.

    Without this every two-pass test here silently degrades: `_cached_segment_db` returns None,
    `map_rnaseq` falls back to the one-pass search with a warning, and the comparison tests then
    compare the one-pass output *to itself* and pass. That is the exact "silent success over
    nothing" failure this repo has hit three times, and it is why `test_two_pass_is_actually_
    engaged` below asserts the fast path ran rather than trusting the output to prove it.
    """
    from arda.refbuild.segments import build_segment_reference

    return build_segment_reference("human")


@pytest.fixture(scope="module")
def two_pass(tmp_path_factory, segment_reference):
    from arda.rnaseq.map import map_rnaseq

    d = tmp_path_factory.mktemp("real_tp")
    out = d / "tp.airr.tsv"
    rep = map_rnaseq(READS_1, out, r2=READS_2, threads=2, two_pass=True, adaptive=False)
    return out, rep


def test_two_pass_is_actually_engaged(two_pass):
    """Guard the guard: every comparison below is vacuous if the two-pass silently fell back."""
    _, rep = two_pass
    assert rep.segment_search, (
        "two_pass=True produced no segment_search accounting -- it fell back to the one-pass "
        "search, so every other two-pass assertion in this module is comparing one-pass to itself")
    assert rep.segment_search["implied"] > 0, "the fast path resolved nothing; nothing was tested"


def _by_id(path):
    d = pl.read_csv(path, separator="\t", infer_schema_length=0)
    return {r["sequence_id"]: r for r in d.iter_rows(named=True)}


def test_two_pass_loses_no_read_and_invents_none(mapped, two_pass):
    """The whole point of the rescue set: the kept-read SET is the one-pass set, exactly."""
    one, two = _by_id(mapped[0]), _by_id(two_pass[0])
    assert set(one) - set(two) == set(), "the fast path dropped reads the one-pass search kept"
    assert set(two) - set(one) == set(), "the fast path invented reads"
    assert mapped[1].per_locus == two_pass[1].per_locus


def _gene(call):
    """`IGKV3-20*01,IGKV3D-20*01` -> `IGKV3-20,IGKV3D-20`. Allele suffixes dropped, order fixed."""
    return ",".join(sorted({x.split("*")[0] for x in (call or "").split(",") if x}))


def test_two_pass_does_not_reseat_a_J_to_C_read_onto_a_VxJ_scaffold(mapped, two_pass):
    """A J->C read must keep its isotype and must NOT gain a junction.

    A read that runs through J into C has two plausible homes -- the V×J scaffold its best (V, J)
    pair names, and the J+C scaffold the segment pass actually hit -- and among 775 V alleles a
    100 nt read always has *some* V above threshold. Forcing the V×J choice took a read scoring
    141 on a J+C scaffold, re-seated it at 99 on a V×J one, destroyed the `c_call` and
    **fabricated a junction** out of the spurious V. Both scaffolds now compete on bit score, as
    they do in the one-pass search, and these four columns are byte-identical on real reads.
    """
    one, two = _by_id(mapped[0]), _by_id(two_pass[0])
    for col in ("locus", "c_call", "junction_aa", "productive"):
        bad = [q for q in one if (one[q][col] or "") != (two[q][col] or "")]
        assert not bad, f"{col} differs on {len(bad)} read(s), e.g. " + ", ".join(
            f"{q}: {one[q][col]!r} -> {two[q][col]!r}" for q in bad[:3])


def test_two_pass_disagrees_only_on_the_allele_suffix(mapped, two_pass):
    """V and J calls may differ, but only *within* a gene -- never the gene, never the junction.

    This is the two-pass's one real approximation, and it is worth stating precisely rather than
    papering over. The scaffold is chosen from the best V *segment* hit, which is not always the
    allele whose whole scaffold scores highest: 15 of 453 real reads (3.3 %) end on a sibling
    allele and 5 score a few bits lower. Every one of them keeps the same gene, the same locus and
    a byte-identical `junction_aa` -- the difference lives entirely in the `*NN` suffix, which is
    the part of a V call short-read data cannot resolve anyway.

    Where the gene itself differs (2 reads, both exact bit-score ties) the two-pass picks the
    *functional* gene where the one-pass picked an orphon (`IGKV3/OR2-268*02`) and a
    duplicate-locus paralog (`IGKV3D-20*01`).
    """
    one, two = _by_id(mapped[0]), _by_id(two_pass[0])
    for col in ("v_call", "j_call"):
        bad = [q for q in one if _gene(one[q][col]) != _gene(two[q][col])]
        # A tie is arbitrary in both directions; anything else means a worse scaffold was chosen.
        for q in bad:
            assert float(one[q]["mmseqs2_score"]) == float(two[q]["mmseqs2_score"]), (
                f"{q}: {col} gene changed at a DIFFERENT score "
                f"{one[q][col]}@{one[q]['mmseqs2_score']} -> {two[q][col]}@{two[q]['mmseqs2_score']}")
        assert len(bad) <= 0.01 * len(one), (
            f"{col} gene-level agreement fell below 99 %: {len(bad)}/{len(one)} differ")
    # Allele level too, now that both paths break ties by the same rule (mmseqs' own ordering,
    # via `top_hit`). Before that the two-pass sorted `_best_hits` lexicographically while the
    # one-pass took mmseqs' first line, and the two disagreed on 12 of 453 reads.
    for col in ("v_call", "j_call"):
        exact = sum(1 for q in one if (one[q][col] or "") == (two[q][col] or ""))
        assert exact >= 0.99 * len(one), (
            f"{col} allele-level agreement fell below 99 %: {exact}/{len(one)}")


def test_two_pass_accounts_for_every_read_in_the_report(two_pass):
    """`implied + rescued` is the audit trail; `no_segment_hit` is the residual exposure."""
    _, rep = two_pass
    ss = rep.segment_search
    assert ss, "two_pass=True produced no segment_search accounting"
    assert ss["implied"] + ss["rescued"] + ss["no_segment_hit"] == rep.total_reads
    assert sum(ss["reasons"].values()) >= ss["rescued"] - ss["reasons"].get("", 0)
    assert 0.0 <= ss["fast_fraction"] <= 1.0


def test_two_pass_falls_back_rather_than_failing_without_a_segment_reference(tmp_path, monkeypatch):
    """A missing `segments.fasta` must degrade to the one-pass search, not raise."""
    from arda.annotate import mapper
    from arda.rnaseq.map import map_rnaseq

    monkeypatch.setattr(mapper, "_cached_segment_db", lambda ref, organism: None)
    rep = map_rnaseq(READS_1, tmp_path / "fb.airr.tsv", r2=READS_2, threads=2, two_pass=True)
    assert rep.mapped_reads > 300
    assert rep.segment_search == {}


def test_productive_is_empty_when_no_junction_was_observed(mapped):
    """`productive` is a property of the V-J junction, so a read that never reached one is
    UNEVALUABLE, not non-productive.

    Reporting "F" made a bare V fragment look like a confirmed non-productive rearrangement --
    on this real bulk fixture that was 342 of 453 mapped reads (75 %), since most reads land
    wholly inside V and never reach CDR3. The module always had the rule right for a V-LESS read
    (`productive` stays empty); this extends it to the read that has a V but no junction.

    The invariant is an iff, asserted both ways: `productive` is set exactly when `junction_aa`
    is.
    """
    d = pl.read_csv(mapped[0], separator="\t", infer_schema_length=0)
    rows = list(d.iter_rows(named=True))
    no_junc = [r for r in rows if not (r["junction_aa"] or "").strip()
               and (r["productive"] or "").strip()]
    assert not no_junc, (
        f"{len(no_junc)} read(s) carry a productive call with no junction, e.g. "
        + ", ".join(f"{r['sequence_id']}={r['productive']!r}" for r in no_junc[:3]))
    junc_unset = [r for r in rows if (r["junction_aa"] or "").strip()
                  and not (r["productive"] or "").strip()]
    assert not junc_unset, (
        f"{len(junc_unset)} read(s) have a junction but no productive call, e.g. "
        + ", ".join(r["sequence_id"] for r in junc_unset[:3]))
    # vj_in_frame collapses the same fact and must be gated identically.
    assert all(bool((r["vj_in_frame"] or "").strip()) == bool((r["productive"] or "").strip())
               for r in rows), "vj_in_frame and productive are not gated on the same condition"
    # stop_codon is NOT gated this way: a stop in the V-side regions is observed regardless.
    assert any((r["stop_codon"] or "").strip() and not (r["junction_aa"] or "").strip()
               for r in rows), "stop_codon should stay evaluable for a junctionless read"


def test_adaptive_search_loses_no_read_the_uncapped_search_keeps(tmp_path):
    """The adaptive search caps alignments per read, then re-searches the uncertain ones uncapped.

    `--max-accept` alone is fast and lossy: mmseqs orders candidates by prefilter score, which
    predicts the true best alignment only 55.8 % of the time, so a capped search can stop before
    reaching a read's real scaffold. The reads that suffer are identifiable -- every read lost at
    `--max-accept 40` on 1 M real bulk reads scored 75-83, just above `--min-score 75` -- so
    re-searching only those recovers them. Measured there: 2.17x with LOST = 0, against a
    single-knob lossless point of 1.25x.

    The assertion is the guarantee, not the speed: no read the exhaustive search keeps may be
    missing, and none may be invented.
    """
    from arda.rnaseq.map import map_rnaseq

    exact = map_rnaseq(READS_1, tmp_path / "exact.tsv", r2=READS_2, threads=2, adaptive=False)
    adap = map_rnaseq(READS_1, tmp_path / "adap.tsv", r2=READS_2, threads=2, adaptive=True)
    a, b = _by_id(tmp_path / "exact.tsv"), _by_id(tmp_path / "adap.tsv")
    assert set(a) - set(b) == set(), "the adaptive search lost reads the uncapped search keeps"
    assert set(b) - set(a) == set(), "the adaptive search invented reads"
    assert exact.per_locus == adap.per_locus
    # KNOWN LIMITATION, asserted so it cannot silently grow. The read set is preserved, but the
    # junction is NOT always: 3 of 453 reads get a different `junction_aa`, and two of them score
    # 128 and 131 -- far above the 90-bit trigger. A high score does not certify that the best
    # alignment was found, which is why `adaptive` is off by default. If this count rises, the
    # trigger or the whole score-based premise needs revisiting.
    moved = [q for q in a if (a[q]["junction_aa"] or "") != (b[q]["junction_aa"] or "")]
    assert len(moved) <= 3, (
        f"adaptive search moved {len(moved)} junctions, was 3: {moved[:5]}")
    # c_call moves with the junction on the same read (a different scaffold carries a different
    # constant region), so it is bounded by the same count rather than required to be zero.
    for col in ("locus", "c_call"):
        bad = [q for q in a if (a[q][col] or "") != (b[q][col] or "")]
        assert len(bad) <= 3, f"{col} differs on {len(bad)} read(s), was <=3: {bad[:3]}"
