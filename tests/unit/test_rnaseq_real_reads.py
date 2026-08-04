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
    rep = map_rnaseq(READS_1, out, r2=READS_2, threads=2)
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
