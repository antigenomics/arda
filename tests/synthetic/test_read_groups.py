"""A sample split across files must give byte-identical results to the same reads in one file.

That is the whole premise of multi-file support, so it is asserted rather than assumed. arda maps
each read group separately and concatenates the Stage-1 AIRR in declared order; that reproduces
single-file row order exactly, because each file is a contiguous run of the library and
`chunked_fragments` never splits a fragment across a chunk boundary -- so `map` output does not
depend on where the boundaries fall. Stages 2-3 then see byte-identical input and run ONCE.

If the identity test ever fails, multi-file samples are silently producing different clonotypes
than the same data concatenated. That is the failure this feature exists to prevent, not a
tolerance to widen.
"""

import json

import pytest

from arda.cluster import split_pairs
from arda.rnaseq import pipeline
from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

READS_1 = "tests/data/rnaseq_real/reads_1.fq.gz"
READS_2 = "tests/data/rnaseq_real/reads_2.fq.gz"

#: Report keys that legitimately differ between one file and several: timing, memory, and the
#: description of the input itself (four uncompressed chunks are not one gzip). Everything the
#: report *claims about the data* must match.
_VOLATILE = frozenset({
    "wall_seconds", "wall_seconds_max", "wall_seconds_sum", "peak_rss_mb", "peak_rss_mb_max",
    "rss_gain_mb", "reads_per_second", "input", "input_bytes", "shards", "read_groups",
    "reference",
})


def _stable(obj):
    if isinstance(obj, dict):
        return {k: _stable(v) for k, v in obj.items() if k not in _VOLATILE}
    if isinstance(obj, list):
        return [_stable(v) for v in obj]
    return obj


@pytest.fixture(scope="module")
def chunks(tmp_path_factory):
    """The fixture pair cut into four contiguous read groups, standing in for four lanes."""
    d = tmp_path_factory.mktemp("chunks")
    return split_pairs(READS_1, d, shards=4, r2=READS_2)


@pytest.fixture(scope="module")
def whole(tmp_path_factory):
    """`arda rnaseq` over the pair as delivered -- the answer everything else must reproduce."""
    out = tmp_path_factory.mktemp("whole")
    rep = pipeline.run([(READS_1, READS_2)], out, "S", threads=2)
    return out, rep


@pytest.fixture(scope="module")
def regrouped(tmp_path_factory, chunks):
    out = tmp_path_factory.mktemp("regrouped")
    rep = pipeline.run(chunks, out, "S", threads=2)
    return out, rep


def test_four_read_groups_reproduce_the_single_file_run_byte_for_byte(whole, regrouped):
    one, _ = whole
    many, _ = regrouped
    for name in ("airr", "assembled_airr", "clones"):
        f = pipeline.OUTPUTS[name].format(prefix="S")
        assert (many / f).read_bytes() == (one / f).read_bytes(), \
            f"{f} differs between one file and four read groups"


def test_the_report_claims_the_same_thing_about_both_runs(whole, regrouped):
    one, _ = whole
    many, _ = regrouped
    assert _stable(json.loads((many / "S.arda.json").read_text())) == \
        _stable(json.loads((one / "S.arda.json").read_text()))


def test_the_parts_directory_does_not_outlive_a_successful_run(regrouped):
    """It is the Stage-1 scratch. Left behind, the next run's `reduce` could glob it."""
    many, _ = regrouped
    assert not (many / "S.parts").exists()


def test_the_read_group_reports_are_summed_not_lost(regrouped, chunks):
    """`reduce` used to drop the library shape, so a sharded run wrote no `sample` stats scope."""
    _, rep = regrouped
    m = rep["map"]
    assert m["read_groups"] == 4
    assert m["paired"] is True
    assert m["input_bytes"] == sum(p.stat().st_size for pair in chunks for p in pair)
    assert m["read_length_min"] == m["read_length_max"] == m["read_length_mean"] == 100


def test_the_library_shape_reaches_the_stats_table(regrouped):
    """These `run`/`map` rows exist only if Stage 1's own fields survived the merge."""
    many, _ = regrouped
    rows = [r.split("\t") for r in (many / "S.stats.tsv").read_text().splitlines()[1:]]
    m = {r[2]: r[3] for r in rows if r[0] == "run" and r[1] == "map"}
    assert m["paired"] == "1"
    assert m["read_groups"] == "4"
    assert int(m["read_length_min"]) == int(m["read_length_max"]) == 100
    assert int(m["total_reads"]) == 1320


# --------------------------------------------------------- why a sample may not simply be split

def test_two_samples_do_not_share_a_clone_set(tmp_path, chunks):
    """Stages 2-3 are global per SAMPLE, and `correct` is exact: the parts sum to the whole.

    With Stage 3 off, splitting the same reads into two samples partitions them cleanly -- the
    counts add up, which is what says `correct` is not double-counting or dropping anything at a
    sample boundary.
    """
    whole = pipeline.run(chunks, tmp_path / "w", "W", threads=2, assemble=False)
    a = pipeline.run(chunks[:2], tmp_path / "s", "A", threads=2, assemble=False)
    b = pipeline.run(chunks[2:], tmp_path / "s", "B", threads=2, assemble=False)
    assert a["correct"]["reads"] + b["correct"]["reads"] == whole["correct"]["reads"]
    assert a["correct"]["clonotypes_out"] + b["correct"]["clonotypes_out"] == \
        whole["correct"]["clonotypes_out"]


def test_splitting_a_sample_costs_the_contigs_that_tile_across_the_split(tmp_path, chunks):
    """...and with Stage 3 ON they do NOT add up, which is why read groups are merged, not run apart.

    Contigs are grown across reads. Reads that tile one long CDR3 in different *samples* never
    meet, so the contig is never built -- precisely the class Stage 3 exists for. Measured on this
    fixture: 57 reads as one sample against 41 + 17 = 58 as two, because the contigs differ.
    """
    whole = pipeline.run(chunks, tmp_path / "w", "W", threads=2)
    a = pipeline.run(chunks[:2], tmp_path / "s", "A", threads=2)
    b = pipeline.run(chunks[2:], tmp_path / "s", "B", threads=2)
    assert a["correct"]["reads"] + b["correct"]["reads"] != whole["correct"]["reads"]


def test_a_rerun_of_the_same_read_groups_is_byte_identical(tmp_path, chunks):
    """Determinism, over the path that concatenates -- not only over a single input file."""
    outs = []
    for i in range(2):
        d = tmp_path / f"run{i}"
        pipeline.run(chunks, d, "S", threads=2)
        outs.append({n: (d / pipeline.OUTPUTS[n].format(prefix="S")).read_bytes()
                     for n in ("airr", "assembled_airr", "clones")})
    assert outs[0] == outs[1]
