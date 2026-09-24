"""The unmapped ledger against real reads, on both delivery paths.

The contract is that it **adds up**: every read arda was handed is in exactly one bucket or in
`mapped_reads`. That cannot be checked against a stub — the buckets are filled by mmseqs returning
nothing, by the constant-region rule, and by the score cutoff, three places with nothing in common
— so it is checked here, on `tests/data/rnaseq_real`, the fixture every other end-to-end claim in
this repo is measured on.
"""

import json

import pytest

from arda.cluster import split_pairs
from arda.rnaseq import pipeline
from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

READS_1 = "tests/data/rnaseq_real/reads_1.fq.gz"
READS_2 = "tests/data/rnaseq_real/reads_2.fq.gz"


@pytest.fixture(scope="module")
def whole(tmp_path_factory):
    out = tmp_path_factory.mktemp("ledger_whole")
    pipeline.run([(READS_1, READS_2)], out, "L", threads=2, prefilter=True)
    return json.loads((out / "L.arda.json").read_text())["map"]


@pytest.fixture(scope="module")
def lanes(tmp_path_factory):
    chunks = split_pairs(READS_1, tmp_path_factory.mktemp("chunks"), shards=4, r2=READS_2)
    out = tmp_path_factory.mktemp("ledger_lanes")
    pipeline.run(chunks, out, "L", threads=2, prefilter=True)
    return json.loads((out / "L.arda.json").read_text())["map"]


def test_every_read_is_in_exactly_one_bucket(whole):
    led = whole["unmapped"]
    buckets = {k: v for k, v in led.items() if k != "accounted"}
    assert whole["mapped_reads"] + sum(buckets.values()) == whole["total_reads"]
    assert led["accounted"] == whole["total_reads"]


def test_the_buckets_are_the_ones_a_reader_can_act_on(whole):
    """Each names a different remedy: drop the prefilter, check the organism, expect
    constant-only mates on a long insert, lower `--min-score`."""
    assert {"prefilter_rejected", "no_hit", "constant_only", "below_min_score"} <= set(
        whole["unmapped"])


def test_the_constant_rule_is_counted_in_READS_not_fragments(whole):
    """Never: `_apply_constant_rule` also drops the constant-only MATE of a fragment it keeps,
    which `constant_only_fragments` never counts — that is the whole point of the donation. Using
    the fragment count here left 45 of 1,320 reads in no bucket at all."""
    assert whole["unmapped"]["constant_only"] >= whole["constant_only_fragments"]


def test_a_lane_split_sample_reports_the_same_ledger_as_one_file(whole, lanes):
    """A bucket wired into the single-file path and not into the shard merge looks fine on a
    laptop and is wrong on every cluster run."""
    assert lanes["unmapped"] == whole["unmapped"]
    assert lanes["unmapped"]["accounted"] == lanes["total_reads"]


def test_the_ledger_reaches_the_qc_table_as_comparable_fractions(tmp_path, whole):
    """Counts do not compare across samples of different depth; the fractions do."""
    from arda.stats import collect

    (tmp_path / "r.json").write_text(json.dumps({"map": whole}))
    rows = collect(report=tmp_path / "r.json")
    sample = {m: v for s, k, m, v in rows if s == "sample"}
    assert float(sample["unmapped_prefilter_rejected_fraction"]) == pytest.approx(
        whole["unmapped"]["prefilter_rejected"] / whole["total_reads"], rel=1e-5)
    run = {m for s, k, m, v in rows if s == "run"}
    assert "unmapped.no_hit" in run
