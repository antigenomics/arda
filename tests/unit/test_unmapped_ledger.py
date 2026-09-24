"""Why a read did not map — the Stage-1 ledger.

`mapped_fraction` alone makes three different diagnoses into one number: "wrong organism",
"this library has almost no receptor content", and "`--min-score` is too strict". The ledger
separates them, and the contract that makes it worth reading is that it **adds up**: every read
arda was handed is in exactly one bucket, or in `mapped_reads`.

A ledger that silently stops balancing is worse than no ledger, so `accounted` is stated by the
report rather than left for the reader to sum — and asserted here on both delivery paths, because
a bucket wired into the single-file path and not into the shard merge would look fine on a laptop
and be wrong on every cluster run.
"""

from __future__ import annotations

from arda.annotate import mapper
from arda.rnaseq.map import RnaseqReport
from arda.rnaseq.pipeline import _merge_map_reports


def _report(**kw) -> dict:
    return _rep(**kw).as_dict()


def _rep(**kw) -> RnaseqReport:
    return RnaseqReport(input="reads.fq", organism="human", **kw)


def test_the_ledger_states_its_own_completeness():
    d = _report(total_reads=100, mapped_reads=40,
                unmapped={"no_hit": 50, "below_min_score": 10})
    assert d["unmapped"]["accounted"] == 100 == d["total_reads"]


def test_a_ledger_that_does_not_balance_says_so_rather_than_looking_fine():
    """The whole point: a bucket someone forgets to wire in shows up as a shortfall."""
    d = _report(total_reads=100, mapped_reads=40, unmapped={"no_hit": 50})
    assert d["unmapped"]["accounted"] == 90
    assert d["unmapped"]["accounted"] != d["total_reads"]


def test_accounted_is_never_summed_twice_when_a_report_is_rebuilt():
    rep = _rep(total_reads=100, mapped_reads=40, unmapped={"no_hit": 60})
    assert rep.as_dict()["unmapped"]["accounted"] == 100
    # `as_dict` is called more than once on a live report — `run()` rewrites the JSON with the
    # whole-run wall time — and an `accounted` that accumulated would drift a little each time.
    assert rep.as_dict()["unmapped"]["accounted"] == 100


def test_a_report_with_no_reads_has_no_ledger_to_state():
    assert _report(total_reads=0, mapped_reads=0)["unmapped"] == {}


# ── the merge ─────────────────────────────────────────────────────────────────────────────────

def test_shards_sum_their_buckets_and_restate_completeness():
    shards = [
        {"total_reads": 100, "mapped_reads": 30, "organism": "human", "paired": True,
         "input": "a.fq", "unmapped": {"no_hit": 60, "below_min_score": 10, "accounted": 100}},
        {"total_reads": 50, "mapped_reads": 20, "organism": "human", "paired": True,
         "input": "b.fq", "unmapped": {"no_hit": 25, "constant_only": 5, "accounted": 50}},
    ]
    merged = _merge_map_reports(shards)
    assert merged["unmapped"] == {"no_hit": 85, "below_min_score": 10, "constant_only": 5,
                                  "accounted": 150}
    assert merged["unmapped"]["accounted"] == merged["total_reads"]


def test_the_merge_recomputes_accounted_rather_than_adding_the_shards_up():
    """Summing `accounted` is right only by coincidence. If a shard's own ledger was short, the
    merged one must be short too — not look balanced because the arithmetic was inherited."""
    shards = [
        {"total_reads": 100, "mapped_reads": 30, "organism": "human", "paired": True,
         "input": "a.fq", "unmapped": {"no_hit": 10, "accounted": 100}},   # shard lies: 30+10=40
    ]
    merged = _merge_map_reports(shards)
    assert merged["unmapped"]["accounted"] == 40
    assert merged["unmapped"]["accounted"] != merged["total_reads"]


def test_a_sample_with_no_ledger_at_all_merges_to_no_ledger():
    shards = [{"total_reads": 10, "mapped_reads": 10, "organism": "human", "paired": False,
               "input": "a.fq"}]
    assert _merge_map_reports(shards)["unmapped"] == {}


# ── the mapper's side ─────────────────────────────────────────────────────────────────────────

def test_the_mapper_takes_the_ledger_and_leaves_it_alone_when_not_asked():
    """`unmapped=None` is the default, so every existing caller is unchanged. The two reasons
    only the mapper can distinguish — `no_hit` and `hit_not_in_reference` — need mmseqs and a
    real reference, so they are asserted end to end in
    `tests/synthetic/test_unmapped_ledger_real.py` rather than against a stub here."""
    import inspect

    sig = inspect.signature(mapper._annotate_chunk)
    assert sig.parameters["unmapped"].default is None
    assert mapper._annotate_chunk([], None, None, "nt", threads=1, sensitivity=5.7,
                                  mm_strand=None) == []
