"""No read may be lost by the fast path.

The segment shortlist exists to make alignment cheap, but arda's claim is near-zero Stage-1
false negatives, so a fast path that silently drops reads is a different tool. Every test here
is about the partition being *total* — each read lands in exactly one of `implied` / `rescue` —
rather than about the speedup.

No DB, no mmseqs: the partition is pure logic over two dicts and a lookup table.
"""

from __future__ import annotations

import pytest

from arda.annotate.shortlist import Shortlist, shortlist

COMBOS = {
    ("TRAV1-2*01", "TRAJ33*01"): "TRA_1",
    ("TRBV20-1*01", "TRBJ2-1*01"): "TRB_7",
    ("IGHV3-21*06", "IGHJ4*02"): "IGH_9",
    ("TRAV23/DV6*01", "TRAJ12*01"): "TRA_5",
    ("TRAV23/DV6*01", "TRDJ1*01"): "TRD_2",
}


def test_a_clean_pair_takes_the_fast_path():
    sl = shortlist({"r1": "TRBV20-1*01"}, {"r1": "TRBJ2-1*01"}, COMBOS)
    assert sl.implied == {"r1": "TRB_7"}
    assert sl.rescue == []
    assert sl.fast_fraction == 1.0


@pytest.mark.parametrize("bv,bj,failed,reason", [
    ({"r": "TRAV1-2*01"}, {}, None, "v_only"),                     # never reached a J
    ({}, {"r": "IGHJ4*02"}, None, "j_only"),                       # J->C read, no V
    ({"r": "TRAV1-2*01"}, {"r": "TRBJ2-1*01"}, None, "no_such_combination"),
    ({"r": "TRBV20-1*01"}, {"r": "TRBJ2-1*01"}, {"r"}, "second_pass_failed"),
])
def test_every_unresolvable_read_is_rescued_with_a_reason(bv, bj, failed, reason):
    sl = shortlist(bv, bj, COMBOS, failed=failed)
    assert sl.implied == {}, "an unresolvable read must not take the fast path"
    assert sl.rescue == ["r"]
    assert sl.reasons == {reason: 1}


def test_trav_dv_crossing_the_locus_is_rescued_not_guessed():
    """The measured failure: 54% of residual V disagreements were TRA->TRD.

    TRAV/DV segments are genuinely dual-use, and a joint argmax over V and J did not resolve
    them. Until a locus prior exists, these are rescued rather than guessed.
    """
    sl = shortlist({"r": "TRAV23/DV6*01"}, {"r": "TRDJ1*01"}, COMBOS)
    assert sl.rescue == ["r"]
    assert sl.reasons == {"trav_dv_locus_ambiguous": 1}
    # ...but a TRAV/DV read that stays inside TRA is fine and must NOT be rescued.
    ok = shortlist({"r": "TRAV23/DV6*01"}, {"r": "TRAJ12*01"}, COMBOS)
    assert ok.implied == {"r": "TRA_5"} and ok.rescue == []


def test_the_partition_is_total_over_a_mixed_population():
    """The invariant: implied + rescue accounts for every read, with no overlap."""
    best_v = {f"v{i}": "TRAV1-2*01" for i in range(30)}
    best_j = {f"j{i}": "IGHJ4*02" for i in range(20)}
    both = {f"b{i}": "IGHV3-21*06" for i in range(50)}
    best_v.update(both)
    best_j.update({f"b{i}": "IGHJ4*02" for i in range(50)})
    failed = {"b0", "b1", "b2"}

    sl = shortlist(best_v, best_j, COMBOS, failed=failed)
    seen = set(sl.implied) | set(sl.rescue)
    assert seen == set(best_v) | set(best_j), "some read is in neither partition"
    assert not (set(sl.implied) & set(sl.rescue)), "a read is in BOTH partitions"
    assert sl.n_total == 100
    assert len(sl.implied) == 47                      # 50 paired minus 3 failed
    assert sum(sl.reasons.values()) == len(sl.rescue)


def test_no_read_is_lost_even_when_nothing_resolves():
    sl = shortlist({f"r{i}": "TRAV1-2*01" for i in range(25)}, {}, COMBOS)
    assert len(sl.rescue) == 25 and sl.implied == {}
    assert sl.fast_fraction == 0.0                    # and no ZeroDivisionError


def test_empty_input_is_not_an_error():
    sl = shortlist({}, {}, COMBOS)
    assert sl.n_total == 0 and sl.fast_fraction == 0.0
    assert sl.as_dict() == {"implied": 0, "rescued": 0, "fast_fraction": 0.0, "reasons": {}}


def test_report_dict_is_serialisable_and_counts_add_up():
    sl = shortlist({"a": "TRBV20-1*01", "b": "TRAV1-2*01"},
                   {"a": "TRBJ2-1*01"}, COMBOS)
    d = sl.as_dict()
    assert d["implied"] + d["rescued"] == sl.n_total
    assert d["reasons"] == {"v_only": 1}
    import json

    json.dumps(d)                                     # must survive the run report


def test_ambiguity_lists_resolve_from_any_member():
    """`v_calls`/`j_calls` can be comma lists; a hit on any member must find the scaffold."""
    from arda.annotate.shortlist import load_combinations

    import polars as pl

    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        pl.DataFrame({"scaffold_id": ["S1"], "locus": ["TRA"],
                      "v_calls": ["TRAV1-2*01,TRAV1-2*02"], "j_calls": ["TRAJ33*01"],
                      "n_pad": ["0"]}).write_csv(fh.name, separator="\t")
        combos = load_combinations(fh.name)
    assert combos[("TRAV1-2*01", "TRAJ33*01")] == "S1"
    assert combos[("TRAV1-2*02", "TRAJ33*01")] == "S1"
    assert shortlist({"r": "TRAV1-2*02"}, {"r": "TRAJ33*01"}, combos).implied == {"r": "S1"}


def test_shortlist_dataclass_defaults_are_independent():
    """A mutable default shared between instances would corrupt a second run in one process."""
    a, b = Shortlist(), Shortlist()
    a.implied["x"] = "S"
    a.rescue.append("y")
    assert b.implied == {} and b.rescue == []
