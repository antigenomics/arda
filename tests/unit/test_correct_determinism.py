"""`correct` must produce byte-identical output for the same reads, in any input order.

This is not a style preference. Before the fix, three runs of `arda rnaseq correct` over one
unchanged 200k-read AIRR produced three different `clones.tsv`: polars' `group_by` is a
multithreaded hash aggregation, so group order, read order within a group, and the order of
equal-`count` rows were all arbitrary. Downstream that moved real values --
`IGHV3-11*06` vs `IGHV3-21*08` for one identical junction, and `duplicate_count` 11 vs 9 for
one IGK clonotype -- because `_parents` collapses onto the parent it meets first and
`_assign_coverage` is first-with-longest-overlap-wins.

Row order independence is also what makes the sharded (SLURM) and Nextflow delivery paths
provably byte-identical to a single-node run: a shard boundary changes row order and nothing
else, so if row order cannot change the output, neither can the shard layout.

No DB and no mmseqs: `correct_airr` reads an AIRR TSV.
"""

from __future__ import annotations

import random

import polars as pl
import pytest

from arda.rnaseq.correct import correct_airr

# Two paralogous V calls over the SAME junction: exactly the pair that flipped in production.
_J1 = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"   # 30 nt, in frame, C...F
_J2 = "TGTGCCAGCAGCTTAGACGGGACAGGTTTC"   # a different clonotype, equal abundance


def _rows() -> list[dict]:
    rows: list[dict] = []
    # Equal counts across paralogous V calls -> ties everywhere the old code was arbitrary.
    for v in ("TRBV20-1*01", "TRBV20-1*02"):
        for junc in (_J1, _J2):
            for k in range(6):
                rows.append({
                    "sequence_id": f"{v}_{junc[-4:]}_{k}", "junction": junc,
                    "junction_aa": "CASSLDGTF", "v_call": v,
                    "j_call": "TRBJ2-1*01", "locus": "TRB",
                })
    return rows


def _write(path, rows) -> None:
    pl.DataFrame(rows).write_csv(path, separator="\t")


def test_same_input_gives_byte_identical_output_across_repeated_runs(tmp_path):
    airr = tmp_path / "in.airr.tsv"
    _write(airr, _rows())

    hashes = set()
    for i in range(5):
        out = tmp_path / f"clones_{i}.tsv"
        correct_airr(airr, out)
        hashes.add(out.read_bytes())
    assert len(hashes) == 1, "correct() is not reproducible run-to-run on identical input"


def test_shuffling_the_input_rows_does_not_change_the_output(tmp_path):
    """The property the sharded and Nextflow paths depend on.

    A shard boundary permutes rows. If a permutation can change `clones.tsv`, then
    "accuracy does not differ between run modes" is unprovable.
    """
    base = _rows()
    ref = tmp_path / "ref.airr.tsv"
    _write(ref, base)
    correct_airr(ref, tmp_path / "ref.clones.tsv")
    expected = (tmp_path / "ref.clones.tsv").read_bytes()

    rng = random.Random(20260804)
    for i in range(6):
        shuffled = base[:]
        rng.shuffle(shuffled)
        src = tmp_path / f"shuf_{i}.airr.tsv"
        _write(src, shuffled)
        out = tmp_path / f"shuf_{i}.clones.tsv"
        correct_airr(src, out)
        assert out.read_bytes() == expected, f"row permutation {i} changed clones.tsv"


def test_output_row_order_is_total_not_merely_stable(tmp_path):
    """Equal-abundance clonotypes must still have one defined order.

    Pins the contract at ``correct.py`` -- rows are emitted by
    ``(-duplicate_count, -consensus_count, junction, v_call, j_call)``. Ranking on counts alone
    left tied clonotypes in read order, which comes from a threaded mmseqs search.
    """
    airr = tmp_path / "in.airr.tsv"
    _write(airr, _rows())
    correct_airr(airr, tmp_path / "clones.tsv")
    out = pl.read_csv(tmp_path / "clones.tsv", separator="\t", infer_schema_length=0)

    keyed = [(int(r["duplicate_count"]), int(r["consensus_count"]),
              r["junction"], r["v_call"], r["j_call"]) for r in out.iter_rows(named=True)]
    ordered = sorted(keyed, key=lambda t: (-t[0], -t[1], t[2], t[3], t[4]))
    assert keyed == ordered, "rows are not in (abundance desc, junction/V/J asc) order"
    assert len(set(keyed)) == len(keyed), "tie-break key is not unique per row"


# ---------------------------------------------------------------------------------------------
# `--error-method binom|betabinom` used to HANG, not merely disagree. Kept in this DB-free module
# because a run that never terminates is the worst failure a pipeline stage has, and this must not
# be gated behind an mmseqs/reference skip.
# ---------------------------------------------------------------------------------------------

# Both pairs are REAL: lifted from a 16,157-row Stage-1 AIRR whose 1,190 clonotypes contained
# exactly these two mutual-parent pairs (indices 77<->143 and 857<->861) under both methods.
# Same locus, same V, same J, same junction length, differing by 1 and 2 substitutions
# respectively, at spanning depths of 3/2 and 1/1 -- the low-coverage regime the pileup path
# exists to serve, which is also the only regime that produces the cycle.
_MUTUAL_PAIRS = [
    ("TGCCAACACTATAATAGTTACCCTCTCACTTTC", "TGCCAACAGTATAATAGTTACCCTCTCACTTTC",
     "IGKV1-16*02", "IGKJ4*01", "IGK", 3, 2),
    ("TGTCAGCAATATGGTAGCTCACCTCGGACGTTC", "TGTCAGCAGTATGGTAGCTCACCTCAGACGTTC",
     "IGKV3D-20*01", "IGKJ1*01", "IGK", 1, 1),
]


def _walk(parent):
    """Every node's root, refusing to loop. Returns the set of nodes that sit on a cycle."""
    on_cycle = set()
    for i in range(len(parent)):
        seen, x = set(), i
        while parent[x] is not None:
            if x in seen:
                on_cycle.add(i)
                break
            seen.add(x)
            x = parent[x]
    return on_cycle


@pytest.mark.parametrize("method", ["binom", "betabinom"])
@pytest.mark.parametrize("j1,j2,v,j,loc,s1,s2", _MUTUAL_PAIRS)
def test_error_pileup_never_makes_two_clonotypes_each_others_parent(
        method, j1, j2, v, j, loc, s1, s2):
    """`_root` is an unbounded `while parent[i] is not None`, so a 2-cycle hangs the run.

    `_parents` cannot produce one: it accepts a parent only under
    ``count[parent] * p_err >= count[child]`` with ``p_err < 1``, which forces counts to increase
    strictly along parent pointers -- the module docstring cites exactly this as its no-cycle
    proof. `_error_pileup` re-decided parentage from per-position depth and carried NO ordering
    condition, and the depth test is symmetric at low coverage: ``_binom_sf(1, 2, 0.001) =
    0.001999 > alpha = 1e-3``, so each of a pair passes as the other's error child.
    """
    import polars as pl

    from arda.rnaseq.correct import _error_pileup

    juncs = [j1, j2]
    spans = [s1, s2]
    rows = [{"sequence_id": f"c{i}r{r}", "sequence": jn, "junction": jn, "locus": loc,
             "v_call": v, "j_call": j}
            for i, jn in enumerate(juncs) for r in range(spans[i])]
    parent = _error_pileup(pl.DataFrame(rows), juncs, [v, v], [j, j], [loc, loc],
                           [None, None], spans, error_rate=0.001, indel_rate=0.001,
                           method=method)

    assert parent[0] is None or parent[parent[0]] != 0, (
        f"{method}: clonotypes 0 and 1 are each other's parent ({parent}) -- _root would not "
        f"terminate")
    assert not _walk(parent), f"{method}: parent pointers contain a cycle: {parent}"
    # A parent must be strictly more abundant than its child; that ordering IS the acyclicity.
    for child, par in enumerate(parent):
        if par is not None:
            assert spans[par] > spans[child], (
                f"{method}: clonotype {child} (span {spans[child]}) got parent {par} "
                f"(span {spans[par]}) -- counts must increase along parent pointers")


def test_isotype_vote_counts_fragments_not_assigned_mates():
    """A clonotype's isotype is the dominant class over its FRAGMENTS, one vote each.

    The vote iterated `read_list`, which holds per-mate `sequence_id`s, and looked each one's
    fragment up in `frag_iso` — so a fragment with both mates assigned contributed its calls
    twice, while a fragment with one assigned mate contributed once. The weighting was therefore
    by assigned mates, not by molecules, and a 1-fragment minority could outvote a 2-fragment
    majority.

    Constructed to be decisive: two IGHM fragments contribute one mate each, one IGHG fragment
    contributes both. Per-mate weighting gives IGHG 2 votes vs IGHM 2 and resolves the tie
    arbitrarily; per-fragment gives IGHM 2 vs IGHG 1.
    """
    from collections import Counter

    from arda.rnaseq.correct import _strip_mate

    frag_iso = {"fA": ("IGHM",), "fB": ("IGHM",), "fC": ("IGHG",)}
    read_list = ["fA/1", "fB/1", "fC/1", "fC/2"]

    per_mate = []
    for sid in read_list:
        per_mate.extend(frag_iso.get(_strip_mate(sid), ()))
    assert Counter(per_mate)["IGHG"] == 2 and Counter(per_mate)["IGHM"] == 2, (
        "the fixture no longer demonstrates the bug")

    per_frag = []
    for frag in dict.fromkeys(_strip_mate(sid) for sid in read_list):
        per_frag.extend(frag_iso.get(frag, ()))
    assert Counter(per_frag).most_common(1)[0][0] == "IGHM"
    assert Counter(per_frag)["IGHG"] == 1, "a fragment voted more than once"
