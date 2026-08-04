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
