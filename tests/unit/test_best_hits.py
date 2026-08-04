"""`_best_hits` must resolve equal-bit-score ties deterministically.

Two paralogous scaffolds can tie exactly on whole-scaffold bit score. Whichever one wins
decides the emitted `v_call`, so if the winner depends on the order mmseqs happened to list
the alignments in, arda's output is not reproducible -- and every byte-identity claim about
the sharded/Nextflow delivery paths is unfounded, because a shard boundary changes nothing
about the tie but everything about the row order feeding it.

No DB, no mmseqs: `_best_hits` reads a TSV, so the tie is expressible directly.
"""

from __future__ import annotations

from pathlib import Path

from arda import mmseqs
from arda.annotate.mapper import _best_hits

_COLS = mmseqs.DEFAULT_FORMAT_OUTPUT.split(",")


def _row(query: str, target: str, bits: float) -> str:
    """One mmseqs hit line; only query/target/bits matter to the tie-break."""
    vals = {
        "query": query, "target": target, "qstart": 1, "qend": 60, "tstart": 1, "tend": 60,
        "qlen": 60, "tlen": 400, "alnlen": 60, "mismatch": 0, "gapopen": 0, "cigar": "60M",
        "qaln": "A" * 60, "taln": "A" * 60, "evalue": 1e-20, "bits": bits, "pident": 100.0,
    }
    return "\t".join(str(vals[c]) for c in _COLS)


def _write(path: Path, rows: list[str]) -> Path:
    path.write_text("\n".join(rows) + "\n")
    return path


def test_equal_bit_scores_resolve_the_same_way_in_either_input_order(tmp_path):
    """The same tie, listed in both orders, must pick the same target."""
    a = _row("read1", "IGHV3-30*01_IGHJ4*02", 88.0)
    b = _row("read1", "IGHV3-30-3*01_IGHJ4*02", 88.0)

    forward = _best_hits(_write(tmp_path / "fwd.tsv", [a, b]))["read1"]["target"]
    reverse = _best_hits(_write(tmp_path / "rev.tsv", [b, a]))["read1"]["target"]

    assert forward == reverse, (
        f"tied hits resolved to {forward!r} vs {reverse!r} depending on input order; "
        "the best-hit tie-break is not deterministic"
    )


def test_a_strictly_higher_score_still_wins_regardless_of_order(tmp_path):
    """The tie-break must not disturb ordinary max-bits selection."""
    lo = _row("read1", "IGHV1-2*01_IGHJ4*02", 70.0)
    hi = _row("read1", "IGHV3-30*01_IGHJ4*02", 91.5)

    for name, rows in (("lo_first", [lo, hi]), ("hi_first", [hi, lo])):
        got = _best_hits(_write(tmp_path / f"{name}.tsv", rows))["read1"]
        assert got["target"] == "IGHV3-30*01_IGHJ4*02", name
        assert got["bits"] == 91.5, name


def test_ties_are_broken_independently_per_query(tmp_path):
    """Several queries, each with its own tie, all resolve stably."""
    rows_fwd, rows_rev = [], []
    for i in range(25):
        q = f"read{i}"
        x = _row(q, f"IGKV1-{i}*01_IGKJ1*01", 60.0)
        y = _row(q, f"IGKV3-{i}*01_IGKJ1*01", 60.0)
        rows_fwd += [x, y]
        rows_rev += [y, x]

    fwd = _best_hits(_write(tmp_path / "many_fwd.tsv", rows_fwd))
    rev = _best_hits(_write(tmp_path / "many_rev.tsv", rows_rev))

    assert {q: h["target"] for q, h in fwd.items()} == {q: h["target"] for q, h in rev.items()}
