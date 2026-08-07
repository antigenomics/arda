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


def _row(query: str, target: str, bits: float,
         tstart: int = 1, tend: int = 60) -> str:
    """One mmseqs hit line; only query/target/bits matter to the tie-break."""
    vals = {
        "query": query, "target": target, "qstart": 1, "qend": 60,
        "tstart": tstart, "tend": tend,
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


def test_a_target_inverted_row_is_not_a_usable_hit(tmp_path):
    """`tstart > tend` must be dropped, not projected.

    mmseqs can express a minus-strand nt hit on the TARGET side. arda tests only the query
    side (`rev = qstart > qend`), so such a row reads as forward, and `transfer_regions` walks
    the target forward from an inverted `tstart` -- sliding the scaffold's cdr3/fwr4 markup by
    `(tlen + 1 - tstart) - tstart` nt and taking the junction window off Cys104.

    Measured consequence, arda 2.5.6 on Jurkat (ERR3003543, a MONOCLONAL T-cell line): the
    window slid 10 nt into V framework 3, opened on a spurious TGT, closed on TGG, satisfied
    `assemble._CANON` (`^C...[FW]$`) and shipped as a 7,408-read clonotype -- 48 % the size of
    the true clone, from which it had taken 5,758 reads.

    The row is not a recoverable reflection: on that read `identity(qaln, taln)` is 0.232
    against the row's own reported `pident` of 91.3, so it is a garbage alignment. Refusing it
    is the only outcome that cannot produce a well-formed junction that is wrong.
    """
    inverted = _row("read1", "TRBV12-3*01_TRBJ1-2*01", 243.0, tstart=180, tend=1)

    assert _best_hits(_write(tmp_path / "inv.tsv", [inverted])) == {}, (
        "a target-inverted alignment row was returned as a usable hit; projecting it slides "
        "the junction window off Cys104 and emits a phantom clonotype"
    )


def test_an_inverted_row_never_beats_a_forward_one(tmp_path):
    """Refusal must survive the tie-break -- an inverted row outscoring a real hit is exactly
    the case where keeping it does damage."""
    good = _row("read1", "TRBV12-3*01_TRBJ1-2*01", 200.0)
    inverted = _row("read1", "TRBV6-6*01_TRBJ1-2*01", 243.0, tstart=180, tend=1)

    for name, rows in (("inv_first", [inverted, good]), ("good_first", [good, inverted])):
        got = _best_hits(_write(tmp_path / f"{name}.tsv", rows))["read1"]
        assert got["target"] == "TRBV12-3*01_TRBJ1-2*01", name
        assert got["tstart"] <= got["tend"], name
