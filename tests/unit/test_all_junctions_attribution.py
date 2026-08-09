"""Under ``--all-junctions``, an ASSEMBLED row outranks the read's own truncated Stage-1 row.

⛔ ``_assign_coverage``'s pass 1 walks the concatenated frame in order -- mapped rows first -- and
stops at the first key it finds. With ``complete_only=False`` the read's own TRUNCATED junction is
itself a clonotype key, so it won that race, ``done`` then blocked the assembled row, and the
contig's clonotype was emitted with ``duplicate_count`` 0 while its reads were credited to a
truncated PREFIX of its own junction.

The contig's complete junction is the better evidence -- it is what the assembly was for -- so the
read is routed there, and a clonotype left with no reads at all is not emitted.

⚠ Neither half touches read conservation: the reads move between clonotypes, and the dropped row has
no reads by definition.
"""

from __future__ import annotations

import polars as pl

from arda.rnaseq.correct import correct_airr

FULL = "TGTGCGAAAGGGGCCCTTCAGAAAACATTACGTTTGGGGGAGTCTATACCCCTAAATCCTTTTGATGTCTGG"
TRUNC = FULL[:45]                      # what the read reached on its own: incomplete


def _row(sid, jn, aa):
    return {"sequence_id": sid, "sequence": FULL, "junction": jn, "junction_aa": aa,
            "v_call": "IGHV3-23*01", "j_call": "IGHJ3*01", "locus": "IGH"}


def _run(tmp_path, complete_only):
    pl.DataFrame([_row("r1/1", TRUNC, "CAKGALQ"), _row("r2/1", TRUNC, "CAKGALQ")]
                 ).write_csv(tmp_path / "m.tsv", separator="\t")
    pl.DataFrame([_row("r1/1", FULL, "CAKGALQKW"), _row("r2/1", FULL, "CAKGALQKW")]
                 ).write_csv(tmp_path / "a.tsv", separator="\t")
    out = tmp_path / f"o{complete_only}.tsv"
    rep = correct_airr(tmp_path / "m.tsv", out, map_d=False,
                       extra_airr=tmp_path / "a.tsv", complete_only=complete_only)
    rows = pl.read_csv(out, separator="\t", infer_schema_length=0).to_dicts()
    return rep, {r["junction_aa"]: int(r["duplicate_count"]) for r in rows}


def test_the_contig_junction_gets_the_reads_not_its_truncated_prefix(tmp_path):
    rep, got = _run(tmp_path, complete_only=False)
    assert got.get("CAKGALQKW") == 2, f"the contig's clonotype should hold both reads, got {got}"
    assert "CAKGALQ" not in got, f"a zero-read truncated clonotype was emitted: {got}"
    assert rep.reads_assigned == 2


def test_complete_only_is_unaffected(tmp_path):
    """The default path never built the truncated clonotype, so both agree."""
    rep, got = _run(tmp_path, complete_only=True)
    assert got == {"CAKGALQKW": 2}
    assert rep.reads_assigned == 2


def test_no_read_is_lost_either_way(tmp_path):
    """⛔ The reads MOVE between clonotypes; they are never discarded."""
    for co in (True, False):
        rep, got = _run(tmp_path, complete_only=co)
        assert rep.reads_assigned == sum(got.values()) == 2, (co, got)
