"""AIRR ``umi_count``: DISTINCT UMIs, not distinct records.

The spec wording is what this column is: *"Number of distinct UMIs represented by this
sequence."* That is why migec's conditional ``.c<k>`` / ``.<m>`` suffixes cannot block it -- they
subdivide the reads that already share one ``<umi>``, so a UMI counted once is counted once
whether it arrived as one record, as several contigs, or as a split pair.

The shapes here are taken from real ``migec assemble`` output (2026-09-24), not invented:

    S.AAAACCCC.1   split barcode, molecule 1   ``.c<k>`` omitted (components == 1)
    S.AAAACCCC.2   split barcode, molecule 2   ``.<m>`` written (molecules > components)
    S.GGGGTTTT     one molecule                both suffixes omitted

`project/design-singlecell.md` S4 has the measurement and why contig mode does *not* produce the
inequality (two contigs of one molecule never overlap, so at most one carries a junction).
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.rnaseq.correct import correct_airr

_JUNCTION = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"   # 30 nt, in frame, C...F


def _write(path, ids: list[str], junction: str = _JUNCTION) -> None:
    pl.DataFrame([{
        "sequence_id": sid, "junction": junction, "junction_aa": "CASSLDGTF",
        "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01", "locus": "TRB",
    } for sid in ids]).write_csv(path, separator="\t", quote_style="never")


def _clones(tmp_path, ids, **kw) -> pl.DataFrame:
    tmp_path.mkdir(parents=True, exist_ok=True)
    airr = tmp_path / "in.airr.tsv"
    _write(airr, ids)
    out = tmp_path / "clones.tsv"
    correct_airr(airr, out, organism="human", map_d=False, **kw)
    return pl.read_csv(out, separator="\t", quote_char=None)


# The split barcode: three records, two UMIs.
_SPLIT = ["S.AAAACCCC.1", "S.AAAACCCC.2", "S.GGGGTTTT"]


def test_a_split_barcode_counts_one_umi_not_two_records(tmp_path):
    """The whole point: `consensus_count` counts records, `umi_count` counts UMIs."""
    df = _clones(tmp_path, _SPLIT, cell_from="migec")
    row = df.row(0, named=True)
    assert row["consensus_count"] == 3
    assert row["umi_count"] == 2


def test_the_column_is_omitted_without_a_dialect(tmp_path):
    """QC trap 4 -- a metric with no input is OMITTED, never written as 0 or 1.

    A bulk run must not read as "one UMI per clonotype", which is what a default-on column of 1s
    would say, and which nothing downstream could distinguish from a real measurement.
    """
    assert "umi_count" not in _clones(tmp_path, _SPLIT).columns


@pytest.mark.parametrize("kw", [
    {"cell_from": "cellranger"},                    # names a cell, never a UMI
    {"cell_from": "prefix"},                        # same
    {"cell_regex": r"^(?P<cell>[A-Z0-9]+)\."},      # a regex defines `cell` and nothing else
])
def test_dialects_that_name_no_umi_omit_the_column(tmp_path, kw):
    assert "umi_count" not in _clones(tmp_path, _SPLIT, **kw).columns


def test_identifiers_that_do_not_parse_omit_the_column(tmp_path):
    """`--cell-from migec` pointed at a bulk run is a mistake, not a 1-UMI cohort."""
    ids = ["SRR123.1", "SRR123.2", "SRR123.3"]
    assert "umi_count" not in _clones(tmp_path, ids, cell_from="migec").columns


def test_umi_count_never_exceeds_consensus_count(tmp_path):
    """A UMI cannot be represented by more UMIs than there are records carrying it."""
    ids = ["S.AAAACCCC.1", "S.AAAACCCC.2", "S.GGGGTTTT", "S.TTTTAAAA", "S.CCCCGGGG.c1"]
    row = _clones(tmp_path, ids, cell_from="migec").row(0, named=True)
    assert row["umi_count"] == 4 <= row["consensus_count"] == 5


def test_the_same_umi_under_two_cells_is_two_umis(tmp_path):
    """Never key on the UMI alone.

    A UMI is short and a barcode space saturates -- that is what migec's occupancy reporting is
    for -- so the same UMI string under two cells is two molecules, and a bare-UMI key would
    silently merge them.
    """
    ids = ["S.AAACCCGGGTTTAAAC.ACGTACGT", "S.TTTGGGCCCAAATTTG.ACGTACGT"]
    row = _clones(tmp_path, ids, cell_from="migec").row(0, named=True)
    assert row["umi_count"] == 2


def test_umi_count_is_appended_last(tmp_path):
    """Never inserted into the DataFrame literal.

    The literal is followed by conditional appends, so a key added to it moves `d_call` whenever
    `umi_count` is on -- the same trap `map.py` names for `cell_id`. A consumer reading the
    shipped set by position must be unaffected by whether this column is present.
    """
    plain = _clones(tmp_path / "a", _SPLIT, flag_chimeras=True)
    withumi = _clones(tmp_path / "b", _SPLIT, cell_from="migec", flag_chimeras=True)
    assert withumi.columns[:len(plain.columns)] == plain.columns
    assert withumi.columns[-1] == "umi_count"
