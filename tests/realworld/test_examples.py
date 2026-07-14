"""The committed ``examples/`` artifacts must still reproduce, and README must not lie.

`examples/example.airr.tsv` sat stale for four release rounds — written with 49 AIRR columns
while the schema grew to 83, and carrying a `d_call` a later E-value gate correctly withdrew —
because nothing regenerated it. This is that check.

The alignment-quality scalars (`mmseqs2_score`/`_evalue`/`_identity`) are excluded from the
comparison: they move with the mmseqs2 version, and nothing in the docs quotes them. Every
biological column — calls, coordinates, regions, `d_support` — is compared exactly.

When this fails after an intentional change, run ``python examples/regenerate.py`` and read
the diff. The diff is the behaviour change.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from arda.annotate.io import read_sequences
from arda.annotate.mapper import annotate_records
from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

EXAMPLES = Path(__file__).resolve().parent.parent.parent / "examples"
VOLATILE = {"mmseqs2_score", "mmseqs2_evalue", "mmseqs2_identity"}


def _reannotate(fasta: Path) -> pl.DataFrame:
    recs = annotate_records(list(read_sequences(fasta)), "human", "nt", threads=4)
    return pl.DataFrame([{k: ("" if v is None else str(v)) for k, v in r.items()} for r in recs])


def _committed(tsv: Path) -> pl.DataFrame:
    return pl.read_csv(tsv, separator="\t", infer_schema_length=0)


def _assert_reproduces(fasta: Path, tsv: Path) -> None:
    got, want = _reannotate(fasta), _committed(tsv)
    assert got.height == want.height, f"{tsv.name}: row count moved"
    missing = set(want.columns) - set(got.columns)
    assert not missing, f"{tsv.name}: committed columns arda no longer emits: {sorted(missing)}"
    added = set(got.columns) - set(want.columns)
    assert not added, (f"{tsv.name}: arda emits {sorted(added)} but the committed artifact "
                       f"predates them -- run `python examples/regenerate.py`")
    for col in want.columns:
        if col in VOLATILE:
            continue
        g = got[col].cast(pl.Utf8).fill_null("").to_list()
        w = want[col].cast(pl.Utf8).fill_null("").to_list()
        assert g == w, f"{tsv.name}: column {col!r} drifted\n  now: {g}\n  was: {w}"


def test_example_airr_still_reproduces():
    _assert_reproduces(EXAMPLES / "example.fasta", EXAMPLES / "example.airr.tsv")


def test_dd_airr_still_reproduces():
    _assert_reproduces(EXAMPLES / "dd.fasta", EXAMPLES / "dd.airr.tsv")


def test_readme_table_for_example_airr():
    """The five rows examples/README.md tabulates, including the IGH D that is NOT called."""
    df = _committed(EXAMPLES / "example.airr.tsv")
    row = {r["sequence_id"]: r for r in df.iter_rows(named=True)}
    assert row["AF043995.1"]["d_call"] == "TRBD1*01"
    assert row["AF043995.1"]["d_support"] == "0.0465"
    # IGH: the best D scores E = 0.286, above the 0.2 gate. A call is not evidence.
    assert not row["PZ235980.1"]["d_call"], "the IGH near-miss must stay a no-call"
    for vj in ("PV083657.1", "PP196642.1", "PQ177856.1"):
        assert not row[vj]["d_call"], "VJ loci have no D gene"


def test_readme_table_for_dd_airr():
    """Both human tandem D-D reads, with the np partition the README claims is exact."""
    df = _committed(EXAMPLES / "dd.airr.tsv")
    row = {r["sequence_id"]: r for r in df.iter_rows(named=True)}

    trd = row["AM408133.1"]
    assert (trd["locus"], trd["d_call"], trd["d2_call"]) == ("TRD", "TRDD2*01", "TRDD3*01")
    assert trd["junction_aa"] == "CALGPRPSYSEELGDTHRADKLIF"
    assert (trd["np1"], trd["np2"], trd["np3"]) == ("CCCCGG", "AGCGAGGAGT", "CCATCGGG")

    igh = row["PX612894.1"]
    assert igh["d_call"] == "IGHD3-9*01"
    assert igh["d2_call"] == "IGHD2/OR15-2a*01,IGHD2/OR15-2b*01", "byte-identical D genes tie"

    for r in (trd, igh):
        seq = r["sequence"]
        d1 = seq[int(r["d_sequence_start"]) - 1 : int(r["d_sequence_end"])]
        d2 = seq[int(r["d2_sequence_start"]) - 1 : int(r["d2_sequence_end"])]
        interior = seq[int(r["v_sequence_end"]) : int(r["j_sequence_start"]) - 1]
        assert r["np1"] + d1 + r["np2"] + d2 + r["np3"] == interior, r["sequence_id"]


def test_dd_reads_reproduce_from_the_junction_alone():
    """Read space (mmseqs -> projection) and junction space (anchors) must agree."""
    from arda.annotate.dmap import map_d_junction

    for r in _committed(EXAMPLES / "dd.airr.tsv").iter_rows(named=True):
        call = map_d_junction(r["junction"], r["v_call"], r["j_call"], "human")
        assert (call.d_call, call.d2_call) == (r["d_call"], r["d2_call"]), r["sequence_id"]


def test_committed_rnaseq_clonotypes_still_reproduce(tmp_path):
    """`arda rnaseq run` on the committed 1035-read FASTQ, compared on a stable projection.

    Row order ties on (duplicate_count, consensus_count), so compare as a sorted set of the
    columns the README quotes rather than row-for-row.
    """
    pytest.importorskip("seqtree")  # `correct` needs the optional rnaseq extra
    from arda.rnaseq.assemble import assemble_contigs
    from arda.rnaseq.correct import correct_airr
    from arda.rnaseq.map import map_rnaseq

    airr, extra, clones = (tmp_path / n for n in ("m.airr.tsv", "a.airr.tsv", "c.tsv"))
    map_rnaseq(EXAMPLES / "rnaseq" / "reads.fq.gz", airr, organism="human", threads=4)
    assemble_contigs(airr, extra, organism="human", threads=4)
    correct_airr(airr, clones, organism="human", extra_airr=extra)

    cols = ["junction_aa", "v_call", "j_call", "locus", "d_call", "d2_call"]
    got = _committed(clones).select(cols).sort(cols)
    want = _committed(EXAMPLES / "rnaseq" / "clones.tsv").select(cols).sort(cols)
    assert got.to_dicts() == want.to_dicts()


def test_rnaseq_example_actually_demonstrates_d():
    """An example that shows no D would document nothing. Lock in what the README claims."""
    df = _committed(EXAMPLES / "rnaseq" / "clones.tsv")
    # 21 -> 20 once 3'-truncated V alleles stopped building scaffolds: reads on a dropped allele
    # reassign to a surviving sibling of the same gene, merging rows that shared a junction and
    # differed only in v_call (CARDIGAGGFGDNFYFFYYMDVW was split across IGHV3-9*01 / IGHV3-43*02).
    # Mapping is untouched -- 925/1035 reads, identical per-locus counts -- only the V calls moved.
    assert df.height == 20
    with_d = df.filter(pl.col("d_call").is_not_null() & (pl.col("d_call") != ""))
    assert with_d.height == 8, "8 of 20 clonotypes carry a d_call"
    dd = df.filter(pl.col("d2_call").is_not_null() & (pl.col("d2_call") != ""))
    assert dd.height == 1 and dd["d_call"][0] == "TRDD2*01" and dd["d2_call"][0] == "TRDD3*01"
    assert set(df["locus"].to_list()) == {"IGH", "IGK", "IGL", "TRA", "TRB", "TRD"}


@pytest.mark.parametrize("path", [
    "example.fasta", "example.airr.tsv", "dd.fasta", "dd.airr.tsv",
    "junctions.tsv", "junctions.markup.tsv", "junctions.report.txt",
    "rnaseq/reads.fq.gz", "rnaseq/clones.tsv", "regenerate.py", "README.md",
])
def test_every_documented_artifact_is_committed(path):
    assert (EXAMPLES / path).exists(), f"examples/{path} is referenced but missing"
