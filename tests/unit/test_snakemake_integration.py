"""The Snakemake workflow's Python prelude, exercised without Snakemake installed.

`integrations/snakemake/arda/Snakefile` is not importable -- `rule` is Snakemake syntax, not
Python -- but everything that can actually be *wrong* lives above the first rule: parsing the
sheet, merging repeated sample ids in row order, resolving relative paths, and asking arda for the
regime's flags rather than restating them. So the prelude is sliced off and executed here.

Without this the integration is covered by nothing: it ships in a directory no test touches, which
is exactly how the Nextflow module sat four releases stale.
"""

from pathlib import Path

import pytest

from arda.cluster import regime_flags

SNAKEFILE = Path(__file__).resolve().parents[2] / "integrations/snakemake/arda/Snakefile"


class _WorkflowError(Exception):
    pass


def _prelude(config):
    """Execute the Snakefile up to its first rule, with Snakemake's globals faked in."""
    src = SNAKEFILE.read_text()
    cut = src.index("\nrule ")
    ns = {"config": config, "WorkflowError": _WorkflowError, "expand": lambda *a, **k: []}
    exec(compile(src[:cut], str(SNAKEFILE), "exec"), ns)  # noqa: S102 — our own file
    return ns


@pytest.fixture
def sheet(tmp_path):
    def make(rows, name="s.tsv"):
        delim = "," if name.endswith(".csv") else "\t"
        p = tmp_path / name
        head = delim.join(["sample", "fastq_1", "fastq_2"])
        p.write_text(head + "\n" + "".join(delim.join(r) + "\n" for r in rows))
        return p
    return make


def test_the_prelude_merges_repeated_samples_in_row_order(sheet, tmp_path):
    s = sheet([("A", "/d/a0_1.fq", "/d/a0_2.fq"),
               ("A", "/d/a1_1.fq", "/d/a1_2.fq"),
               ("B", "/d/b0_1.fq", "/d/b0_2.fq")])
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
    assert ns["SAMPLES"] == ["A", "B"]
    assert ns["READ_GROUPS"]["A"] == [("/d/a0_1.fq", "/d/a0_2.fq"),
                                      ("/d/a1_1.fq", "/d/a1_2.fq")]
    assert ns["READ_GROUPS"]["B"] == [("/d/b0_1.fq", "/d/b0_2.fq")]


def test_part_names_are_zero_padded_so_sorted_is_numeric(sheet, tmp_path):
    """`arda cluster reduce` merges parts in NAME order and that order IS read order."""
    ns = _prelude({"samples": str(sheet([("A", "/d/x_1.fq", "")])), "outdir": str(tmp_path)})
    names = [Path(ns["part"]("A", i)).name for i in (2, 10)]
    assert names == ["shard_00002.airr.tsv", "shard_00010.airr.tsv"]
    assert sorted(names) == names


def test_a_blank_fastq_2_is_single_end_not_an_empty_path(sheet, tmp_path):
    ns = _prelude({"samples": str(sheet([("A", "/d/x_1.fq", "")])), "outdir": str(tmp_path)})
    assert ns["READ_GROUPS"]["A"] == [("/d/x_1.fq", None)]


def test_relative_paths_resolve_against_the_sheet(sheet, tmp_path):
    s = sheet([("A", "reads/x_1.fq", "reads/x_2.fq")])
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
    assert ns["READ_GROUPS"]["A"] == [(str(s.parent / "reads/x_1.fq"),
                                       str(s.parent / "reads/x_2.fq"))]


def test_csv_is_read_as_csv(sheet, tmp_path):
    s = sheet([("A", "/d/x_1.fq", "/d/x_2.fq")], name="s.csv")
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path)})
    assert ns["READ_GROUPS"]["A"] == [("/d/x_1.fq", "/d/x_2.fq")]


def test_the_flags_come_from_arda_not_from_a_restatement(sheet, tmp_path):
    """Never: the workflow must not carry its own copy of the regime presets.

    The two speed levers do not compose, and any denoising preset but `fast` reads a Stage-1
    column that only `--junction-quality` writes. A second copy of that knowledge drifts, and the
    symptom is a silent 2-4x slowdown or a gate that never runs.
    """
    cfg = {"samples": str(sheet([("A", "/d/x_1.fq", "")])), "outdir": str(tmp_path),
           "regime": "amplicon", "map_threads": 16, "reduce_threads": 4}
    ns = _prelude(cfg)
    assert ns["MAP_FLAGS"] == regime_flags("amplicon", threads=16)[0]
    assert ns["REDUCE_FLAGS"] == regime_flags("amplicon", threads=4)[1]
    assert "--junction-quality" in ns["MAP_FLAGS"]
    assert "--ec-mode amplicon" in ns["REDUCE_FLAGS"]


def test_a_missing_sheet_or_a_bad_regime_is_refused_up_front(sheet, tmp_path):
    with pytest.raises(_WorkflowError, match="samples="):
        _prelude({})
    with pytest.raises(_WorkflowError, match="regime must be"):
        _prelude({"samples": str(sheet([("A", "/d/x_1.fq", "")])), "outdir": str(tmp_path),
                  "regime": "bulk"})
    with pytest.raises(_WorkflowError, match="needs both"):
        _prelude({"samples": str(sheet([("", "/d/x_1.fq", "")])), "outdir": str(tmp_path)})


def test_an_empty_sheet_is_refused(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("sample\tfastq_1\tfastq_2\n")
    with pytest.raises(_WorkflowError, match="no rows"):
        _prelude({"samples": str(p), "outdir": str(tmp_path)})
