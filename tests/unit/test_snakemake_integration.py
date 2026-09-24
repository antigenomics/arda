"""The Snakemake workflow's Python prelude, exercised without Snakemake installed.

`integrations/snakemake/arda/Snakefile` is not importable -- `rule` is Snakemake syntax, not
Python -- but everything that can actually be *wrong* lives above the first rule: reading the
sheet through arda rather than through a second copy of the parsing rules, turning arda's
`ValueError` into Snakemake's `WorkflowError`, and asking arda for the regime's flags. So the
prelude is sliced off and executed here.

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
    exec(compile(src[:cut], str(SNAKEFILE), "exec"), ns)  # noqa: S102 -- our own file
    return ns


@pytest.fixture
def sheet(tmp_path):
    """Write a sheet whose FASTQs exist: arda validates its inputs at parse time."""
    def make(rows, name="s.tsv", reads_dir=None):
        delim = "," if name.endswith(".csv") else "\t"
        d = tmp_path / reads_dir if reads_dir else tmp_path
        d.mkdir(parents=True, exist_ok=True)
        out = []
        for sid, r1, r2 in rows:
            for f in (r1, r2):
                if f:
                    (d / f).write_text("@r\nACGT\n+\nIIII\n")
            pre = f"{reads_dir}/" if reads_dir else ""
            out.append((sid, f"{pre}{r1}" if r1 else "", f"{pre}{r2}" if r2 else ""))
        p = tmp_path / name
        head = delim.join(["sample", "fastq_1", "fastq_2"])
        p.write_text(head + "\n" + "".join(delim.join(r) + "\n" for r in out))
        return p
    return make


def test_the_prelude_merges_repeated_samples_in_row_order(sheet, tmp_path):
    s = sheet([("A", "a0_1.fq", "a0_2.fq"),
               ("A", "a1_1.fq", "a1_2.fq"),
               ("B", "b0_1.fq", "b0_2.fq")])
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
    assert ns["SAMPLES"] == ["A", "B"]
    assert ns["READ_GROUPS"]["A"] == [(str(tmp_path / "a0_1.fq"), str(tmp_path / "a0_2.fq")),
                                      (str(tmp_path / "a1_1.fq"), str(tmp_path / "a1_2.fq"))]
    assert ns["READ_GROUPS"]["B"] == [(str(tmp_path / "b0_1.fq"), str(tmp_path / "b0_2.fq"))]


def test_part_names_are_zero_padded_so_sorted_is_numeric(sheet, tmp_path):
    """`arda cluster reduce` merges parts in NAME order and that order IS read order."""
    ns = _prelude({"samples": str(sheet([("A", "x_1.fq", "")])), "outdir": str(tmp_path)})
    names = [Path(ns["part"]("A", i)).name for i in (2, 10)]
    assert names == ["shard_00002.airr.tsv", "shard_00010.airr.tsv"]
    assert sorted(names) == names


def test_a_blank_fastq_2_is_single_end_not_an_empty_path(sheet, tmp_path):
    ns = _prelude({"samples": str(sheet([("A", "x_1.fq", "")])), "outdir": str(tmp_path)})
    assert ns["READ_GROUPS"]["A"] == [(str(tmp_path / "x_1.fq"), None)]


def test_relative_paths_resolve_against_the_sheet(sheet, tmp_path):
    s = sheet([("A", "x_1.fq", "x_2.fq")], reads_dir="reads")
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
    assert ns["READ_GROUPS"]["A"] == [(str(tmp_path / "reads/x_1.fq"),
                                       str(tmp_path / "reads/x_2.fq"))]


def test_csv_is_read_as_csv(sheet, tmp_path):
    s = sheet([("A", "x_1.fq", "x_2.fq")], name="s.csv")
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path)})
    assert ns["READ_GROUPS"]["A"] == [(str(tmp_path / "x_1.fq"), str(tmp_path / "x_2.fq"))]


def test_the_flags_come_from_arda_not_from_a_restatement(sheet, tmp_path):
    """Never: the workflow must not carry its own copy of the regime presets.

    The two speed levers do not compose, and any denoising preset but `fast` reads a Stage-1
    column that only `--junction-quality` writes. A second copy of that knowledge drifts, and the
    symptom is a silent 2-4x slowdown or a gate that never runs.
    """
    cfg = {"samples": str(sheet([("A", "x_1.fq", "")])), "outdir": str(tmp_path),
           "regime": "amplicon", "map_threads": 16, "reduce_threads": 4}
    ns = _prelude(cfg)
    assert ns["MAP_FLAGS"] == regime_flags("amplicon", threads=16)[0]
    assert ns["REDUCE_FLAGS"] == regime_flags("amplicon", threads=4)[1]
    assert "--junction-quality" in ns["MAP_FLAGS"]
    assert "--ec-mode amplicon" in ns["REDUCE_FLAGS"]


def test_the_sheet_is_parsed_by_arda_so_the_cli_rules_apply_here_too(sheet, tmp_path):
    """Never: the workflow must not carry its own copy of the SHEET rules either.

    The prelude used to re-implement the parse, and the copy accepted a sample whose read groups
    mix paired with single-end -- which `arda rnaseq --samples` refuses, because half a sample's
    reads silently losing their mate is not a configuration. One parser, one answer.
    """
    s = sheet([("A", "p_1.fq", "p_2.fq"), ("A", "s_1.fq", "")])
    with pytest.raises(_WorkflowError, match="mixes paired and single-end"):
        _prelude({"samples": str(s), "outdir": str(tmp_path)})


def test_a_missing_sheet_or_a_bad_regime_is_refused_up_front(sheet, tmp_path):
    with pytest.raises(_WorkflowError, match="samples="):
        _prelude({})
    with pytest.raises(_WorkflowError, match="regime must be"):
        _prelude({"samples": str(sheet([("A", "x_1.fq", "")])), "outdir": str(tmp_path),
                  "regime": "bulk"})
    with pytest.raises(_WorkflowError, match="empty sample id"):
        _prelude({"samples": str(sheet([("", "x_1.fq", "")])), "outdir": str(tmp_path)})


def test_an_empty_sheet_is_refused(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("sample\tfastq_1\tfastq_2\n")
    with pytest.raises(_WorkflowError, match="no rows"):
        _prelude({"samples": str(p), "outdir": str(tmp_path)})


# --- the nf-core/airrflow samplesheet ---------------------------------------------------------
# One sheet has to drive this workflow, the Nextflow module and `arda cluster` alike. The prelude
# does not parse it -- `arda.samples.read_sheet` does -- so these assert the WIRING: that the
# airrflow dialect arrives intact, and that a per-sample `species` reaches the flags.

def _airrflow_sheet(tmp_path, rows):
    cols = ["sample_id", "subject_id", "species", "pcr_target_locus", "tissue", "sex", "age",
            "biomaterial_provider", "single_cell", "filename_R1", "filename_R2"]
    body = []
    for r in rows:
        for key in ("filename_R1", "filename_R2"):
            if r.get(key):
                (tmp_path / r[key]).write_text("@r\nACGT\n+\nIIII\n")
        body.append("\t".join(str(r.get(c, "")) for c in cols))
    p = tmp_path / "airrflow.tsv"
    p.write_text("\t".join(cols) + "\n" + "\n".join(body) + "\n")
    return p


def _r(**kw):
    base = dict(sample_id="S1", subject_id="D1", species="human", pcr_target_locus="ig",
                tissue="blood", sex="F", age="1", biomaterial_provider="lab",
                single_cell="FALSE", filename_R1="a_1.fq", filename_R2="a_2.fq")
    base.update(kw)
    return base


def test_an_airrflow_samplesheet_drives_the_workflow_unchanged(tmp_path):
    s = _airrflow_sheet(tmp_path, [_r(), _r(sample_id="S2", filename_R1="b_1.fq",
                                       filename_R2="b_2.fq")])
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
    assert ns["SAMPLES"] == ["S1", "S2"]
    assert ns["READ_GROUPS"]["S1"] == [(str(tmp_path / "a_1.fq"), str(tmp_path / "a_2.fq"))]


def test_the_organism_comes_from_the_sheet_per_sample_not_from_one_constant(tmp_path):
    """A cohort may legitimately mix organisms; `--config organism=` is only the fallback."""
    s = _airrflow_sheet(tmp_path, [_r(species="human"),
                                   _r(sample_id="M1", species="mouse",
                                      filename_R1="m_1.fq", filename_R2="m_2.fq")])
    ns = _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
    assert ns["ORGANISM_OF"] == {"S1": "human", "M1": "mouse"}
    assert "--organism mouse" in ns["flags_for"]("M1", 0, 8)
    assert "--organism human" in ns["flags_for"]("S1", 0, 8)


def test_a_sheet_without_species_falls_back_to_the_config(sheet, tmp_path):
    ns = _prelude({"samples": str(sheet([("A", "x_1.fq", "")])), "outdir": str(tmp_path),
                   "organism": "mouse"})
    assert ns["ORGANISM_OF"] == {"A": "mouse"}


def test_a_single_cell_sheet_is_refused_rather_than_folded_into_one_repertoire(tmp_path):
    """Never: `arda cells` takes a per-molecule UMI consensus, not a read pair."""
    s = _airrflow_sheet(tmp_path, [_r(single_cell="TRUE")])
    with pytest.raises(_WorkflowError, match="single_cell=TRUE"):
        _prelude({"samples": str(s), "outdir": str(tmp_path / "out")})
