"""A sample is declared, never inferred.

These tests exist because the alternative -- deriving the grouping from filenames -- is what every
other tool does and what this one deliberately does not. `RNA-SAMPLE_ID:12:00XX919:3_1.fastq.gz`
has no field a rule can call the sample, and a rule that guesses wrong turns one repertoire into
four without saying so.
"""

import pytest

from arda.samples import Sample, load, read_sheet


@pytest.fixture
def fq(tmp_path):
    """Make named empty FASTQs and return a path-maker."""
    def make(name: str):
        p = tmp_path / name
        p.write_text("@r1\nACGT\n+\nIIII\n")
        return p
    return make


# ---------------------------------------------------------------- option form

def test_one_pair_with_out_prefix_is_the_historical_invocation(fq):
    a, b = fq("s_1.fq"), fq("s_2.fq")
    got = load(r1=[a], r2=[b], out_prefix="S1")
    assert got == [Sample("S1", ((a, b),))]


def test_single_end_needs_no_r2(fq):
    a = fq("s_1.fq")
    got = load(r1=[a], out_prefix="S1")
    assert got == [Sample("S1", ((a, None),))]


def test_five_pairs_five_ids_are_five_samples(fq):
    pairs = [(fq(f"x{i}_1.fq"), fq(f"x{i}_2.fq")) for i in range(5)]
    got = load(r1=[p[0] for p in pairs], r2=[p[1] for p in pairs],
               ids=[f"S{i}" for i in range(5)])
    assert [s.id for s in got] == ["S0", "S1", "S2", "S3", "S4"]
    assert all(len(s.pairs) == 1 for s in got)


def test_repeated_ids_merge_in_declared_order(fq):
    """The user's case: the first three read groups are one sample, the last two are their own."""
    pairs = [(fq(f"x{i}_1.fq"), fq(f"x{i}_2.fq")) for i in range(5)]
    got = load(r1=[p[0] for p in pairs], r2=[p[1] for p in pairs],
               ids=["A", "A", "A", "B", "C"])
    assert [s.id for s in got] == ["A", "B", "C"]
    assert got[0].pairs == tuple(pairs[:3])          # declared order, never sorted
    assert got[1].pairs == (pairs[3],)
    assert got[2].pairs == (pairs[4],)


def test_interleaved_ids_still_group_by_first_appearance(fq):
    pairs = [(fq(f"x{i}_1.fq"), fq(f"x{i}_2.fq")) for i in range(4)]
    got = load(r1=[p[0] for p in pairs], r2=[p[1] for p in pairs], ids=["A", "B", "A", "B"])
    assert [s.id for s in got] == ["A", "B"]
    assert got[0].pairs == (pairs[0], pairs[2])
    assert got[1].pairs == (pairs[1], pairs[3])


def test_more_than_one_pair_without_ids_is_refused_not_guessed(fq):
    a, b = fq("RNA-S:12:00XX919:1_1.fq"), fq("RNA-S:12:00XX919:2_1.fq")
    with pytest.raises(ValueError, match="will not guess"):
        load(r1=[a, b], out_prefix="S")


def test_mismatched_option_counts_name_the_counts(fq):
    a, b, c = fq("a_1.fq"), fq("a_2.fq"), fq("b_1.fq")
    with pytest.raises(ValueError, match=r"--r2 given 1 time\(s\) against 2 --r1"):
        load(r1=[a, c], r2=[b], ids=["A", "B"])
    with pytest.raises(ValueError, match=r"--id given 1 time\(s\) against 2 --r1"):
        load(r1=[a, c], ids=["A"])


def test_out_prefix_and_id_together_are_ambiguous(fq):
    with pytest.raises(ValueError, match="redundant and ambiguous"):
        load(r1=[fq("a_1.fq")], ids=["A"], out_prefix="B")


def test_a_missing_file_names_the_option_index(fq, tmp_path):
    with pytest.raises(ValueError, match=r"--r1 #2: no such file"):
        load(r1=[fq("a_1.fq"), tmp_path / "gone.fq"], ids=["A", "B"])


def test_an_id_that_would_become_a_path_is_refused(fq):
    with pytest.raises(ValueError, match="becomes an output filename"):
        load(r1=[fq("a_1.fq")], ids=["a/b"])


def test_a_sample_may_not_mix_paired_and_single_end(fq):
    a, b, c = fq("a_1.fq"), fq("a_2.fq"), fq("c_1.fq")
    # positionally matched, so an empty R2 slot cannot be expressed -- build it through a sheet
    with pytest.raises(ValueError, match="mixes paired and single-end"):
        _sheet_load(a, b, c)


def _sheet_load(a, b, c):
    sheet = a.parent / "mixed.tsv"
    sheet.write_text(f"sample\tfastq_1\tfastq_2\nS\t{a}\t{b}\nS\t{c}\t\n")
    return read_sheet(sheet)


# ----------------------------------------------------------------- sheet form

def test_tsv_sheet_merges_repeated_samples(fq):
    rows = [(fq(f"l{i}_1.fq"), fq(f"l{i}_2.fq")) for i in range(4)]
    sheet = rows[0][0].parent / "s.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\n"
                     + "".join(f"A\t{a}\t{b}\n" for a, b in rows[:3])
                     + f"B\t{rows[3][0]}\t{rows[3][1]}\n")
    got = read_sheet(sheet)
    assert [s.id for s in got] == ["A", "B"]
    assert got[0].pairs == tuple(rows[:3])


def test_csv_sheet_and_nf_core_extra_columns_are_accepted(fq, caplog):
    a, b = fq("a_1.fq"), fq("a_2.fq")
    sheet = a.parent / "s.csv"
    sheet.write_text(f"sample,fastq_1,fastq_2,strandedness,seq_platform\n"
                     f"CONTROL_REP1,{a},{b},auto,ILLUMINA\n")
    got = read_sheet(sheet)
    assert got == [Sample("CONTROL_REP1", ((a, b),))]
    assert "strandedness, seq_platform" in caplog.text


def test_relative_paths_resolve_against_the_sheet_not_the_cwd(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    for n in ("a_1.fq", "a_2.fq"):
        (data / n).write_text("@r\nACGT\n+\nIIII\n")
    sheet = data / "s.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\nA\ta_1.fq\ta_2.fq\n")
    monkeypatch.chdir(tmp_path)              # cwd is NOT where the files are
    got = read_sheet(sheet)
    assert got[0].pairs == ((data / "a_1.fq", data / "a_2.fq"),)


def test_a_typo_in_the_r2_header_is_an_error_not_a_single_end_run(fq):
    a, b = fq("a_1.fq"), fq("a_2.fq")
    sheet = a.parent / "s.tsv"
    sheet.write_text(f"sample\tfastq1\tfastq2\nA\t{a}\t{b}\n")
    with pytest.raises(ValueError, match="missing fastq_1"):
        read_sheet(sheet)


def test_sheet_is_exclusive_with_the_options(fq):
    a = fq("a_1.fq")
    sheet = a.parent / "s.tsv"
    sheet.write_text(f"sample\tfastq_1\nA\t{a}\n")
    with pytest.raises(ValueError, match="exclusive"):
        load(r1=[a], ids=["A"], sheet=sheet)
    with pytest.raises(ValueError, match="cannot name them all"):
        load(sheet=sheet, out_prefix="X")


def test_an_empty_sheet_is_refused(tmp_path):
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\n")
    with pytest.raises(ValueError, match="no rows"):
        read_sheet(sheet)



# ── the two label columns ─────────────────────────────────────────────────────────────────────

def test_project_and_batch_are_read_as_labels_and_no_longer_warned(tmp_path, caplog):
    """They change nothing about a run; `arda qc batch` groups a cohort by them."""
    for name in ("a_1.fq", "a_2.fq", "b_1.fq", "b_2.fq"):
        (tmp_path / name).write_text("")
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfastq_1\tfastq_2\tproject\tbatch\n"
                     "PT01\ta_1.fq\ta_2.fq\tTRIAL9\tRUN3\n"
                     "PT02\tb_1.fq\tb_2.fq\tTRIAL9\tRUN3\n")
    with caplog.at_level("WARNING"):
        got = read_sheet(sheet)
    assert [(s.id, s.project, s.batch) for s in got] == [("PT01", "TRIAL9", "RUN3"),
                                                         ("PT02", "TRIAL9", "RUN3")]
    assert "unknown sample-sheet column" not in caplog.text


def test_a_sheet_without_them_still_parses_and_the_labels_are_empty(tmp_path):
    (tmp_path / "a_1.fq").write_text("")
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfastq_1\nPT01\ta_1.fq\n")
    (one,) = read_sheet(sheet)
    assert (one.project, one.batch) == ("", "")


def test_a_samples_lanes_take_their_labels_from_its_first_row(tmp_path):
    """A sample's lanes belong to one batch by definition; a sheet saying otherwise has a typo."""
    for name in ("l1_1.fq", "l2_1.fq"):
        (tmp_path / name).write_text("")
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfastq_1\tbatch\n"
                     "PT01\tl1_1.fq\tRUN3\n"
                     "PT01\tl2_1.fq\tRUN4\n")
    (one,) = read_sheet(sheet)
    assert one.batch == "RUN3" and len(one.pairs) == 2
