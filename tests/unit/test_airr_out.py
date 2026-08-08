"""Unit tests for AIRR output assembly (column order, empty/None handling)."""

import pytest

from arda.annotate.airr_out import (
    write_airr, airr_header, format_rows, _format_rows_py, _format_rows_cpp,
)
from arda.annotate.transfer import AIRR_COLUMNS, _empty_record


def test_airr_header_is_tab_joined_columns():
    assert airr_header() == "\t".join(AIRR_COLUMNS)


def test_write_airr_empty_writes_header_only(tmp_path):
    out = tmp_path / "empty.tsv"
    write_airr([], out)
    assert out.read_text() == "\t".join(AIRR_COLUMNS) + "\n"


def test_write_airr_orders_columns_and_roundtrips(tmp_path):
    rec = {c: "" for c in AIRR_COLUMNS}
    rec["sequence_id"] = "q1"
    rec["v_call"] = "IGHV1-2*01"
    out = tmp_path / "one.tsv"
    write_airr([rec], out)
    lines = out.read_text().splitlines()
    assert lines[0] == "\t".join(AIRR_COLUMNS)          # stable column order
    assert lines[1].split("\t")[0] == "q1"
    assert lines[1].split("\t")[AIRR_COLUMNS.index("v_call")] == "IGHV1-2*01"


def test_format_rows_renders_none_as_empty_field():
    rec = {c: None for c in AIRR_COLUMNS}
    rec["sequence_id"] = "q1"
    text = format_rows([rec])
    assert text.endswith("\n")
    fields = text.rstrip("\n").split("\t")
    assert fields[0] == "q1"
    # A None value becomes an empty field, not the string "None".
    assert "None" not in fields


def test_format_rows_empty_input_is_empty_string():
    assert format_rows([]) == ""


def test_format_rows_missing_key_renders_as_empty_field_not_an_error():
    """The Python version used ``dict.get``, so a record that skipped ``_empty_record`` must not
    raise. The C++ port inherits this from ``PyDict_GetItem`` returning NULL, but it is a contract
    worth pinning: an earlier rewrite of this function broke exactly this case."""
    text = format_rows([{"sequence_id": "q1"}])
    fields = text.rstrip("\n").split("\t")
    assert len(fields) == len(AIRR_COLUMNS)
    assert fields[0] == "q1"
    assert all(f == "" for f in fields[1:])


def test_the_cpp_formatter_is_byte_identical_to_the_python_one():
    """``_markup.format_rows`` replaces a loop that was 6.3 % of an amplicon run's wall (0.365 s ->
    0.171 s, in-run A/B). It is only allowed to exist if it is indistinguishable from the reference
    implementation, so assert that on every shape the real pipeline produces: filled records,
    all-None records, missing keys, non-string values and non-ASCII."""
    if _format_rows_cpp is None:                      # pragma: no cover - extension not built
        pytest.skip("_markup extension not built")

    full = _empty_record("q1", "ACGT" * 30)
    full.update(locus="TRA", v_call="TRAV25*02", j_call="TRAJ23*02",
                junction="TGTGCAGGGAAAGCTTATCTTC", junction_aa="CAGGKLIF",
                mmseqs2_score=141.0, mmseqs2_identity=97.3, rev_comp="T",
                v_sequence_end=65, j_sequence_start=69)
    nones = dict.fromkeys(AIRR_COLUMNS)
    nones["sequence_id"] = "q2"
    sparse = {"sequence_id": "q3", "v_call": "IGHV1-2*01"}
    unicode_rec = _empty_record("q4", "ACGT")
    unicode_rec["cdr3_aa"] = "CAWSÉF"
    numeric = _empty_record("q5", "ACGT")
    numeric.update(mmseqs2_evalue=1e-30, fwr1_start=7, productive=True)

    for batch in ([full], [nones], [sparse], [unicode_rec], [numeric],
                  [full, nones, sparse, unicode_rec, numeric], []):
        assert _format_rows_cpp(batch, tuple(AIRR_COLUMNS)) == _format_rows_py(batch)


def test_the_cpp_string_helpers_match_the_python_ones():
    """``_common_prefix``/``_common_suffix`` are per-character Python loops called 138,065 times
    per 100 k-read amplicon run (once per V allele in ``v_anchor_prefix``, twice per read in
    ``_anchored_vj_bounds``), so they moved into ``_markup``. They decide how much of a junction a
    called germline templates -- and the Cys104 junction gate is a threshold on exactly that
    number -- so an off-by-one here silently changes which junctions arda emits."""
    import random

    from arda.annotate import transfer as T

    if T._common_prefix is T._common_prefix_py:       # pragma: no cover - extension not built
        pytest.skip("_markup extension not built")

    random.seed(1)
    cases = [("", ""), ("A", ""), ("", "A"), ("ACGT", "ACGT"), ("TGT", "TGC")]
    for _ in range(5000):
        cases.append((
            "".join(random.choice("ACGT") for _ in range(random.randint(0, 25))),
            "".join(random.choice("ACGT") for _ in range(random.randint(0, 25))),
        ))
    for a, b in cases:
        assert T._common_prefix(a, b) == T._common_prefix_py(a, b), (a, b)
        assert T._common_suffix(a, b) == T._common_suffix_py(a, b), (a, b)
