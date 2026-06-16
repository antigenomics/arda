"""Unit tests for AIRR output assembly (column order, empty/None handling)."""

from arda.annotate.airr_out import write_airr, airr_header, format_rows
from arda.annotate.transfer import AIRR_COLUMNS


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
