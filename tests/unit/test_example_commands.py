"""Example regeneration must use installed commands and the current markup schema."""

import importlib.util
from pathlib import Path

import polars as pl
import typer.main
from typer.testing import CliRunner

from arda.cli import app


def test_junction_example_uses_supported_markup_options(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[2] / "examples/regenerate.py"
    spec = importlib.util.spec_from_file_location("regenerate_examples", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls = []
    monkeypatch.setattr(module, "EX", tmp_path)
    monkeypatch.setattr(module, "_sh", lambda *args: calls.append(args))
    module.build_junctions()
    available = {flag for p in typer.main.get_command(app).commands["markup"].params for flag in p.opts}
    assert {arg for arg in calls[0][4:] if arg.startswith("--")} <= available


def test_markup_uses_the_same_unquoted_tsv_dialect_as_other_writers(tmp_path):
    source, out = tmp_path / "input.tsv", tmp_path / "out.tsv"
    source.write_text("cdr3\tv\tj\tspecies\nCASSPLGQAYEQYF\tTRBV5-1*01\tTRBJ2-7*01\thuman\n")
    result = CliRunner().invoke(app, ["markup", "-i", str(source), "-o", str(out)])
    assert result.exit_code == 0, result.output
    frame = pl.read_csv(out, separator="\t", quote_char=None)
    assert frame["errors"][0] in (None, "")
    assert frame["cdr3fix"][0].startswith('{"cdr3":')
