"""CLI output contract on shipped real tandem-D examples."""
from pathlib import Path

import polars as pl
from typer.testing import CliRunner

from arda.cli import app
from tests.conftest import requires_mmseqs, requires_human_db


@requires_mmseqs
@requires_human_db
def test_cli_tandem_switch_preserves_other_annotation(tmp_path):
    fasta = Path(__file__).resolve().parents[2] / 'examples/dd.fasta'
    baseline = None
    for name, flags in [('default', []), ('single', ['--no-map-dd']), ('none', ['--no-map-d'])]:
        target = tmp_path / f'{name}.tsv'
        result = CliRunner().invoke(app, ['annotate', '-i', str(fasta), '-o', str(target), *flags])
        assert result.exit_code == 0, result.output
        rows = pl.read_csv(target, separator='\t', quote_char=None, infer_schema_length=0).fill_null('')
        stable = rows.select(c for c in rows.columns if not c.startswith(('d_', 'd2_', 'np')))
        if baseline is None:
            baseline = stable
            assert rows.filter(pl.col('d2_call') != '').height > 0
        assert stable.equals(baseline)
        if name != 'default':
            assert rows.filter(pl.col('d2_call') != '').height == 0
        if name == 'none':
            assert rows.filter(pl.col('d_call') != '').height == 0
