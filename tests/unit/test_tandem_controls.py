"""Issue #157: the tandem search control survives public and scheduler boundaries."""
import inspect

import pytest
import typer.main

from arda.cli import app
from arda.cluster import regime_flags


def test_every_d_cli_command_exposes_independent_tandem_control():
    def check(command):
        options = {flag for param in command.params for flag in
                   (*getattr(param, 'opts', ()), *getattr(param, 'secondary_opts', ()))}
        if '--map-d' in options:
            assert {'--map-dd', '--no-map-dd'} <= options, command.name
        for child in getattr(command, 'commands', {}).values():
            check(child)

    check(typer.main.get_command(app))


@pytest.mark.parametrize('regime', ['rnaseq', 'amplicon'])
@pytest.mark.parametrize('map_d', [False, True])
def test_cluster_disables_tandem_in_both_stages(regime, map_d):
    flags = regime_flags(regime, map_d=map_d, map_dd=False)
    for stage in flags:
        assert '--no-map-dd' in stage
        assert ('--no-map-d' in stage.split()) == (not map_d)


def test_public_annotation_defaults_preserve_tandem():
    from arda.adapter import annotate_sequences
    from arda.annotate.mapper import annotate_file, annotate_records
    from arda.annotate.dmap import map_d_junction
    for function in (annotate_sequences, annotate_file, annotate_records, map_d_junction):
        assert inspect.signature(function).parameters['map_dd'].default is True


def test_native_batch_matches_scalar_d_selection_and_exclusion():
    import random
    from arda.annotate import transfer as t

    rng = random.Random(157)
    for _ in range(100):
        query = ''.join(rng.choices('ACGT', k=60))
        germlines = [(str(i), ''.join(rng.choices('ACGT', k=20))) for i in range(12)]
        germlines.append(('tied', germlines[0][1]))
        hits = [(di, score, qs-1, qe-1, ds-1, de-1)
                for _, di, score, qs, qe, ds, de, _ in
                t._markup.d_candidates_batch([query], [s for _, s in germlines], 1)]
        first = t._best_d(query, germlines, 4)
        assert t._best_d(query, germlines, 4, alignments=hits) == first
        if first:
            excluded = (first[3], first[4])
            assert t._best_d(query, germlines, 4, exclude=excluded, alignments=hits) == t._best_d(
                query, germlines, 4, exclude=excluded)


@pytest.mark.parametrize('map_d,map_dd', [(True, True), (True, False), (False, False)])
def test_pipeline_finish_d_controls_reach_clonotypes(tmp_path, map_d, map_dd):
    from pathlib import Path
    import polars as pl
    from arda.rnaseq.pipeline import finish

    airr = Path(__file__).resolve().parents[2] / 'examples/dd.airr.tsv'
    finish(airr, tmp_path, 'sample', map_d=map_d, map_dd=map_dd, write_qc=False)
    clones = pl.read_csv(tmp_path / 'sample.clones.tsv', separator='\t', quote_char=None).fill_null('')
    assert clones.height > 0
    if not map_d:
        assert 'd_call' not in clones.columns and 'd2_call' not in clones.columns
    elif not map_dd:
        assert clones.filter(pl.col('d2_call') != '').height == 0
        assert clones.filter(pl.col('d_call') != '').height > 0
    else:
        assert clones.filter(pl.col('d2_call') != '').height > 0
