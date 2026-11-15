"""Observed pair unions retain selected-base qualities and both observed tails."""
from __future__ import annotations

import random
import pytest

from arda.refbuild.translate import reverse_complement
from arda.rnaseq.map import read_pairs


def _fastq(path, sequence, quality):
    path.write_text(f'@pair\n{sequence}\n+\n{quality}\n')
    return path


@pytest.mark.parametrize('limit', [None, 1])
def test_overlap_merge_retains_the_quality_of_each_selected_base(tmp_path, limit):
    rng = random.Random(19)
    fragment = ''.join(rng.choices('ACGT', k=80))
    left, right = fragment[:55], fragment[30:]
    # Two disagreements outside the exact anchor: choose R1 at one, R2 at the other.
    right = right[:15] + next(b for b in 'ACGT' if b != right[15]) + right[16:]
    left = left[:48] + next(b for b in 'ACGT' if b != left[48]) + left[49:]
    q1 = '5' * 55
    q2_forward = 'I' * 15 + '!' + 'I' * (len(right) - 16)
    r1 = _fastq(tmp_path/'r1.fq', left, q1)
    r2 = _fastq(tmp_path/'r2.fq', reverse_complement(right), q2_forward[::-1])
    records = list(read_pairs(r1, r2, reconstruct=True, with_qual=True, limit=limit))
    expected_quality = q1[:30] + 'I' * 15 + '5' + 'I' * (len(right) - 16)
    assert records == [('pair', fragment, expected_quality)]


@pytest.mark.parametrize('limit', [None, 1])
def test_contained_mate_does_not_truncate_the_observed_r1_tail(tmp_path, limit):
    rng = random.Random(29)
    fragment = ''.join(rng.choices('ACGT', k=90))
    r1 = _fastq(tmp_path/'r1.fq', fragment, '5'*90)
    r2 = _fastq(tmp_path/'r2.fq', reverse_complement(fragment[20:70]), 'I'*50)
    assert list(read_pairs(r1, r2, reconstruct=True, with_qual=True, limit=limit)) == [
        ('pair', fragment, '5'*20 + 'I'*50 + '5'*20)]


def test_rna_defaults_and_scheduled_regime_preserve_independent_mate_opt_out():
    import inspect
    from typer.main import get_command
    from arda.cli import app
    from arda.cluster import regime_flags
    from arda.rnaseq.pipeline import run

    assert inspect.signature(run).parameters['reconstruct'].default is True
    commands = get_command(app).commands
    param = next(p for p in commands['rnaseq'].params if p.name == 'reconstruct')
    assert param.default is True and '--no-reconstruct' in param.secondary_opts
    assert '--reconstruct' in regime_flags('rnaseq')[0]
    assert '--reconstruct' not in regime_flags('rnaseq', reconstruct=False)[0]
    assert '--reconstruct' not in regime_flags('amplicon')[0]


def test_reconstruct_fasta_keeps_unknown_quality_unknown(tmp_path):
    fragment = 'ACGTACGTACGTTTGGCCAATTGGCCAAGGTTCCAAGGTTCCAAGATCGATCGATCG'
    r1, r2 = tmp_path/'r1.fa', tmp_path/'r2.fa'
    r1.write_text('>pair\n'+fragment[:40]+'\n')
    r2.write_text('>pair\n'+reverse_complement(fragment[-40:])+'\n')
    assert list(read_pairs(r1, r2, reconstruct=True, with_qual=True)) == [('pair', fragment, '')]


def test_read_group_mean_uses_exact_lengths_before_rounding():
    from arda.rnaseq.pipeline import _merge_map_reports

    parts = [dict(major_only=True, total_reads=3, read_length_min=100,
                  read_length_max=101, read_length_mean=100.3, read_length_sum=301),
             dict(major_only=True, total_reads=2, read_length_min=105,
                  read_length_max=106, read_length_mean=105.5, read_length_sum=211)]
    merged = _merge_map_reports(parts)
    assert merged['read_length_mean'] == 102.4
    assert merged['read_length_sum'] == 512
