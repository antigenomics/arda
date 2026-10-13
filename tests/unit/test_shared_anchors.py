"""Shared V alleles retain every locus declared by the shipped reference."""
import csv

import pytest

from arda import cdr3fix


@pytest.mark.parametrize('organism', ['human', 'mouse', 'rat', 'rabbit', 'rhesus_monkey'])
def test_shared_anchor_loci_and_candidates(organism, tmp_path, monkeypatch):
    path = cdr3fix._anchor_path(organism)
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    memberships = {}
    for row in rows:
        memberships.setdefault((row['segment'], row['allele']), set()).add(row['locus'])
    shared = {key: loci for key, loci in memberships.items() if loci == {'TRA', 'TRD'}}
    if organism in ('human', 'mouse'):
        assert len(shared) == (15 if organism == 'human' else 42)
    reversed_path = tmp_path / 'anchors.tsv'
    with reversed_path.open('w') as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0], delimiter='\t')
        writer.writeheader()
        writer.writerows(reversed(rows))
    try:
        for source in (path, reversed_path):
            monkeypatch.setattr(cdr3fix, '_anchor_path', lambda org: source)
            cdr3fix.load_anchors.cache_clear()
            anchors = cdr3fix.load_anchors(organism)
            for key, loci in memberships.items():
                assert set(anchors[key].loci) == loci, key
            for locus in ('TRA', 'TRD'):
                owners = cdr3fix._run_owners(anchors, 'V', locus)
                candidates = {allele for values in owners.values() for _, allele in values}
                usable = {allele for (seg, allele) in shared
                          if anchors[(seg, allele)].status == 'ok'
                          and anchors[(seg, allele)].functionality != 'P'}
                assert usable <= candidates, usable - candidates
    finally:
        cdr3fix.load_anchors.cache_clear()


@pytest.mark.parametrize('v_call', ['TRAV38-2/DV8*01', 'TRAV13-4/DV7*01'])
def test_shared_v_defers_to_j_for_chain(v_call):
    assert cdr3fix.resolve_locus(v_call, 'TRAJ1*01') == 'TRA'
    assert cdr3fix.resolve_locus(v_call, 'TRDJ1*01') == 'TRD'
    assert cdr3fix.resolve_locus(v_call) == 'TRA'


@pytest.mark.parametrize('organism,v_call,j_call', [
    ('human', 'TRAV38-2', 'TRDJ1*01'),
    ('human', 'DV8', 'TRDJ1*01'),
    ('mouse', 'TRAV13-4', 'TRDJ1*01'),
    ('mouse', 'DV7', 'TRDJ1*01'),
])
def test_shared_alias_chain_across_consumers(organism, v_call, j_call):
    from arda.annotate.dmap import map_d_junction
    from arda.scenarios import germlines_for
    anchors = cdr3fix.load_anchors(organism)
    v_id = cdr3fix.resolve_allele(v_call, 'V', anchors)
    va, ja = anchors[('V', v_id)], anchors[('J', j_call)]
    aa = va.templated_aa + 'GG' + ja.templated_aa
    nt = va.germline_nt + 'GGGGGG' + ja.germline_nt
    assert cdr3fix.markup_cdr3(aa, v_call, j_call, organism).locus == 'TRD'
    germlines = germlines_for(v_call, j_call, organism)
    assert germlines is not None and germlines.locus == 'TRD'
    assert germlines.d_germlines
    assert map_d_junction(nt, v_call, j_call, organism).locus == 'TRD'
