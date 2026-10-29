"""Sequencing DNA IUPAC is always replaced with A, before downstream work."""
import polars as pl
import pytest

from arda.annotate.io import normalize_nt, normalize_nt_frame, read_sequences
from arda.rnaseq.map import read_pairs
from arda.singlecell import read_records


def test_every_dna_iupac_and_provenance_columns():
    assert normalize_nt('acgtRYSWKMBDHVNryswkmbdhvn') == 'ACGT' + 'A' * 22
    frame = pl.DataFrame({'sequence': ['acgtn'], 'junction': ['TGNNNA'],
                          'junction_aa': ['CX'], 'sequence_alignment': ['aN-g'],
                          'germline_alignment': ['aN-g'], 'sequence_id': ['NRY'],
                          'junction_quality': ['NRYSWK']})
    out = normalize_nt_frame(frame)
    assert out['sequence'][0] == 'ACGTA'
    assert out['junction'][0] == 'TGAAAA'
    assert out['junction_aa'][0] == '*K'
    assert out['sequence_alignment'][0] == 'AA-G'
    for col in ('germline_alignment', 'sequence_id', 'junction_quality'):
        assert out[col].equals(frame[col])


@pytest.mark.parametrize('fmt', ['fasta', 'fastq'])
def test_file_and_single_cell_ingress(tmp_path, fmt):
    path = tmp_path / ('reads.fa' if fmt == 'fasta' else 'reads.fq')
    path.write_text('>r cD:i:2\nacgtnry\n' if fmt == 'fasta' else '@r cD:i:2\nacgtnry\n+\nNRYSWKM\n')
    assert list(read_sequences(path))[0] == ('r', 'ACGTAAA')
    assert list(read_records(path))[0][1] == 'ACGTAAA'
    assert list(read_pairs(path))[0] == ('r', 'ACGTAAA')
    assert list(read_pairs(path, limit=1))[0] == ('r', 'ACGTAAA')
    assert list(read_sequences(path, seqtype='aa'))[0][1] == 'acgtnry'
    if fmt == 'fastq':
        assert list(read_sequences(path, with_qual=True))[0][2] == 'NRYSWKM'


@pytest.mark.parametrize('seqtype,expected', [('nt', 'ACGTAAA'), ('aa', 'acgtnry')])
def test_annotation_normalizes_before_search(tmp_path, monkeypatch, seqtype, expected):
    from arda.annotate import mapper
    seen = []
    monkeypatch.setattr(mapper, '_prep', lambda *args: (None, None, 1, 7, None))
    def createdb(query, *args, **kwargs):
        seen.append(query.read_text())
        raise RuntimeError('search boundary reached')
    monkeypatch.setattr(mapper.mmseqs, 'createdb', createdb)
    with pytest.raises(RuntimeError, match='search boundary reached'):
        mapper.annotate_records([('r', 'acgtnry')], seqtype=seqtype)
    assert seen == [f'>r\n{expected}\n']


def test_paired_and_direct_cell_ingress(tmp_path):
    from arda.singlecell import assemble_cell
    from arda.rnaseq.map import merge_pair
    r1, r2 = tmp_path/'r1.fq', tmp_path/'r2.fq'
    for path, seq in ((r1, 'acgtnry'), (r2, 'rynnacg')):
        path.write_text(f'@r\n{seq}\n+\nIIIIIII\n')
    expected = [('r/1', 'ACGTAAA'), ('r/2', 'AAAAACG')]
    assert list(read_pairs(r1, r2)) == expected
    assert list(read_pairs(r1, r2, limit=1)) == expected
    ambiguous = 'ACGTNRYACGTACGTACGT'
    normal = normalize_nt(ambiguous)
    kw = dict(k=3, min_overlap=3, min_molecules=1)
    assert assemble_cell('cell', [ambiguous], **kw) == assemble_cell('cell', [normal], **kw)
    assert merge_pair(ambiguous, 'ACGTACGTACGT', min_overlap=8) == merge_pair(normal, 'ACGTACGTACGT', min_overlap=8)
