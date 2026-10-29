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


def test_airr_reader_normalizes_nt_but_preserves_protein(tmp_path):
    from arda.annotate.airr_out import read_airr
    path = tmp_path / 'input.tsv'
    pl.DataFrame({'sequence': ['NRYSWKM'], 'sequence_alignment': ['NRYSWKM'],
                  'germline_alignment': ['NRYSWKM']}).write_csv(path, separator='\t')
    assert read_airr(path)['sequence'][0] == 'A' * 7
    protein = read_airr(path, seqtype='aa')
    assert protein['sequence'][0] == 'NRYSWKM'
    assert read_airr(path)['germline_alignment'].equals(protein['germline_alignment'])


@pytest.mark.parametrize('sequence', ['ACGX', 'ACG?', 'ACG-', 'ACG.', 'ACGU'])
def test_raw_invalid_dna_symbols_are_rejected(sequence, tmp_path):
    with pytest.raises(ValueError, match='invalid symbol'):
        normalize_nt(sequence)
    path = tmp_path/'bad.fa'
    path.write_text(f'>r\n{sequence}\n')
    with pytest.raises(ValueError, match='invalid symbol'):
        list(read_sequences(path))
    with pytest.raises(ValueError, match='invalid symbol'):
        list(read_pairs(path))
    with pytest.raises(ValueError, match='invalid symbol'):
        list(read_records(path))
    with pytest.raises(ValueError, match='invalid symbol'):
        normalize_nt_frame(pl.DataFrame({'sequence': [sequence]}))


def test_alignment_gaps_are_explicit_and_invalid_symbols_still_reject():
    assert normalize_nt('n.-ry', aligned=True) == 'A.-AA'
    frame = pl.DataFrame({'sequence_alignment': ['n.-ry']})
    assert normalize_nt_frame(frame)['sequence_alignment'][0] == 'A.-AA'
    with pytest.raises(ValueError, match='invalid symbol'):
        normalize_nt('A-X', aligned=True)
    with pytest.raises(ValueError, match='invalid symbol'):
        normalize_nt_frame(pl.DataFrame({'sequence_alignment': ['A-X']}))


def test_direct_projection_checks_normalized_strand_before_alignment(monkeypatch):
    from arda.annotate import transfer
    seen = []
    def project(*args, **kwargs):
        seen.append(args[0])
        raise RuntimeError('projection reached')
    monkeypatch.setattr(transfer._markup, 'transfer_regions', project)
    from types import SimpleNamespace
    ref = SimpleNamespace(starts=[1], ends=[3], v_sequence_end=0, j_sequence_start=0)
    hit = dict(qaln='TTT', taln='NNN', qstart=1, qend=3, tstart=1)
    # Normalize submitted N->A before reverse-complementing: coding T is valid.
    with pytest.raises(RuntimeError, match='projection reached'):
        transfer.transfer_hit('r', 'TTT', hit, ref, rev_comp=True, submitted_seq='NNN')
    assert seen == ['TTT']
    with pytest.raises(ValueError, match='disagrees'):
        transfer.transfer_hit('r', 'AAA', dict(hit, qaln='AAA'), ref,
                              rev_comp=True, submitted_seq='NNN')
    with pytest.raises(ValueError, match='before alignment'):
        transfer.transfer_hit('r', 'NNN', dict(hit, qaln='NNN'), ref,
                              rev_comp=True, submitted_seq='NNN')
