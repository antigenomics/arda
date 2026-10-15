"""Shipped and cached indexes obey the same source/payload identity contract."""
from arda import mmseqs
from arda.annotate import mapper
from tests.unit.test_index_concurrency import _fake_createdb


def test_version_key_bridges_separators_not_versions():
    assert mmseqs.version_key('18.8cc5c') == mmseqs.version_key('18-8cc5c')
    assert mmseqs.version_key('17-b804f') != mmseqs.version_key('18-8cc5c')


def test_committed_index_requires_identity_not_legacy_timestamp(tmp_path, monkeypatch):
    monkeypatch.setattr(mapper, 'vdj_dir', lambda org: tmp_path)
    monkeypatch.setattr(mmseqs, 'version', lambda: '18-8cc5c')
    monkeypatch.setattr(mmseqs, 'createdb', _fake_createdb)
    fasta = tmp_path / 'alleles.fasta'
    fasta.write_text('>x\nACGT\n')
    db = tmp_path / 'mmseqs' / 'nt' / 'db'
    db.parent.mkdir(parents=True)
    db.write_text('LEGACY')
    (db.parent / 'VERSION').write_text('18.8cc5c')
    assert mapper._committed_index('x', 'nt', fasta) is None
    verified = mapper._createdb_atomic(fasta, db, 2)
    assert mapper._committed_index('x', 'nt', fasta) == verified
    monkeypatch.setattr(mmseqs, 'version', lambda: '18.8cc5c')
    assert mapper._committed_index('x', 'nt', fasta) == verified
    monkeypatch.setattr(mmseqs, 'version', lambda: '17-b804f')
    assert mapper._committed_index('x', 'nt', fasta) is None


def test_public_identity_hashes_exact_resolved_search_files(tmp_path, monkeypatch):
    import hashlib
    from pathlib import Path
    from arda import reference_index_identity

    monkeypatch.setattr(mapper, 'vdj_dir', lambda org: tmp_path)
    monkeypatch.setattr(mapper, 'data_dir', lambda: tmp_path / 'data')
    monkeypatch.setattr(mmseqs, 'version', lambda: '18-8cc5c')
    monkeypatch.setattr(mmseqs, 'createdb', _fake_createdb)
    fasta = tmp_path / 'alleles.aa.fasta'
    fasta.write_text('>x\nACGT\n')
    identity = reference_index_identity('human')
    db = Path(identity['database_path'])
    assert identity['source_fasta'] == str(fasta.resolve())
    assert identity['source_fasta_sha256'] == hashlib.sha256(fasta.read_bytes()).hexdigest()
    assert identity['files_sha256'] == {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in db.parent.glob('db*')
    }
    assert identity == reference_index_identity('human')
