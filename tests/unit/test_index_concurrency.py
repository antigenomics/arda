"""Index generations must remain coherent through concurrent publication and source changes."""
import multiprocessing as mp
from pathlib import Path

import pytest

from arda.annotate import mapper


def _fake_createdb(fasta, out, dbtype=2):
    out = Path(out)
    for suffix in ('', '.index', '.dbtype', '_h', '_h.index', '_h.dbtype'):
        Path(str(out) + suffix).write_bytes(Path(fasta).read_bytes() + suffix.encode())


def _build(args):
    fasta, db = map(Path, args)
    mapper.mmseqs.createdb = _fake_createdb
    mapper.mmseqs.version = lambda: '18-8cc5c'
    return str(mapper._createdb_atomic(fasta, db, 2))


def test_concurrent_builders_publish_one_immutable_generation(tmp_path):
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>x\nACGT\n')
    db = tmp_path / 'cache' / 'db'
    with mp.Pool(4) as pool:
        paths = pool.map(_build, [(str(fasta), str(db))] * 4)
    assert len(set(paths)) == 1
    assert Path(paths[0]).read_bytes() == fasta.read_bytes()
    assert not list(db.parent.glob('.db.lock'))
    assert not list(db.parent.glob('generations/.build-*'))


def test_source_identity_overrules_mtime_and_preserves_old_readers(tmp_path, monkeypatch):
    monkeypatch.setattr(mapper.mmseqs, 'createdb', _fake_createdb)
    monkeypatch.setattr(mapper.mmseqs, 'version', lambda: '18-8cc5c')
    fasta, db = tmp_path / 'ref.fa', tmp_path / 'cache' / 'db'
    fasta.write_text('>x\nOLD\n')
    old = mapper._createdb_atomic(fasta, db, 2)
    # A legacy cache, even newer than the source, remains untouched and unused.
    db.write_text('STALE')
    Path(str(db) + '.99').write_text('STALE SPLIT')
    fasta.write_text('>x\nNEW\n')
    new = mapper._createdb_atomic(fasta, db, 2)
    assert old != new
    assert old.read_text() == '>x\nOLD\n'
    assert new.read_text() == '>x\nNEW\n'
    assert db.read_text() == 'STALE'
    assert not Path(str(new) + '.99').exists()


@pytest.mark.parametrize('damage', ['missing', 'corrupt', 'extra'])
def test_payload_integrity_is_validated_before_reuse(tmp_path, monkeypatch, damage):
    monkeypatch.setattr(mapper.mmseqs, 'createdb', _fake_createdb)
    monkeypatch.setattr(mapper.mmseqs, 'version', lambda: '18-8cc5c')
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>x\nACGT\n')
    db = tmp_path / 'cache' / 'db'
    old = mapper._createdb_atomic(fasta, db, 2)
    sidecar = Path(str(old) + '_h.index')
    if damage == 'missing':
        sidecar.unlink()
    elif damage == 'corrupt':
        sidecar.write_text('corrupt')
    else:
        Path(str(old) + '.0').write_text('unexpected split')
    new = mapper._createdb_atomic(fasta, db, 2)
    assert new != old
    assert mapper._createdb_atomic(fasta, db, 2) == new


def test_split_payloads_and_source_change_during_build(tmp_path, monkeypatch):
    monkeypatch.setattr(mapper.mmseqs, 'version', lambda: '18-8cc5c')
    fasta, db = tmp_path / 'ref.fa', tmp_path / 'cache' / 'db'
    fasta.write_text('>x\nACGT\n')
    def split_build(fasta, out, dbtype=2):
        _fake_createdb(fasta, out, dbtype)
        Path(out).rename(Path(str(out) + '.0'))
        Path(str(out) + '_h').rename(Path(str(out) + '_h.0'))
    monkeypatch.setattr(mapper.mmseqs, 'createdb', split_build)
    result = mapper._createdb_atomic(fasta, db, 2)
    assert Path(str(result) + '.0').exists()
    assert mapper._createdb_atomic(fasta, db, 2) == result
    def changing_build(fasta, out, dbtype=2):
        _fake_createdb(fasta, out, dbtype)
        fasta.write_text('>x\nCHANGED\n')
    monkeypatch.setattr(mapper.mmseqs, 'createdb', changing_build)
    with pytest.raises(RuntimeError, match='changed during index build'):
        mapper._createdb_atomic(fasta, db, 2, force=True)
    assert Path(str(result) + '.0').exists()
