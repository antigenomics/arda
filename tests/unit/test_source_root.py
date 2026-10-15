"""An installed ARDA nested in another project's venv must use its own cache."""
from arda import paths


def test_source_root_rejects_unrelated_ancestor_and_honors_explicit_home(tmp_path, monkeypatch):
    monkeypatch.delenv('ARDA_HOME', raising=False)
    (tmp_path / 'database').mkdir()
    (tmp_path / 'pyproject.toml').write_text('[project]\nname="tcren"\n')
    installed = tmp_path / '.venv/lib/python/site-packages/arda/paths.py'
    monkeypatch.setattr(paths, '__file__', str(installed))
    paths._source_root.cache_clear()
    try:
        assert paths._source_root() is None
        monkeypatch.setenv('ARDA_HOME', str(tmp_path))
        paths._source_root.cache_clear()
        assert paths._source_root() == tmp_path
        monkeypatch.delenv('ARDA_HOME')
        monkeypatch.setattr(paths, '__file__', str(tmp_path / 'src/arda/paths.py'))
        paths._source_root.cache_clear()
        assert paths._source_root() == tmp_path
    finally:
        paths._source_root.cache_clear()
