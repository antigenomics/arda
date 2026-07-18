"""The precompiled mmseqs index shipped in ``database/`` is used only if it matches the running mmseqs.

The match is on the version *marker*, and it was an exact string ``==``. But ``mmseqs version`` prints
the same release+commit with different punctuation across builds -- the official static binary says
``18-8cc5c``, the bioconda build ``18.8cc5c``. So the moment arda's toolchain moved from conda's
mmseqs to the static one, every committed marker (written ``18.8cc5c``) stopped matching the runtime
(``18-8cc5c``): the shipped indexes were silently never used and every run rebuilt a private cache.

These pin the fix: a separator-insensitive key bridges the cosmetic difference without accepting a
genuinely different version.
"""

from pathlib import Path

from arda import mmseqs
from arda.annotate import mapper


def test_version_key_bridges_separators_not_versions():
    assert mmseqs.version_key("18.8cc5c") == mmseqs.version_key("18-8cc5c")
    assert mmseqs.version_key("  18-8CC5C\n") == "18-8cc5c"          # strip + lowercase
    assert mmseqs.version_key("17-b804f") != mmseqs.version_key("18-8cc5c")


def _fake_index(tmp_path: Path, marker: str) -> Path:
    d = tmp_path / "mmseqs" / "nt"
    d.mkdir(parents=True)
    (d / "db").write_text("index")
    (d / "VERSION").write_text(marker + "\n")
    return tmp_path


def test_committed_index_accepts_a_separator_variant_marker(tmp_path, monkeypatch):
    """A ``18.8cc5c`` marker must satisfy a ``18-8cc5c`` runtime -- the exact case the bug rejected."""
    monkeypatch.setattr(mapper, "vdj_dir", lambda org: _fake_index(tmp_path, "18.8cc5c"))
    monkeypatch.setattr(mapper.mmseqs, "version", lambda: "18-8cc5c")
    assert mapper._committed_index("x", "nt") is not None


def test_committed_index_rejects_a_different_version(tmp_path, monkeypatch):
    monkeypatch.setattr(mapper, "vdj_dir", lambda org: _fake_index(tmp_path, "17-b804f"))
    monkeypatch.setattr(mapper.mmseqs, "version", lambda: "18-8cc5c")
    assert mapper._committed_index("x", "nt") is None
