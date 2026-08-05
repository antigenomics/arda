"""Binary discovery must reject an mmseqs that cannot use the shipped indexes.

An mmseqs index is only reusable by the release that compiled it. Taking whatever `mmseqs`
sits on PATH therefore had a silent cost: `database/`'s precompiled DBs were discarded and a
private cache rebuilt on every fresh environment, with no error and no log line. Observed on a
cluster whose `~/bin/mmseqs` (a bare-git-hash build) shadowed conda's version-matched one.

These tests use fake `mmseqs` shims -- no real binary, no network, no DB.
"""

from __future__ import annotations

import os
import stat
import textwrap

import pytest

from arda import mmseqs


@pytest.fixture(autouse=True)
def _clear_caches():
    """`mmseqs_binary` and `_version_of` are lru_cached; discovery tests must not share state."""
    mmseqs.mmseqs_binary.cache_clear()
    mmseqs._version_of.cache_clear()
    yield
    mmseqs.mmseqs_binary.cache_clear()
    mmseqs._version_of.cache_clear()


def _shim(path, version_string: str):
    """A minimal executable that answers `mmseqs version`."""
    path.write_text(textwrap.dedent(f"""\
        #!/bin/sh
        [ "$1" = version ] && echo "{version_string}" && exit 0
        exit 1
    """))
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_version_key_bridges_punctuation_but_not_releases():
    assert mmseqs.version_key("18.8cc5c") == mmseqs.version_key("18-8cc5c")
    assert mmseqs.version_key(" 18-8CC5C\n") == mmseqs.version_key("18.8cc5c")
    assert mmseqs.version_key("17-b804f") != mmseqs.version_key("18-8cc5c")


# The static release asset prints its FULL commit hash; bioconda and the committed index
# marker print release+short-commit. Release 18 is commit 8cc5c..., so these are one build.
_FULL_18 = "8cc5ce367b5638c4306c2d7cfc652dd099a4643f"
_OTHER = "76da68ad7577378410c075049e18666fcc94f8d1"   # a genuinely different build


@pytest.mark.parametrize("a,b,expected", [
    ("18.8cc5c", "18-8cc5c", True),        # conda vs static punctuation
    (_FULL_18, "18.8cc5c", True),          # full hash vs release+short -- the same build
    ("18.8cc5c", _FULL_18, True),          # and symmetric
    (_OTHER, "18.8cc5c", False),           # different commit -> different build
    ("17-b804f", "18-8cc5c", False),       # different release
])
def test_versions_compatible(a, b, expected):
    assert mmseqs.versions_compatible(a, b) is expected


def test_the_binary_arda_itself_fetches_is_not_rejected():
    """Regression: a pure `version_key` compare rejected arda's OWN auto-fetched binary.

    The osx-universal asset reports a full commit hash while the shipped index marker says
    `18.8cc5c`. Under the old comparison every macOS install would have discarded the
    precompiled index, warned, and rebuilt a private cache -- for the binary arda downloaded
    for itself.
    """
    assert mmseqs.versions_compatible(_FULL_18, "18.8cc5c")


def test_a_full_hash_candidate_is_accepted_end_to_end(tmp_path, monkeypatch):
    static = _shim(tmp_path / "mmseqs", _FULL_18)
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: None)
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18.8cc5c")
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: str(static))
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")

    def _boom():
        raise AssertionError("must not re-fetch: the PATH binary is the same build")

    monkeypatch.setattr(mmseqs, "_auto_fetch", _boom)
    assert mmseqs.mmseqs_binary() == str(static)


def test_explicit_override_wins_and_is_not_version_checked(tmp_path, monkeypatch):
    """$ARDA_MMSEQS is the user's call — never second-guess it."""
    wrong = _shim(tmp_path / "mine", "1-deadbee")
    monkeypatch.setenv("ARDA_MMSEQS", str(wrong))
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18-8cc5c")
    assert mmseqs.mmseqs_binary() == str(wrong)


def test_a_matching_path_binary_is_accepted(tmp_path, monkeypatch):
    good = _shim(tmp_path / "mmseqs", "18-8cc5c")
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: None)
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18.8cc5c")  # conda spelling
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: str(good))
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")
    assert mmseqs.mmseqs_binary() == str(good)


def test_a_mismatched_path_binary_triggers_the_fetch(tmp_path, monkeypatch):
    """The bug this exists for: PATH has an mmseqs, but not the one the index needs."""
    stale = _shim(tmp_path / "mmseqs", "76da68ad7577378410c075049e18666fcc94f8d1")
    fetched = _shim(tmp_path / "fetched", "18-8cc5c")
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.delenv("ARDA_NO_AUTO_FETCH", raising=False)
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: None)
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18-8cc5c")
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: str(stale))
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")
    monkeypatch.setattr(mmseqs, "_auto_fetch", lambda: str(fetched))

    assert mmseqs.mmseqs_binary() == str(fetched), "a mismatched PATH binary was accepted"


def test_mismatch_with_no_fetch_available_falls_back_but_warns(tmp_path, monkeypatch):
    """Degrade, but never silently — the warning must name the consequence."""
    stale = _shim(tmp_path / "mmseqs", "17-b804f")
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setenv("ARDA_NO_AUTO_FETCH", "1")
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: None)
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18-8cc5c")
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: str(stale))
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")

    with pytest.warns(RuntimeWarning, match="rebuilt"):
        assert mmseqs.mmseqs_binary() == str(stale)


def test_no_committed_index_means_no_version_filter(tmp_path, monkeypatch):
    """A plain `pip install` ships no index (the reference asset omits them).

    With nothing to be compatible with, any working mmseqs is correct — the first run
    just builds its own cache. Filtering here would force a pointless download.
    """
    whatever = _shim(tmp_path / "mmseqs", "99-zzzzz")
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: None)
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: None)
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: str(whatever))
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")

    def _boom():
        raise AssertionError("auto-fetch must not run when there is no index to match")

    monkeypatch.setattr(mmseqs, "_auto_fetch", _boom)
    assert mmseqs.mmseqs_binary() == str(whatever)


def test_the_bundled_wheel_is_preferred_over_path(tmp_path, monkeypatch):
    bundled = _shim(tmp_path / "bundled", "18-8cc5c")
    on_path = _shim(tmp_path / "mmseqs", "18-8cc5c")
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: str(bundled))
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18-8cc5c")
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: str(on_path))
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")
    assert mmseqs.mmseqs_binary() == str(bundled)


def test_nothing_found_at_all_raises_with_actionable_advice(tmp_path, monkeypatch):
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setenv("ARDA_NO_AUTO_FETCH", "1")
    monkeypatch.setattr(mmseqs, "_bundled_binary", lambda: None)
    monkeypatch.setattr(mmseqs, "committed_index_version", lambda: "18-8cc5c")
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _n: None)
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path / "nonexistent")

    # The advice has to be advice that WORKS. It used to say `pip install 'arda-mapper[mmseqs]'`,
    # and that extra pointed at a distribution never published to PyPI -- so the one actionable
    # line in the error was itself a hard failure (`No matching distribution found`). Assert both
    # halves: a route that exists, and not that one.
    with pytest.raises(mmseqs.MMseqsError) as exc:
        mmseqs.mmseqs_binary()
    msg = str(exc.value)
    assert "ARDA_NO_AUTO_FETCH" in msg and "$ARDA_MMSEQS" in msg
    assert "arda-mapper[mmseqs]" not in msg


def test_committed_index_version_never_triggers_a_reference_download(monkeypatch):
    """Discovery must not reach `database_dir()`, which auto-fetches over the network."""
    import arda.paths as paths

    def _boom(*_a, **_k):
        raise AssertionError("committed_index_version() must not resolve database_dir()")

    monkeypatch.setattr(paths, "database_dir", _boom)
    mmseqs.committed_index_version()  # must not raise


def test_version_of_returns_none_for_a_binary_that_will_not_run(tmp_path):
    assert mmseqs._version_of(str(tmp_path / "does-not-exist")) is None
