"""IgBLAST discovery and auto-fetch (offline).

Lives in the DB-free suite deliberately. The bug these cover is that ``arda igblast`` was
unreachable from a plain ``pip install``, and a test that skips without a reference database
would have skipped on exactly the installation where the bug lives.
"""

import pytest

from arda import _igblast_fetch, igblast


def _fake_release(root, *, complete=True):
    """A directory shaped like an installed IgBLAST release."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "igblastn").write_text("#!/bin/sh\n")
    (root / "internal_data" / "human").mkdir(parents=True, exist_ok=True)
    (root / "internal_data" / "human" / "human_V.nin").write_text("")
    if complete:
        (root / _igblast_fetch.VERSION_FILE).write_text("1.22.0\n")
    return root


@pytest.mark.parametrize(
    "system, expected",
    [("Linux", "x64-linux"), ("Darwin", "x64-macosx"), ("Windows", "x64-win64")],
)
def test_platform_suffix(monkeypatch, system, expected):
    monkeypatch.delenv("ARDA_IGBLAST_ASSET", raising=False)
    monkeypatch.setattr(_igblast_fetch.platform, "system", lambda: system)
    assert _igblast_fetch.platform_suffix() == expected


def test_platform_suffix_env_override(monkeypatch):
    monkeypatch.setenv("ARDA_IGBLAST_ASSET", "x64-linux")
    monkeypatch.setattr(_igblast_fetch.platform, "system", lambda: "SunOS")
    assert _igblast_fetch.platform_suffix() == "x64-linux"


def test_explicit_root_wins_and_is_validated(monkeypatch, tmp_path):
    root = _fake_release(tmp_path / "mine")
    monkeypatch.setenv("ARDA_IGBLAST", str(root))
    igblast.igblast_root.cache_clear()
    assert igblast.igblast_root() == root
    assert igblast.tool("igblastn") == root / "igblastn"
    assert igblast.igdata_env()["IGDATA"] == str(root)
    assert igblast.has_internal_annotation("human", "IG")
    assert not igblast.has_internal_annotation("human", "TR")

    monkeypatch.setenv("ARDA_IGBLAST", str(tmp_path / "empty"))
    igblast.igblast_root.cache_clear()
    with pytest.raises(igblast.IgBlastError, match="does not contain igblastn"):
        igblast.igblast_root()
    igblast.igblast_root.cache_clear()


def test_a_source_checkout_bin_is_used_before_fetching(monkeypatch, tmp_path):
    """setup.sh's layout must keep working — and must not trigger a download."""
    root = _fake_release(tmp_path / "bin")
    monkeypatch.delenv("ARDA_IGBLAST", raising=False)
    monkeypatch.setattr(igblast, "bin_dir", lambda: root)
    monkeypatch.setattr(igblast, "fetch", lambda *a, **k: pytest.fail("must not fetch"))
    igblast.igblast_root.cache_clear()
    assert igblast.igblast_root() == root
    igblast.igblast_root.cache_clear()


def test_readiness_is_the_marker_not_the_executable(tmp_path):
    """The gate must not be `igblastn.exists()`.

    Two shipped bugs in this repo had exactly that shape: an artifact whose presence gated its
    own completeness became visible mid-build, and every other process reported success over a
    half-written tree. Here that would mean an IgBLAST with executables but no `internal_data`,
    whose failure mode is the misleading "ships no internal annotation for organism 'human'".
    """
    partial = _fake_release(tmp_path / "partial", complete=False)
    assert (partial / "igblastn").exists()
    assert not _igblast_fetch.is_complete(partial)
    (partial / _igblast_fetch.VERSION_FILE).write_text("1.22.0\n")
    assert _igblast_fetch.is_complete(partial)


def test_fetch_is_a_no_op_when_already_complete(tmp_path, monkeypatch):
    root = _fake_release(tmp_path / "igblast")
    monkeypatch.setattr(_igblast_fetch, "find_tarball", lambda s: pytest.fail("must not download"))
    assert _igblast_fetch.fetch(root) == root


def test_no_auto_fetch_refuses_loudly(tmp_path, monkeypatch):
    """Air-gapped runs must get an error naming the fix, not a silent empty result."""
    monkeypatch.setenv("ARDA_NO_AUTO_FETCH", "1")
    with pytest.raises(RuntimeError, match="ARDA_NO_AUTO_FETCH"):
        _igblast_fetch.fetch(tmp_path / "igblast")


def test_staging_is_a_sibling_of_the_destination(tmp_path, monkeypatch):
    """Assert WHERE the staging directory lives, not merely that the result is complete.

    `shutil.move` across a filesystem boundary silently degrades to a recursive copy, which is
    how the reference fetch broke before — and a same-filesystem test cannot reproduce that.
    What is testable, and what actually prevents it, is that staging is a sibling of the
    destination, so the final move is always a rename within one filesystem.
    """
    dest = tmp_path / "cache" / "igblast"
    seen = {}

    def fake_lay_out(release_root, staging):
        seen["staging"] = staging
        _fake_release(staging, complete=False)

    monkeypatch.setattr(_igblast_fetch, "find_tarball", lambda s: "ncbi-igblast-1.22.0-x64-linux.tar.gz")
    monkeypatch.setattr(_igblast_fetch, "_download", lambda url, d, **k: d.write_bytes(b""))
    monkeypatch.setattr(_igblast_fetch, "_safe_extract", lambda t, i: i)
    monkeypatch.setattr(_igblast_fetch, "lay_out", fake_lay_out)

    assert _igblast_fetch.fetch(dest) == dest
    assert seen["staging"].parent == dest.parent, "staging must be a sibling of the destination"
    assert _igblast_fetch.is_complete(dest)
    assert _igblast_fetch.installed_version(dest) == "1.22.0"
    assert not seen["staging"].exists(), "staging must not survive a successful install"
