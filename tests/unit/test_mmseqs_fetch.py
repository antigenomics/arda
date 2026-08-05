"""Unit tests for mmseqs binary discovery and static-binary auto-fetch (offline)."""

from pathlib import Path

import pytest

from arda import _mmseqs_fetch, mmseqs


@pytest.mark.parametrize(
    "system, machine, expected",
    [
        ("Darwin", "arm64", "mmseqs-osx-universal.tar.gz"),
        ("Darwin", "x86_64", "mmseqs-osx-universal.tar.gz"),
        ("Linux", "x86_64", "mmseqs-linux-avx2.tar.gz"),
        ("Linux", "aarch64", "mmseqs-linux-arm64.tar.gz"),
    ],
)
def test_default_asset_per_platform(monkeypatch, system, machine, expected):
    monkeypatch.delenv("ARDA_MMSEQS_ASSET", raising=False)
    monkeypatch.setattr(_mmseqs_fetch.platform, "system", lambda: system)
    monkeypatch.setattr(_mmseqs_fetch.platform, "machine", lambda: machine)
    assert _mmseqs_fetch.default_asset() == expected


def test_default_asset_env_override(monkeypatch):
    monkeypatch.setenv("ARDA_MMSEQS_ASSET", "mmseqs-linux-sse41.tar.gz")
    assert _mmseqs_fetch.default_asset() == "mmseqs-linux-sse41.tar.gz"


def test_default_asset_unsupported_platform(monkeypatch):
    monkeypatch.delenv("ARDA_MMSEQS_ASSET", raising=False)
    monkeypatch.setattr(_mmseqs_fetch.platform, "system", lambda: "Windows")
    monkeypatch.setattr(_mmseqs_fetch.platform, "machine", lambda: "amd64")
    with pytest.raises(RuntimeError, match="Unsupported platform"):
        _mmseqs_fetch.default_asset()


def test_binary_discovery_prefers_env(monkeypatch):
    monkeypatch.setenv("ARDA_MMSEQS", "/custom/mmseqs")
    mmseqs.mmseqs_binary.cache_clear()
    assert mmseqs.mmseqs_binary() == "/custom/mmseqs"
    mmseqs.mmseqs_binary.cache_clear()


def test_no_auto_fetch_raises_without_network(monkeypatch, tmp_path):
    # Nothing on env / bin / PATH, auto-fetch disabled -> clean error, no download.
    monkeypatch.delenv("ARDA_MMSEQS", raising=False)
    monkeypatch.setenv("ARDA_NO_AUTO_FETCH", "1")
    monkeypatch.setattr(mmseqs, "bin_dir", lambda: tmp_path)
    monkeypatch.setattr(mmseqs.shutil, "which", lambda _: None)

    def _boom(*a, **k):  # auto-fetch must NOT be attempted
        raise AssertionError("auto-fetch attempted despite ARDA_NO_AUTO_FETCH")

    monkeypatch.setattr(mmseqs, "_auto_fetch", _boom)
    mmseqs.mmseqs_binary.cache_clear()
    with pytest.raises(mmseqs.MMseqsError, match="not found"):
        mmseqs.mmseqs_binary()
    mmseqs.mmseqs_binary.cache_clear()


# --- concurrency ---------------------------------------------------------------------------
#
# arda runs concurrently against ONE cache by design (a Nextflow process per sample, a SLURM
# array task per shard), and on a cold cache they all call `fetch` in the same instant. The old
# implementation copied straight onto `bin/mmseqs`, so `dest.exists()` -- the gate saying the
# work was done -- went true on the FIRST BYTE, and a concurrent reader could exec a truncated
# binary. Nothing raised. These pin the properties that make that impossible.


def _fake_archive(tmp: Path, payload: bytes = b"#!/bin/sh\ntrue\n") -> None:
    """Stand in for `_download` + the real tarball: write an archive shaped like MMseqs2's."""
    import tarfile

    src = tmp / "mmseqs" / "bin"
    src.mkdir(parents=True, exist_ok=True)
    (src / "mmseqs").write_bytes(payload)
    with tarfile.open(tmp / "asset.tar.gz", "w:gz") as tf:
        tf.add(src / "mmseqs", arcname="mmseqs/bin/mmseqs")


def _install_fakes(monkeypatch, payload=b"#!/bin/sh\ntrue\n"):
    monkeypatch.setattr(_mmseqs_fetch, "default_asset", lambda: "asset.tar.gz")

    def fake_download(url, dest, **kw):
        _fake_archive(dest.parent, payload)

    monkeypatch.setattr(_mmseqs_fetch, "_download", fake_download)


def test_the_binary_is_published_by_rename_not_written_in_place(tmp_path, monkeypatch):
    """`bin/mmseqs` must appear atomically, already executable.

    Asserted via os.replace rather than by racing threads, because a race test that happens to
    pass proves nothing. The invariant is structural: the only call that creates `dest` is a
    rename, and the source was already chmod'ed.
    """
    import os
    import stat as st

    _install_fakes(monkeypatch)
    seen = {}
    real_replace = os.replace

    def spy_replace(src, dst):
        seen["src"], seen["dst"] = Path(src), Path(dst)
        seen["mode_before_publish"] = Path(src).stat().st_mode
        real_replace(src, dst)

    monkeypatch.setattr(_mmseqs_fetch.os, "replace", spy_replace)
    dest = _mmseqs_fetch.fetch(tmp_path / "bin")

    assert seen["dst"] == dest, "dest must be created by the rename, not written in place"
    assert seen["src"].parent.parent == dest.parent, "staging must share dest's filesystem"
    assert seen["mode_before_publish"] & st.S_IXUSR, "must be executable BEFORE it is published"
    assert dest.stat().st_mode & st.S_IXUSR


def test_a_failed_download_leaves_no_binary_for_anyone_to_exec(tmp_path, monkeypatch):
    """A crash mid-install must not leave a partial `bin/mmseqs` behind.

    This is the failure that had no error attached to it: the next process would find the file,
    treat its existence as "already installed", and run it.
    """
    monkeypatch.setattr(_mmseqs_fetch, "default_asset", lambda: "asset.tar.gz")

    def dying_download(url, dest, **kw):
        dest.write_bytes(b"partial")
        raise RuntimeError("network died mid-download")

    monkeypatch.setattr(_mmseqs_fetch, "_download", dying_download)
    with pytest.raises(RuntimeError, match="network died"):
        _mmseqs_fetch.fetch(tmp_path / "bin")
    assert not (tmp_path / "bin" / "mmseqs").exists()
    # and the lock is released, so the next process is not stuck behind a corpse
    assert not (tmp_path / "bin" / ".mmseqs.lock").exists()


def test_a_second_caller_does_not_download_again(tmp_path, monkeypatch):
    """The lock's `done` check: whoever we queued behind already did the work."""
    _install_fakes(monkeypatch)
    dest = _mmseqs_fetch.fetch(tmp_path / "bin")
    assert dest.exists()

    def must_not_run(*a, **k):
        raise AssertionError("re-downloaded an mmseqs that was already installed")

    monkeypatch.setattr(_mmseqs_fetch, "_download", must_not_run)
    assert _mmseqs_fetch.fetch(tmp_path / "bin") == dest


def test_an_interrupted_copy_never_leaves_a_partial_binary_at_dest(tmp_path, monkeypatch):
    """The actual old bug: copy2 wrote INTO `bin/mmseqs`, so a partial one stayed behind.

    Distinct from the failed-download case, which the old code survived by accident (it died
    before copying). Here the copy itself dies after writing bytes -- which is what a full disk,
    an OOM kill or a SLURM timeout does -- and the old implementation left a truncated,
    executable-looking `bin/mmseqs` that the next process would find, trust and run.
    """
    _install_fakes(monkeypatch)
    real_copy2 = _mmseqs_fetch.shutil.copy2

    def dying_copy2(src, dst, **kw):
        Path(dst).write_bytes(b"#!/bin/sh\ntru")   # a plausible-looking truncation
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(_mmseqs_fetch.shutil, "copy2", dying_copy2)
    with pytest.raises(OSError, match="No space left"):
        _mmseqs_fetch.fetch(tmp_path / "bin")
    assert not (tmp_path / "bin" / "mmseqs").exists(), \
        "a partial binary was published; the next process will exec it"

    # And the cache is not poisoned: a later, healthy run still installs correctly.
    monkeypatch.setattr(_mmseqs_fetch.shutil, "copy2", real_copy2)
    dest = _mmseqs_fetch.fetch(tmp_path / "bin")
    assert dest.exists() and dest.read_bytes() == b"#!/bin/sh\ntrue\n"
