"""Unit tests for mmseqs binary discovery and static-binary auto-fetch (offline)."""

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
