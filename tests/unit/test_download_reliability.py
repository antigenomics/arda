"""Small reference/tool downloads must fail promptly and reject incomplete HTTP payloads."""
import io
import urllib.error

import pytest

from arda import _mmseqs_fetch as fetch


def test_permanent_http_error_is_not_retried(tmp_path, monkeypatch):
    calls = []

    def open_url(req, **kw):
        calls.append(kw)
        raise urllib.error.HTTPError(req.full_url, 404, "missing", {}, None)

    monkeypatch.setattr(fetch.urllib.request, "urlopen", open_url)
    monkeypatch.setattr(fetch.time, "sleep", lambda _: None)
    with pytest.raises(RuntimeError, match="404"):
        fetch._download("https://example.test/missing", tmp_path / "download")
    assert len(calls) == 1
    assert calls[0]["timeout"] <= 30


def test_truncated_http_body_never_counts_as_success(tmp_path, monkeypatch):
    class Response(io.BytesIO):
        headers = {"Content-Length": "4"}

    calls = []

    def open_url(req, **kw):
        calls.append(kw)
        return Response(b"ab")

    monkeypatch.setattr(fetch.urllib.request, "urlopen", open_url)
    monkeypatch.setattr(fetch.time, "sleep", lambda _: None)
    with pytest.raises(RuntimeError, match="[Ii]ncomplete|[Tt]runcated"):
        fetch._download("https://example.test/archive", tmp_path / "download")
    assert len(calls) == 3
