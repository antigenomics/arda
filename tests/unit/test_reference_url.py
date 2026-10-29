"""The reference asset URL must not be derived from the package version.

The reference is data; the package version changes every release. Deriving one from the other
coupled them silently, so bumping the version ahead of a GitHub release made every cold-cache
install die on `HTTP Error 404` for an asset that was never going to exist at that tag.

Found by running eight concurrent arda processes against a fresh cache on a cluster — the same
experiment that was supposed to be testing the build lock. It would have hit the first real
`pip install` of any release whose asset was not published yet.

No network: these only inspect URL construction.
"""

from __future__ import annotations

from arda import __version__
from arda._database_fetch import _REFERENCE_TAG, _candidate_urls, reference_url


def test_default_url_uses_the_reference_tag_not_the_package_version():
    url = reference_url()
    assert f"/v{_REFERENCE_TAG}/" in url
    if _REFERENCE_TAG != __version__:
        assert f"/v{__version__}/" not in url, (
            "the reference URL is being derived from the package version again; bumping the "
            "version would 404 every cold-cache install until a release asset exists"
        )


def test_candidates_fall_back_to_the_package_version():
    """A release that DID ship its own asset before the tag was updated still resolves."""
    urls = _candidate_urls()
    assert urls[0] == reference_url()
    if _REFERENCE_TAG != __version__:
        assert any(f"/v{__version__}/" in u for u in urls)


def test_an_explicit_version_is_the_only_candidate():
    """A caller naming a version must not silently be given a different one."""
    urls = _candidate_urls("2.5.6")
    assert urls == [reference_url("2.5.6")]
    assert len(urls) == 1
    assert "/v2.5.6/" in urls[0]


def test_env_override(monkeypatch):
    monkeypatch.setenv("ARDA_REFERENCE_TAG", "v9.9.9")
    assert "/v9.9.9/" in reference_url()


def test_leading_v_is_not_doubled(monkeypatch):
    monkeypatch.setenv("ARDA_REFERENCE_TAG", "v1.2.3")
    assert "/vv1.2.3/" not in reference_url()
    assert "/v1.2.3/" in reference_url()


def test_explicit_environment_tag_is_not_silently_replaced(monkeypatch):
    monkeypatch.setenv("ARDA_REFERENCE_TAG", "v9.9.9")
    assert _candidate_urls() == [reference_url()]


def test_pinned_data_manifest_covers_current_completeness_resources():
    import hashlib
    import json
    from pathlib import Path
    from arda import _database_fetch as dbf
    manifest=json.loads(Path(dbf.__file__).with_name('_reference_manifest.json').read_text())
    assert manifest['tag'] == dbf._REFERENCE_TAG
    base=Path(__file__).resolve().parents[2]/'database'
    for relative,digest in manifest['files_sha256'].items():
        assert hashlib.sha256((base/relative).read_bytes()).hexdigest() == digest
    for organism in dbf._ORGANISMS:
        assert f'vdj/{organism}/d_germlines.tsv' in manifest['files_sha256']
        assert 'd_germlines.tsv' in dbf._REQUIRED
