"""Cold-cache races, incomplete downloads and refreshes must never expose partial references."""
import multiprocessing as mp
import shutil
import tarfile
import time
from pathlib import Path

import pytest

from arda import _database_fetch as dbf

_ORGANISMS = ("human", "mouse", "rat", "rabbit", "rhesus_monkey")
@pytest.fixture(autouse=True)
def synthetic_reference_identity(monkeypatch):
    # These synthetic archives are not the published dataset. Identity tests override this map.
    monkeypatch.setattr(dbf, "_PINNED_FILES", {})


_FILES = ("alleles.fasta", "alleles.aa.fasta", "markup.tsv", "markup.aa.tsv",
          "cdr3_anchors.tsv", "combinations.tsv", "d_germlines.fasta")


def _make_tarball(path, label="original"):
    stage = path.parent / ("stage-" + label)
    for org in _ORGANISMS:
        d = stage / "vdj" / org
        d.mkdir(parents=True, exist_ok=True)
        for name in _FILES + (("d_prior.tsv",) if org in ("human", "mouse") else ()):
            (d / name).write_text((label + "\n") * 200)
    with tarfile.open(path, "w:gz") as tf:
        tf.add(stage / "vdj", arcname="vdj")


def _fetch(dest):
    dbf._PINNED_FILES = {}  # spawned workers use the synthetic archive, not real pinned data
    asset = Path(dest).parent / "asset.tar.gz"

    def download(url, out):
        with (Path(dest).parent / "downloads.log").open("a") as log:
            log.write("download\n")
        time.sleep(0.05)
        shutil.copyfile(asset, out)

    dbf._download = download
    return str(dbf.fetch_database(Path(dest)))


def test_concurrent_fetchers_download_once_and_share_one_complete_generation(tmp_path):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    with mp.get_context("spawn").Pool(4) as pool:
        roots = pool.map(_fetch, [str(dest)] * 4)
    assert len(set(roots)) == 1
    assert (tmp_path / "downloads.log").read_text().splitlines() == ["download"]
    assert len(list((dest / ".references").iterdir())) == 1
    root = Path(roots[0])
    assert all((root / "vdj" / o / f).is_file() for o in _ORGANISMS for f in _FILES)
    assert dbf.installed_database(dest) == root
    assert not list(dest.glob(".arda_db_*"))


def test_empty_legacy_cache_is_repaired(tmp_path, monkeypatch):
    dest = tmp_path / "database"
    (dest / "vdj").mkdir(parents=True)
    _make_tarball(tmp_path / "asset.tar.gz")
    calls = []

    def download(url, out):
        calls.append(url)
        shutil.copyfile(tmp_path / "asset.tar.gz", out)

    monkeypatch.setattr(dbf, "_download", download)
    root = dbf.fetch_database(dest)
    assert calls
    assert (root / "vdj/human/markup.tsv").stat().st_size > 1000


def test_incomplete_archive_is_rejected_before_publication(tmp_path, monkeypatch):
    (tmp_path / "empty/vdj").mkdir(parents=True)
    asset = tmp_path / "asset.tar.gz"
    with tarfile.open(asset, "w:gz") as tf:
        tf.add(tmp_path / "empty/vdj", arcname="vdj")
    monkeypatch.setattr(dbf, "_download", lambda url, out: shutil.copyfile(asset, out))
    with pytest.raises(RuntimeError, match="reference|Reference"):
        dbf.fetch_database(tmp_path / "database")
    assert not (tmp_path / "database/REFERENCE.json").exists()


def test_failed_refresh_preserves_current_generation(tmp_path, monkeypatch):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    old = Path(_fetch(str(dest)))
    content = (old / "vdj/human/alleles.fasta").read_bytes()
    replace = dbf.os.replace

    def fail_pointer(src, target):
        if Path(target) == dest / "REFERENCE.json" or Path(target) == dest / "vdj":
            raise OSError("simulated publication failure")
        return replace(src, target)

    monkeypatch.setattr(dbf.os, "replace", fail_pointer)
    with pytest.raises(OSError, match="publication failure"):
        dbf.fetch_database(dest, force=True)
    assert (old / "vdj/human/alleles.fasta").read_bytes() == content
    assert dbf.installed_database(dest) == old


def test_refresh_keeps_reader_paths_immutable(tmp_path):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    old = Path(_fetch(str(dest)))
    content = (old / "vdj/human/alleles.fasta").read_bytes()
    _make_tarball(tmp_path / "asset.tar.gz", "replacement")
    new = dbf.fetch_database(dest, force=True)
    assert old != new
    assert (old / "vdj/human/alleles.fasta").read_bytes() == content
    assert (new / "vdj/human/alleles.fasta").read_bytes() != content
    assert dbf.installed_database(dest) == new


def test_corrupted_cached_file_is_repaired(tmp_path):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    old = Path(_fetch(str(dest)))
    (old / "vdj/human/markup.tsv").write_text("broken\n")
    new = dbf.fetch_database(dest)
    assert old != new
    assert (new / "vdj/human/markup.tsv").stat().st_size > 1000


def test_staging_is_on_the_destination_filesystem(tmp_path, monkeypatch):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    seen = []
    mkdtemp = dbf.tempfile.mkdtemp

    def spy(*a, **kw):
        path = mkdtemp(*a, **kw)
        seen.append(Path(path))
        return path

    monkeypatch.setattr(dbf.tempfile, "mkdtemp", spy)
    _fetch(str(dest))
    assert seen and all(p.parent == dest for p in seen)


def test_known_legacy_reference_is_adopted_without_network(tmp_path, monkeypatch):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    shutil.copytree(tmp_path / "stage-original/vdj", dest / "vdj")
    expected = {p.relative_to(dest).as_posix(): dbf._sha256(p)
                for p in (dest / "vdj").rglob("*") if p.is_file()}
    monkeypatch.setattr(dbf, "_PINNED_FILES", expected, raising=False)
    monkeypatch.setattr(dbf, "_download", lambda *a: pytest.fail("complete cache must stay offline"))
    root = dbf.fetch_database(dest)
    assert root != dest
    assert dbf.installed_database(dest) == root
    assert (dest / "vdj/human/alleles.fasta").is_file(), "old readers retain their path"


def test_explicit_version_cannot_reuse_a_different_generation(tmp_path, monkeypatch):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    old = Path(_fetch(str(dest)))
    calls = []

    def download(url, out):
        calls.append(url)
        shutil.copyfile(tmp_path / "asset.tar.gz", out)

    monkeypatch.setattr(dbf, "_download", download)
    new = dbf.fetch_database(dest, version="different")
    assert old != new
    assert calls == [dbf.reference_url("different")]
    assert dbf.installed_database(dest, version="different") == new
    assert dbf.installed_database(dest) is None


def test_pinned_asset_payload_must_match_known_file_hashes(tmp_path, monkeypatch):
    _make_tarball(tmp_path / "asset.tar.gz")
    monkeypatch.setattr(dbf, "_PINNED_FILES", {"vdj/human/alleles.fasta": "0" * 64})
    monkeypatch.setattr(dbf, "_download", lambda url, out: shutil.copyfile(tmp_path / "asset.tar.gz", out))
    with pytest.raises(RuntimeError, match="[Cc]hecksum"):
        dbf.fetch_database(tmp_path / "database")


def test_legacy_copy_is_revalidated_before_publication(tmp_path, monkeypatch):
    dest = tmp_path / "database"
    _make_tarball(tmp_path / "asset.tar.gz")
    shutil.copytree(tmp_path / "stage-original/vdj", dest / "vdj")
    expected = {p.relative_to(dest).as_posix(): dbf._sha256(p)
                for p in (dest / "vdj").rglob("*") if p.is_file()}
    monkeypatch.setattr(dbf, "_PINNED_FILES", expected)
    copy = dbf.shutil.copy2

    def changed_during_copy(src, target):
        result = copy(src, target)
        if str(target).endswith("human/markup.tsv"):
            Path(target).write_text("changed during staging\n")
        return result

    monkeypatch.setattr(dbf.shutil, "copy2", changed_during_copy)
    with pytest.raises(RuntimeError, match="[Cc]hecksum"):
        dbf.fetch_database(dest)
    assert not (dest / "REFERENCE.json").exists()
