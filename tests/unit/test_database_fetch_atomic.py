"""``vdj/`` must never be visible half-populated -- its existence is the gate.

`paths.database_dir` decides "is the reference already here?" by testing whether ``vdj/`` is a
directory. So anything that makes ``vdj/`` appear before it is complete does not cause a slow
download or a crash -- it causes every *other* arda process to search an incomplete reference and
report success over nothing.

Two ways that happened, both fixed here:

* **No lock.** The Nextflow module runs one process per sample, a SLURM array one per task; on first
  use in a fresh cache they all fetch at once, into the same path.
* **A cross-filesystem move.** The old code extracted into ``/tmp`` and called ``shutil.move``, which
  is a rename only *within* one filesystem and silently degrades to a recursive **copy** across them
  -- and ``/tmp`` and ``~/.cache`` normally are two different filesystems. So ``vdj/`` got populated
  file by file, in full view of everyone else.
"""

import multiprocessing as mp
import os
import tarfile
import time
from pathlib import Path

import pytest

from arda import _database_fetch as dbf

# Enough files that a file-by-file copy has a wide window to be caught in the act.
_ORGANISMS = [f"org{i}" for i in range(12)]


def _make_tarball(path: Path) -> None:
    """A stand-in reference archive: vdj/<org>/alleles.fasta, as the real asset is laid out."""
    stage = path.parent / "stage"
    for org in _ORGANISMS:
        d = stage / "vdj" / org
        d.mkdir(parents=True, exist_ok=True)
        (d / "alleles.fasta").write_text(f">{org}\nACGT\n" * 200)
    with tarfile.open(path, "w:gz") as tf:
        tf.add(stage / "vdj", arcname="vdj")


def _complete(vdj: Path) -> bool:
    return all((vdj / o / "alleles.fasta").is_file() for o in _ORGANISMS)


def _fetch(dest: str) -> None:
    """One arda process fetching the reference into a shared cache."""
    tarball = Path(dest).parent / "asset.tar.gz"

    def fake_download(url, out):
        time.sleep(0.05)  # a real download is not instant; give the racers room to interleave
        out.write_bytes(tarball.read_bytes())

    dbf._download = fake_download
    dbf.fetch_database(Path(dest))


def _watch(dest: str) -> str:
    """Another arda process: the moment `vdj/` exists it is considered usable -- is it?"""
    vdj = Path(dest) / "vdj"
    for _ in range(4000):
        if vdj.is_dir():
            # No sleep: the whole point is what a reader sees the instant the gate opens.
            return "complete" if _complete(vdj) else f"PARTIAL: {len(list(vdj.iterdir()))} of {len(_ORGANISMS)}"
        time.sleep(0.002)
    return "never appeared"


def test_a_reader_never_sees_a_partially_extracted_reference(tmp_path):
    """Catches any file-by-file materialization of ``vdj/`` (verified: swap the `os.replace` for a
    `copytree` and this reports ``PARTIAL: 5 of 12``).

    It does **not** by itself reproduce the shipped bug, and that is worth knowing rather than
    assuming: the old `shutil.move` is a plain rename when source and destination share a
    filesystem, which they do under pytest's tmp_path. It only degrades to a copy across a
    filesystem boundary -- i.e. in production, where the staging dir was ``/tmp`` and the cache is
    ``~/.cache``. `test_staging_is_a_sibling_of_the_target` is the deterministic guard for that.
    """
    dest = tmp_path / "database"
    dest.mkdir()
    _make_tarball(tmp_path / "asset.tar.gz")

    watcher = mp.Pool(1)
    async_result = watcher.apply_async(_watch, (str(dest),))
    fetcher = mp.Process(target=_fetch, args=(str(dest),))
    fetcher.start()
    saw = async_result.get(timeout=60)
    fetcher.join()
    watcher.close()

    assert saw == "complete", (
        f"a concurrent arda process saw {saw}. `vdj/` became visible before it was fully "
        "extracted, and its existence is what every other process uses to decide the reference "
        "is ready."
    )


def test_concurrent_fetchers_download_once_and_leave_one_good_tree(tmp_path):
    dest = tmp_path / "database"
    dest.mkdir()
    _make_tarball(tmp_path / "asset.tar.gz")

    with mp.Pool(4) as pool:
        pool.map(_fetch, [str(dest)] * 4)

    assert _complete(dest / "vdj"), "the reference tree is incomplete after 4 concurrent fetchers"
    leftovers = [p.name for p in dest.iterdir() if p.name != "vdj"]
    assert not leftovers, f"staging/lock dirs left behind: {leftovers}"


def test_force_replaces_the_tree_without_ever_exposing_a_partial_one(tmp_path):
    dest = tmp_path / "database"
    dest.mkdir()
    _make_tarball(tmp_path / "asset.tar.gz")
    _fetch(str(dest))

    stale = dest / "vdj" / "org0" / "alleles.fasta"
    stale.write_text(">stale\n")
    dbf.fetch_database(dest, force=True)

    assert _complete(dest / "vdj")
    assert stale.read_text() != ">stale\n", "--force did not actually replace the old tree"
    assert not [p.name for p in dest.iterdir() if p.name != "vdj"], "staging dirs left behind"


@pytest.mark.skipif(os.name == "nt", reason="POSIX rename semantics")
def test_staging_is_a_sibling_of_the_target(tmp_path):
    """The staging dir must live in `dest`, or `os.replace` is a cross-device copy, not a rename.

    This is the actual defect: `/tmp` and `~/.cache` are usually different filesystems, so the old
    `shutil.move` silently copied `vdj/` into place file by file instead of renaming it.
    """
    dest = tmp_path / "database"
    dest.mkdir()
    _make_tarball(tmp_path / "asset.tar.gz")
    seen: list[Path] = []

    real_mkdtemp = dbf.tempfile.mkdtemp

    def spy(*a, **kw):
        d = real_mkdtemp(*a, **kw)
        seen.append(Path(d))
        return d

    dbf.tempfile.mkdtemp = spy
    try:
        _fetch(str(dest))
    finally:
        dbf.tempfile.mkdtemp = real_mkdtemp

    assert seen, "no staging directory was created"
    assert seen[0].parent == dest, (
        f"staged in {seen[0].parent}, not in {dest}: os.replace across filesystems is a copy, "
        "which is exactly the bug"
    )
