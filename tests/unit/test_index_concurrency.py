"""Several arda processes may build the same mmseqs index at once. That must not corrupt it.

Not hypothetical. arda's Nextflow module launches one process **per sample** and a SLURM array one
per task, so on first use in a fresh environment every one of them finds no index and tries to build
it into the same path.

The failure mode is silent and total. `mmseqs createdb` creates its ``db`` file the instant it starts
writing, so every *other* arda process saw ``db.exists() == True``, skipped building, and searched a
**half-written database**: ``0/200000 reads mapped, loci={}`` -- no error, clean exit code, and a
whole 27-dataset benchmark run of zeros. (``build_index`` was worse: it unlinked the files a
concurrent reader was mid-search on, and mmseqs died with ``Cannot open index file db_h.index.1``.)

The invariant these tests pin: **any process that can see ``db`` must find it complete.**
"""

import multiprocessing as mp
import os
import time
from pathlib import Path

from arda import mmseqs
from arda.annotate import mapper

_PAYLOAD = "X" * 4000


def _slow_write(path: Path) -> None:
    """Stand in for `mmseqs createdb`: the file exists from the first byte, complete much later."""
    with open(path, "w") as fh:
        for i in range(0, len(_PAYLOAD), 100):
            fh.write(_PAYLOAD[i : i + 100])
            fh.flush()
            time.sleep(0.01)


def _build(db: str) -> None:
    mmseqs.createdb = lambda fasta, out, dbtype=2: _slow_write(Path(out))
    mapper._createdb_atomic(Path("/dev/null"), Path(db), 2)


def _read_when_visible(db: str) -> int:
    """A second arda process: waits for the index to appear, then opens it."""
    p = Path(db)
    for _ in range(2000):
        if p.exists():
            time.sleep(0.02)
            return p.stat().st_size
        time.sleep(0.005)
    return -1


def test_a_reader_never_sees_a_half_built_index(tmp_path):
    db = tmp_path / "db"
    writer = mp.Process(target=_build, args=(str(db),))
    pool = mp.Pool(1)
    writer.start()
    size = pool.apply(_read_when_visible, (str(db),))
    writer.join()
    pool.close()
    # Without the atomic move this reads ~300 of 4000 bytes: the database exists, is incomplete, and
    # every subsequent search silently returns no hits.
    assert size == len(_PAYLOAD), (
        f"a concurrent reader saw a {size}-byte index (complete = {len(_PAYLOAD)}): "
        "the db became visible before it was finished"
    )


def _build_tagged(args) -> str:
    db, tag = args
    mmseqs.createdb = lambda fasta, out, dbtype=2: Path(out).write_text(tag * 40)
    mapper._createdb_atomic(Path("/dev/null"), Path(db), 2)
    return tag


def test_concurrent_builders_leave_one_coherent_index(tmp_path):
    db = tmp_path / "db"
    with mp.Pool(4) as pool:
        pool.map(_build_tagged, [(str(db), c) for c in "ABCD"])

    got = db.read_text()
    assert got in {c * 40 for c in "ABCD"}, "the index is a mix of several builders' writes"
    assert not list(tmp_path.glob(".db.lock")), "lock left behind"
    assert not list(tmp_path.glob(".db.tmp.*")), "temp dir left behind"


def test_db_is_moved_into_place_last(tmp_path, monkeypatch):
    """Readers test ``db.exists()``, so ``db`` must land after every sibling it depends on.

    `monkeypatch`, not a bare assignment: patching `mmseqs.createdb` module-globally leaks into every
    later test in the process.
    """
    db = tmp_path / "db"
    moved: list[str] = []

    def createdb(fasta, out, dbtype=2):
        out = Path(out)
        for suffix in (".index", "_h", "_h.index", ".dbtype", ".lookup"):
            (out.parent / (out.name + suffix)).write_text("x")
        out.write_text("db")

    real_replace = mapper.os.replace

    def spy(src, dst):
        moved.append(Path(dst).name)
        return real_replace(src, dst)

    monkeypatch.setattr(mmseqs, "createdb", createdb)
    monkeypatch.setattr(mapper.os, "replace", spy)
    mapper._createdb_atomic(Path("/dev/null"), db, 2)

    assert moved[-1] == "db", f"`db` must be moved into place last; order was {moved}"


def test_a_stale_cache_is_rebuilt(tmp_path, monkeypatch):
    """A db OLDER than its source fasta is stale and must be rebuilt -- "done" means current, not present.

    `_cached_target_db` decides staleness by mtime, then calls `_createdb_atomic` to rebuild. When the
    build guard gated on bare existence, the stale file it found counted as done and the rebuild
    silently no-op'd, so every run kept searching the previous scaffolds. Found as shifted mouse
    markup: a 352-nt TRB scaffold cached under a reference since rebuilt to 346 nt, projecting the
    346-nt coords through a 352-nt alignment and sliding the junction off Cys104.
    """
    fasta = tmp_path / "alleles.fasta"
    fasta.write_text(">x\nACGT\n")
    db = tmp_path / "db"
    db.write_text("STALE")
    old = fasta.stat().st_mtime - 100          # the cache predates its source -> stale
    os.utime(db, (old, old))

    monkeypatch.setattr(mmseqs, "createdb", lambda f, out, dbtype=2: Path(out).write_text("FRESH"))
    mapper._createdb_atomic(fasta, db, 2)
    assert db.read_text() == "FRESH", "a cache older than its source fasta was not rebuilt"
