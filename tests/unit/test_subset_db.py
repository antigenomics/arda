"""`_subset_db` builds the hand-made query sub-DB the two-pass fast path aligns against.

Every defect this file guards has already shipped at least once, and all of them are silent —
a malformed sub-DB does not raise, it makes reads fall through to the full-reference rescue (the
measured rate on amplicon was 16–30 %, which is most of what `--two-pass` gives back in speed) or,
worse, produces alignment rows that cannot be attributed to any read.
"""
from __future__ import annotations

import pytest

from arda import mmseqs
from arda.annotate.mapper import _db_keys, _subset_db

requires_mmseqs = pytest.mark.skipif(
    not mmseqs.available() if hasattr(mmseqs, "available") else False,
    reason="mmseqs not resolvable")


def _make_db(tmp_path, names_seqs):
    fa = tmp_path / "q.fa"
    fa.write_text("".join(f">{n}\n{s}\n" for n, s in names_seqs))
    db = tmp_path / "qdb"
    mmseqs.createdb(fa, db, dbtype=2)
    return db


@requires_mmseqs
def test_createsubdb_writes_its_own_lookup(tmp_path):
    """The file `convertalis` resolves query keys to NAMES through.

    A comment in `_subset_db` long claimed `createsubdb` does not write it. Acting on that and
    writing it by hand overwrote mmseqs' own and broke the rescue path — 490 ids reported "absent
    from lookup". Pin the real behaviour so the same wrong assumption cannot be made twice.
    """
    db = _make_db(tmp_path, [("r1", "ACGTACGTAC"), ("r2", "TTTTGGGGTT"), ("r3", "CCCCAAAACC")])
    sub = _subset_db(db, ["r1", "r3"], tmp_path / "sub")
    assert (tmp_path / "sub.lookup").exists(), "createsubdb no longer writes .lookup"
    sub_keys = _db_keys(sub)
    # It copies the source lookup WHOLE rather than subsetting it -- `r2` is listed even though it
    # is not in the sub-DB. Harmless, because lookup is keyed by the preserved numeric key and a
    # key that is never returned is never asked about; recorded so the surplus entries are not
    # mistaken for a bug later. What matters is that the SELECTED ids resolve.
    assert {"r1", "r3"} <= set(sub_keys)
    assert sub_keys["r1"] == _db_keys(db)["r1"]


@requires_mmseqs
def test_the_sub_db_carries_headers_so_query_names_survive(tmp_path):
    """`createsubdb` acts on ONE database. Without the matching `_h` call, `convertalis` cannot
    print query names and emits numeric keys — which key the result dict by something the caller
    never asked about."""
    db = _make_db(tmp_path, [("readA", "ACGTACGTAC"), ("readB", "TTTTGGGGTT")])
    _subset_db(db, ["readA"], tmp_path / "sub")
    for suffix in ("", ".index", ".dbtype", "_h", "_h.index"):
        assert (tmp_path / f"sub{suffix}").exists(), f"sub{suffix} missing"


@requires_mmseqs
def test_keys_are_preserved_from_the_source(tmp_path):
    """The caller keeps using the SOURCE name->key mapping after subsetting, so `createsubdb` must
    not renumber. If it ever did, every later lookup would silently address the wrong read."""
    db = _make_db(tmp_path, [(f"r{i}", "ACGT" * 5) for i in range(6)])
    src_keys = _db_keys(db)
    sub = _subset_db(db, ["r1", "r4"], tmp_path / "sub")
    sub_keys = _db_keys(sub)
    assert sub_keys["r1"] == src_keys["r1"]
    assert sub_keys["r4"] == src_keys["r4"]


@requires_mmseqs
def test_an_unknown_id_raises_rather_than_being_dropped(tmp_path):
    """This function sits on the no-read-lost path: silently skipping an id it cannot find is
    exactly how the whole fast path once fell through to rescue while still producing correct
    output — fast to miss, because nothing failed."""
    db = _make_db(tmp_path, [("r1", "ACGTACGTAC")])
    with pytest.raises(mmseqs.MMseqsError, match="absent from"):
        _subset_db(db, ["r1", "nope"], tmp_path / "sub")
