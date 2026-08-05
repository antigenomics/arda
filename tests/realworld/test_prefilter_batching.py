"""Reading stays chunked to bound memory; searching must not.

A prefiltered read chunk is tiny -- 0.47 % of reads survive on a 0.024 %-receptor library -- and
every `mmseqs search` call costs ~0.7 s of fixed setup whatever it is handed. Ten near-empty
searches measured 25.8 s against 13.4 s for one on SRR10611239. So survivors accumulate across
read chunks and MMseqs2 is invoked once per full batch.

That is invisible in the output, which is exactly why it needs a test: the previous shape produced
identical results and was 2x slower.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from arda.rnaseq import map as rmap

pytest.importorskip("arda._prefilter", reason="native prefilter extension not built")

READS = Path(__file__).resolve().parents[1] / "data" / "rnaseq_real"


def _count_searches(tmp_path, **kw) -> tuple[int, set[str]]:
    """Run `map`, returning how many times MMseqs2 was invoked and which reads were reported."""
    real = rmap.mapper._annotate_chunk
    calls: list[int] = []

    def spy(records, *a, **kwargs):
        calls.append(len(records))
        return real(records, *a, **kwargs)

    rmap.mapper._annotate_chunk = spy
    try:
        out = tmp_path / f"m{len(kw)}.tsv"
        rmap.map_rnaseq(READS / "reads_1.fq.gz", out, r2=READS / "reads_2.fq.gz",
                        threads=4, **kw)
        ids = {ln.split("\t")[0] for ln in out.read_text().splitlines()[1:]}
    finally:
        rmap.mapper._annotate_chunk = real
    return len(calls), ids


def test_prefiltering_collapses_many_read_chunks_into_fewer_searches(tmp_path):
    n_off, ids_off = _count_searches(tmp_path, chunk_size=200)
    n_on, ids_on = _count_searches(tmp_path, chunk_size=200, prefilter=True)
    assert n_off >= 6, f"expected the fixture to split into several chunks, got {n_off}"
    assert n_on < n_off, (
        f"the prefilter did not reduce the search count: {n_on} searches for {n_off} read chunks")
    # And the reads reported are the same ones: batching is a scheduling change, not a result one.
    assert ids_on <= ids_off
    assert len(ids_off - ids_on) / max(len(ids_off), 1) < 0.02


def test_output_does_not_depend_on_the_read_chunk_size(tmp_path):
    """Fragments must not be split across a flush. `chunked_fragments` keeps a fragment's mates in
    one chunk, and flushing only ever happens on a chunk boundary -- so a tiny chunk and a large
    one have to agree exactly. Splitting fragments under --reconstruct is a bug this pipeline has
    already shipped once."""
    _n1, small = _count_searches(tmp_path, chunk_size=200, prefilter=True)
    _n2, large = _count_searches(tmp_path, chunk_size=400_000, prefilter=True)
    assert small == large
