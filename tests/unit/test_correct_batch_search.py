"""Correction submits whole query batches and preserves filtered reference indices."""
import pytest

from arda.rnaseq import correct


def test_parent_and_rate_search_use_batches(monkeypatch):
    original = correct.seqtree.Index
    calls = []

    class BatchedIndex:
        @staticmethod
        def build(*args, **kwargs):
            native = original.build(*args, **kwargs)

            class Search:
                def search(self, *args, **kwargs):
                    raise AssertionError('per-junction search')

                def search_batch(self, queries, params, threads):
                    calls.append((list(queries), threads))
                    return native.search_batch(queries, params, threads)

            return Search()

    monkeypatch.setattr(correct.seqtree, 'Index', BatchedIndex)
    junctions = ['AAA', 'AAC', 'NNN', '']
    counts = [10000, 10, 1, 1]
    v, j = ['V'] * 4, ['J'] * 4
    parents = correct._parents(junctions, counts, v, j, max_subs=1, max_indel=0,
        error_rate=.01, indel_rate=.001, require_vj=True, threads=1)
    assert parents == [None, 0, None, None]
    rate, _ = correct.estimate_error_rate(junctions, counts, v, j, threads=1)
    assert rate == pytest.approx(10 / 10000 / 3)
    assert calls == [(['AAA', 'AAC'], 1), (['AAA'], 1)]
