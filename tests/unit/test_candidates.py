"""Candidate evidence must survive allele, offset and mismatch ambiguity."""
from types import SimpleNamespace

import pytest

from arda import _markup
from arda.annotate import candidates as mod


def test_d_candidates_keep_alleles_offsets_and_mismatches():
    rows = _markup.d_candidates_batch(["AAACAAA"], ["AAAAAAA", "AAACAAA"])
    assert {row[1] for row in rows} == {0, 1}
    assert len(rows) == 2  # one evidence span per allele, never offset expansion
    assert any(di == 0 and qe - qs + 1 > score for _, di, score, qs, qe, _, _, _ in rows)


@pytest.mark.parametrize("extra", ["", "reverse", "unknown", "qlen", "bounds"])
def test_batched_candidate_api_preserves_all_hits(tmp_path, monkeypatch, extra):
    entries = {"V|a": SimpleNamespace(locus="TRB", v_call="TRBV1*01,TRBV1*02", j_call=""),
               "V|b": SimpleNamespace(locus="TRB", v_call="TRBV2*01", j_call="")}
    monkeypatch.setattr(mod, "load_reference", lambda *a, **k: SimpleNamespace(entries=entries, d_germlines={}))
    monkeypatch.setattr(mod, "_cached_segment_db", lambda *a: tmp_path / "segments")
    monkeypatch.setattr(mod.mmseqs, "createdb", lambda *a, **k: tmp_path / "query")
    calls = []
    monkeypatch.setattr(mod.mmseqs, "search", lambda *a, **k: calls.append(k) or tmp_path / "result")

    def convert(*args, **kwargs):
        args[3].write_text("q1\tV|a\t1\t7\t2\t8\t7\t10\t7\t1\t0\t7M\tAAACAAA\tAAAAAAA\t0.001\t10\t85.7\n"
                           "q1\tV|b\t2\t7\t1\t6\t7\t10\t6\t0\t0\t6M\tAACAAA\tAACAAA\t0.001\t9\t100\n")
    original_convert = convert
    def convert(*args, **kwargs):
        original_convert(*args, **kwargs)
        if extra:
            fields = args[3].read_text().splitlines()[0].split("\t")
            if extra == "reverse": fields[2], fields[3] = fields[3], fields[2]
            if extra == "unknown": fields[0] = "missing"
            if extra == "qlen": fields[6] = "100"
            if extra == "bounds": fields[5] = "100"
            with args[3].open("a") as handle:
                handle.write("\t".join(fields) + "\n")
    monkeypatch.setattr(mod.mmseqs, "convertalis", convert)
    if extra in ("unknown", "qlen", "bounds"):
        with pytest.raises(RuntimeError, match="invalid MMseqs"):
            mod.alignment_candidates([("q1", "aaacaaa"), ("q2", "AAAA")])
        return
    result = mod.alignment_candidates([("q1", "aaacaaa"), ("q2", "AAAA")])
    assert result["limits"]["rejected_reverse_hits"] == int(extra == "reverse")
    assert result["queries"]["sequence"][0] == "AAACAAA"
    assert len(calls) == 1
    hits = result["alignments"].to_dicts()
    assert len(hits) == 2
    assert hits[0]["alleles"] == ["TRBV1*01", "TRBV1*02"]
    assert [(h["qstart"], h["tstart"], h["mismatch"]) for h in hits] == [(1, 2, 1), (2, 1, 0)]
    assert result["queries"].height == 2
    assert {h["sequence_id"] for h in hits} == {"q1"}
    assert result["limits"]["max_seqs"] == 300
    serial = mod.alignment_candidates([("q1", "aaacaaa"), ("q2", "AAAA")], threads=1)
    parallel = mod.alignment_candidates([("q1", "aaacaaa"), ("q2", "AAAA")], threads=4)
    assert serial["alignments"].equals(parallel["alignments"])


@pytest.mark.parametrize("records", [[("q", "AAA"), ("q", "AAA")], [("q x", "AAA")], [("q", "XYZ")]])
def test_candidate_input_validation(records):
    with pytest.raises(ValueError):
        mod.alignment_candidates(records)


def test_empty_schemas_and_bound():
    result = mod.alignment_candidates([])
    assert result["queries"].schema == mod._QUERY_SCHEMA
    assert result["alignments"].schema == mod._SCHEMA
    with pytest.raises(ValueError, match="batch_size"):
        mod.alignment_candidates(((str(i), "AAA") for i in range(100)), batch_size=2)


def test_d_batch_thread_invariant():
    queries = ["AAACAAAGGGTTT" * 5, "TTTCCCAAATGG"] * 32
    germlines = ["AAACAAA", "AAAAGAA", "GGGTTT"]
    assert _markup.d_candidates_batch(queries, germlines, threads=1) == \
        _markup.d_candidates_batch(queries, germlines, threads=4)
    with pytest.raises(ValueError, match="threads"):
        _markup.d_candidates_batch(queries, germlines, threads=0)


@pytest.mark.parametrize("mode", ["lazy", "budget"])
def test_candidate_import_budget_subprocess(mode):
    import os
    import subprocess
    import sys

    program = """
import sys
import arda
assert 'polars' not in sys.modules
""" if mode == "lazy" else """
import polars as pl
assert pl.thread_pool_size() == 4
import arda
try:
    arda.alignment_candidates([], threads=1)
except ValueError as exc:
    assert 'Polars is already initialized' in str(exc)
else:
    raise AssertionError('oversized existing pool accepted')
"""
    result = subprocess.run([sys.executable, "-c", program],
                            env=dict(os.environ, POLARS_MAX_THREADS="4"),
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


def test_candidate_rejects_search_override(monkeypatch):
    monkeypatch.setenv("ARDA_MMSEQS_SEARCH_OPTS", "--strand 2")
    with pytest.raises(ValueError, match="ARDA_MMSEQS_SEARCH_OPTS"):
        mod.alignment_candidates([])
