"""Unit tests for the library adapter's input normalization.

``annotate_sequences`` accepts either raw strings or ``(id, sequence)`` pairs and
must hand ``annotate_records`` a list of ``(id, sequence)`` tuples. We stub the
heavy mapper to assert only the normalization contract.
"""

import arda.annotate.mapper as mapper_mod
from arda.adapter import annotate_sequences


def test_raw_strings_get_synthetic_ids(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        mapper_mod, "annotate_records",
        lambda pairs, **kw: captured.update(pairs=pairs, kw=kw) or [])
    annotate_sequences(["ACGT", "TTTT"], organism="mouse", seqtype="nt", map_d=False)
    assert captured["pairs"] == [("seq0", "ACGT"), ("seq1", "TTTT")]
    assert captured["kw"] == {"organism": "mouse", "seqtype": "nt", "map_d": False}


def test_id_sequence_pairs_are_preserved(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        mapper_mod, "annotate_records",
        lambda pairs, **kw: captured.update(pairs=pairs) or [])
    annotate_sequences([("read1", "ACGT"), ("read2", "TTTT")])
    assert captured["pairs"] == [("read1", "ACGT"), ("read2", "TTTT")]
