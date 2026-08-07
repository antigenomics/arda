"""The `v_only` calibration dump: raw scores out, nothing else changed.

`_dump_vonly` exists to answer one question offline -- can a SEGMENT-scale threshold reproduce the
kept set `--min-score 75` produces on whole-scaffold MMseqs2 bit scores? -- for the 43 % of amplicon
mates that hit a V and never reach a J. It is a measurement harness, so the properties that matter
are that it is inert when unset, that it selects exactly the `v_only` population, and that it emits
the scores rather than any decision made from them.
"""

from __future__ import annotations

import pytest

from arda.annotate.mapper import _VONLY_COLS, _dump_vonly


@pytest.fixture
def scene():
    """Four rescued reads, one of each reason. Only `v1` is `v_only`."""
    best_v = {"v1": "IGHV3-23*01", "both": "IGHV1-2*02"}
    best_j = {"j1": "IGHJ4*02", "both": "IGHJ6*03"}
    seg_rows = {
        ("v1", "V"): {"bits": 240.0, "qstart": 1, "qend": 120, "tstart": 176},
        ("both", "V"): {"bits": 210.0, "qstart": 3, "qend": 108, "tstart": 190},
    }
    best = {
        "v1": {"target": "IGH_991", "bits": 131.0},
        "j1": {"target": "IGH_JC_3", "bits": 88.0},
        "both": {"target": "IGH_412", "bits": 205.0},
    }
    seqs = {"v1": "A" * 150, "j1": "C" * 150, "both": "G" * 150, "c1": "T" * 150}
    return ["both", "c1", "j1", "v1"], best_v, best_j, seg_rows, best, seqs


def _rows(path):
    lines = path.read_text().splitlines()
    return [dict(zip(_VONLY_COLS, line.split("\t"))) for line in lines[1:]]


def test_writes_nothing_at_all_when_the_env_var_is_unset(scene, tmp_path, monkeypatch):
    # The whole point of an env-var harness is that a normal run is untouched. If this ever writes
    # by default it becomes a per-chunk file append on every production run.
    monkeypatch.delenv("ARDA_VONLY_DUMP", raising=False)
    _dump_vonly(*scene)
    assert list(tmp_path.iterdir()) == []


def test_selects_v_only_and_nothing_else(scene, tmp_path, monkeypatch):
    # `j_only`, `c_only` and reads with BOTH sides also pass through `rescue`; a dump that caught
    # them would inflate the population whose cost this change claims to remove (93.4 %, not 100 %).
    out = tmp_path / "vonly.tsv"
    monkeypatch.setenv("ARDA_VONLY_DUMP", str(out))
    _dump_vonly(*scene)
    assert [r["read_id"] for r in _rows(out)] == ["v1"]


def test_emits_both_scores_on_their_own_scales(scene, tmp_path, monkeypatch):
    # The two numbers are the experiment. `seg_bits` is ungapped MATCH 2 / MISMATCH -3; the
    # scaffold bits are MMseqs2's. Emitting only one, or a ratio, would presuppose the answer.
    out = tmp_path / "vonly.tsv"
    monkeypatch.setenv("ARDA_VONLY_DUMP", str(out))
    _dump_vonly(*scene)
    (row,) = _rows(out)
    assert float(row["seg_bits"]) == 240.0
    assert float(row["scaffold_bits"]) == 131.0
    assert row["scaffold"] == "IGH_991"
    assert row["v_allele"] == "IGHV3-23*01"


def test_emits_the_length_terms_a_per_locus_fit_needs(scene, tmp_path, monkeypatch):
    # An ungapped score grows with aligned length, and V genes differ in length across loci, so a
    # single global cut on raw `seg_bits` is not defensible. The fit needs to normalise, which it
    # can only do if the coordinates come out with the score.
    out = tmp_path / "vonly.tsv"
    monkeypatch.setenv("ARDA_VONLY_DUMP", str(out))
    _dump_vonly(*scene)
    (row,) = _rows(out)
    assert (int(row["seg_qstart"]), int(row["seg_qend"])) == (1, 120)
    assert int(row["read_len"]) == 150


def test_appends_across_chunks_with_exactly_one_header(scene, tmp_path, monkeypatch):
    # `map` calls this once per chunk (200k reads), so a second call must extend the same table.
    # A header per chunk would break every reader silently -- the extra lines parse as data.
    out = tmp_path / "vonly.tsv"
    monkeypatch.setenv("ARDA_VONLY_DUMP", str(out))
    _dump_vonly(*scene)
    _dump_vonly(*scene)
    text = out.read_text()
    assert text.count("read_id") == 1
    assert len(_rows(out)) == 2


def test_a_v_only_read_the_rescue_could_not_place_still_appears(tmp_path, monkeypatch):
    # A read the full reference returned nothing for is the most interesting row in the table: it
    # is a v_only read with a segment score and NO scaffold score, i.e. one the current gate drops.
    # Skipping it would bias the fit towards reads that already pass.
    out = tmp_path / "vonly.tsv"
    monkeypatch.setenv("ARDA_VONLY_DUMP", str(out))
    _dump_vonly(["v1"], {"v1": "IGHV3-23*01"}, {},
                {("v1", "V"): {"bits": 61.0, "qstart": 1, "qend": 40, "tstart": 250}},
                {}, {"v1": "A" * 150})
    (row,) = _rows(out)
    assert float(row["seg_bits"]) == 61.0
    assert row["scaffold_bits"] == "" and row["scaffold"] == ""
