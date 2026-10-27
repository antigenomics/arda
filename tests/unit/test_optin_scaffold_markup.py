"""Recover only exact, anchor-backed J markup in opt-in germline scaffolds."""

from arda.refbuild.build import _complete_scaffold_j
from arda.refbuild.combinations import Scaffold


def test_exact_j_geometry_including_noncanonical_anchor():
    # Cys118 is legitimate: the aux coordinate, not an F/W motif, defines FR4.
    v, j = "ATGTGTGCT", "AAATGCGGTGCTGGT"
    sc = Scaffold("TRA_0", "TRA", v + j, ["TRDV1*01"], ["TRAJ35*01"])
    rec = {"v_call": "TRDV1*01", "fwr1_start": "1", "fwr3_end": "6"}
    assert _complete_scaffold_j(rec, sc, v, j, (2, 0))
    assert rec["cdr3_start"] == 7 and rec["cdr3_end"] == 12
    assert rec["fwr4_start"] == 13 and rec["fwr4_end"] == 24
    assert rec["junction"] == "TGTGCTAAATGC"
    assert rec["junction_aa"] == "CAKC"
    for region in ("cdr3", "fwr4"):
        assert rec[region] == sc.sequence[rec[f"{region}_start"] - 1:rec[f"{region}_end"]]


def test_missing_aux_or_mismatched_germline_cannot_invent_markup():
    v, j = "ATGTGTGCT", "AAATTTGGTGCTGGT"
    sc = Scaffold("TRA_0", "TRA", v + j, ["TRDV1*01"], ["TRAJ1*01"])
    rec = {"v_call": "TRDV1*01", "fwr1_start": "1", "fwr3_end": "6"}
    assert not _complete_scaffold_j(rec.copy(), sc, v, j, None)
    assert not _complete_scaffold_j(rec.copy(), sc, v, "CCC" + j, (2, 0))
    assert not _complete_scaffold_j(dict(rec, fwr3_end="3"), sc, v, j, (2, 0))


def test_partial_v_is_normalized_exactly_as_in_scaffold_construction():
    from arda.refbuild.combinations import build_locus_scaffolds
    from arda.refbuild.loci import loci_for

    v = "A" + "GCTAAA" + "GCT" * 7 + "TGT"
    j = "AAATTTGGTGCTGGT"
    sc, = build_locus_scaffolds(loci_for(allow_chimeras=True)[0],
                                {"TRDV1*01": v}, {"TRAJ1*01": j}, {"TRAJ1*01": 0})
    rec = {"fwr1_start": "1", "fwr3_end": str(len(v) - 1)}
    assert _complete_scaffold_j(rec, sc, v, j, (2, 0))
    assert rec["j_sequence_start"] == len(v)


def test_optin_build_uses_anchor_frame_and_completes_every_j(tmp_path, monkeypatch):
    import importlib
    import logging

    import polars as pl

    from arda.refbuild.loci import loci_for

    build = importlib.import_module("arda.refbuild.build")
    v = "GCT" * 9 + "TGT"
    js = {"TRAJ31*01": "AAATTTGGTGCTGGT", "TRAJ41*02": "AAATTTGGTGGTGGT"}
    monkeypatch.setattr(build.imgt, "load_functional_alleles", lambda sd, group, stem:
                        {"TRDV1*01": v} if stem == "TRDV" else js if stem == "TRAJ" else {})
    monkeypatch.setattr(build.combinations, "load_v_fwr3_stops", lambda _: {"TRDV1*01": 30})
    monkeypatch.setattr(build.combinations, "load_j_fr4_offsets", lambda _: {"TRAJ31*01": (2, 0)})
    monkeypatch.setattr(build, "_scaffold_fasta_path", lambda *args: tmp_path / "sc.fasta")
    monkeypatch.setattr("arda.igblast.has_internal_annotation", lambda *args: True)

    def annotate(path, *args, **kwargs):
        return pl.DataFrame([
            dict(sequence_id=sid, v_call="TRDV1*01", fwr1_start="1", fwr1_end="3",
                 cdr1_start="4", cdr1_end="6", fwr2_start="7", fwr2_end="9",
                 cdr2_start="10", cdr2_end="12", fwr3_start="13", fwr3_end="30",
                 **{region: v[start:end] for region, start, end in
                    [("fwr1", 0, 3), ("cdr1", 3, 6), ("fwr2", 6, 9),
                     ("cdr2", 9, 12), ("fwr3", 12, 30)]})
            for sid, seq in build.imgt.read_fasta(path)
        ])

    monkeypatch.setattr(build.airr_extract, "annotate_scaffolds", annotate)
    rows, aa, combos, fasta, proteins = build._process_locus(
        "human", "test", loci_for(allow_chimeras=True)[0], {"TRAJ31*01": 1}, logging.getLogger())
    assert len(rows) == 2
    assert {r["j_call"] for r in rows} == set(js)
    for row in rows:
        sequence = dict(fasta)[row["scaffold_id"]]
        for region in build.REGION_NAMES:
            assert row[region] == sequence[row[f"{region}_start"] - 1:row[f"{region}_end"]]
        assert row["junction"] == sequence[row["fwr3_end"] - 3:row["fwr4_start"] + 2]


def test_optin_appends_scaffolds_without_renumbering_or_reframing_existing_targets(tmp_path, monkeypatch):
    """New vocabulary must not reshuffle equal-score legacy alignments or their anchors."""
    import importlib
    import logging
    import pytest
    from arda.refbuild.loci import loci_for
    from arda.refbuild.combinations import build_locus_scaffolds

    build = importlib.import_module("arda.refbuild.build")
    vs = {"TRAV1*01": "GCT" * 9 + "TGT"}
    ds = {"TRDV1*01": "GCC" * 9 + "TGT"}
    js = {"TRAJ31*01": "AAATTTGGTGCTGGT"}
    frames = {"TRAJ31*01": 1}  # historical aux frame differs from the trusted anchor frame
    expected = build_locus_scaffolds(loci_for()[0], vs, js, frames)
    monkeypatch.setattr(build.imgt, "load_functional_alleles", lambda sd, group, stem:
                        ds if stem == "TRDV" else js if stem == "TRAJ" else vs.copy())
    monkeypatch.setattr(build.combinations, "load_v_fwr3_stops", lambda _: {a: 30 for a in vs | ds})
    monkeypatch.setattr(build.combinations, "load_j_fr4_offsets", lambda _: {"TRAJ31*01": (2, 0)})
    monkeypatch.setattr(build, "_scaffold_fasta_path", lambda *args: tmp_path / "sc.fasta")
    monkeypatch.setattr("arda.igblast.has_internal_annotation", lambda *args: True)
    captured = []

    class ReachedAnnotation(BaseException):
        pass

    def annotate(path, *args, **kw):
        captured.extend(build.imgt.read_fasta(path))
        raise ReachedAnnotation

    monkeypatch.setattr(build.airr_extract, "annotate_scaffolds", annotate)
    with pytest.raises(ReachedAnnotation):
        build._process_locus("human", "test", loci_for(allow_chimeras=True)[0], frames,
                             logging.getLogger())
    assert captured[:len(expected)] == [(s.scaffold_id, s.sequence) for s in expected]
    assert len(captured) == len(expected) + 1
    assert captured[-1][0] == f"TRA_{len(expected)}"


def test_rebuilding_requires_an_isolated_source_checkout(tmp_path, monkeypatch):
    import importlib
    import pytest
    build = importlib.import_module("arda.refbuild.build")
    monkeypatch.setattr(build, "_source_root", lambda: None, raising=False)
    monkeypatch.setattr(build, "vdj_dir", lambda _: tmp_path / "output")
    monkeypatch.setattr(build.imgt, "download_reference", lambda: pytest.fail("must fail before downloading"))
    with pytest.raises(RuntimeError, match="source checkout"):
        build.build_species("human")
