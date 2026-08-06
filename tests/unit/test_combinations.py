"""Unit tests for scaffold enumeration and frame-preserving padding."""

from arda.refbuild.loci import LOCI
from arda.refbuild import combinations
from arda.refbuild.translate import translate


def _locus(name):
    return next(l for l in LOCI if l.name == name)


def test_vj_padding_preserves_frame():
    # Synthetic V (frame 0, len 30) + J with coding frame 2.
    v = {"V1": "ATG" * 10}            # len 30
    j = {"J1": "CC" + "TGGGGGCAGGGG"}  # J coding starts at offset 2 -> WGQG...
    frames = {"J1": 2}
    sc = combinations.build_locus_scaffolds(_locus("TRA"), v, j, frames)
    assert len(sc) == 1
    s = sc[0]
    # (len_V + n_pad + jframe) must be a multiple of 3
    assert (30 + s.n_pad + 2) % 3 == 0
    # whole scaffold reads in frame 0 with the J coding in frame
    prot = translate(s.sequence, 0)
    assert "W" in prot  # the conserved J tryptophan appears in frame


def test_vdj_gets_d_spacer_vj_does_not():
    v = {"V1": "ATG" * 10}
    j = {"J1": "TGGGGGCAGGGG"}
    frames = {"J1": 0}
    vj = combinations.build_locus_scaffolds(_locus("TRA"), v, j, frames)[0]
    vdj = combinations.build_locus_scaffolds(_locus("IGH"), v, j, frames)[0]
    assert vdj.n_pad - vj.n_pad == combinations.DEFAULT_D_SPACER_NT


def test_trd_configured_to_share_trav_dv_v_genes():
    # Regression guard (DB-free, so it runs in CI unlike test_locus_disambiguation which needs a DB):
    # TRA and TRD share V genes filed under TRAV as ".../DV..."; without them a δ rearrangement on a
    # shared V gene has no TRD scaffold and is miscalled TRA. The locus follows J, so a TRAV/DV + TRDJ
    # scaffold must exist and be labelled TRD.
    trd = _locus("TRD")
    assert trd.v_shared == ("TRAV", "/DV")
    v = {"TRAV14/DV4*01": "ATG" * 10}
    j = {"TRDJ1*01": "TGGGGGCAGGGG"}
    sc = combinations.build_locus_scaffolds(trd, v, j, {"TRDJ1*01": 0})
    assert len(sc) == 1 and sc[0].locus == "TRD"
    assert sc[0].v_calls == ["TRAV14/DV4*01"]


def test_tra_configured_to_share_trdv_v_genes():
    """The INVERSE of the rule above, which was missing and cost real junctions.

    TRDV1/2/3 are dedicated delta V genes but lie *inside* the TRA locus on chr14, between the TRAV
    genes and the TRAJ cluster, so an alpha rearrangement can join one to a TRAJ. Without a
    TRDV x TRAJ scaffold such a read has nowhere to land: measured on the PRJNA371303 TRA amplicon,
    IgBLAST calls **147 of 9,300** truth reads `TRDV1*01` + a TRAJ with a real junction, and arda
    calls every one of those J genes correctly while emitting `v_call = null` and therefore no
    junction -- `junction_aa` accuracy **0.0952** on that stratum against 0.9049 on TRA overall.

    DB-free, so it runs in CI. The locus follows the J, so the scaffold must be labelled TRA.
    """
    tra = _locus("TRA")
    assert tra.v_shared == ("TRDV", ""), "TRA must pull the whole TRDV stem"
    v = {"TRDV1*01": "ATG" * 10}
    j = {"TRAJ7*01": "TGGGGGCAGGGG"}
    sc = combinations.build_locus_scaffolds(tra, v, j, {"TRAJ7*01": 0})
    assert len(sc) == 1 and sc[0].locus == "TRA"
    assert sc[0].v_calls == ["TRDV1*01"]


def test_dedup_collapses_identical_scaffolds():
    # Two V alleles with identical sequence collapse to one scaffold.
    v = {"V1": "ATG" * 10, "V2": "ATG" * 10}
    j = {"J1": "TGGGGGCAGGGG"}
    sc = combinations.build_locus_scaffolds(_locus("TRA"), v, j, {"J1": 0})
    assert len(sc) == 1
    assert sc[0].v_calls == ["V1", "V2"]


def test_load_j_frames_parses_aux_skipping_comments(tmp_path, monkeypatch):
    aux_dir = tmp_path / "optional_file"
    aux_dir.mkdir()
    (aux_dir / "human_gl.aux").write_text(
        "# comment line\n"
        "IGHJ1*01\t1\t13\n"
        "IGHJ2*01\t2\t10\n"
        "malformed_line_without_frame\n"
        "\n"
    )
    monkeypatch.setattr(combinations, "bin_dir", lambda: tmp_path)
    frames = combinations.load_j_frames("human")
    assert frames == {"IGHJ1*01": 1, "IGHJ2*01": 2}


def test_load_j_frames_missing_file_returns_empty(tmp_path, monkeypatch):
    monkeypatch.setattr(combinations, "bin_dir", lambda: tmp_path)
    assert combinations.load_j_frames("nonexistent") == {}
