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


def test_tra_does_NOT_share_the_trdv_stem():
    """The V-gene sharing between TRA and TRD is ASYMMETRIC, and that is the biology.

    * **TRDV1/2/3 are dedicated delta V genes** and rearrange to TRDJ. A ``TRDV1 + TRAJ`` scaffold
      is not a rearrangement that occurs, so building one invites reads onto a chimera.
    * **TRAV/DV genes pair with either**, and which J they took is what defines the locus. Both
      directions are already covered: ``TRAV/DV x TRAJ`` comes free with the TRA build (IMGT files
      those genes under TRAV), and ``TRAV/DV x TRDJ`` is what TRD's ``v_shared`` adds.

    Pinned because the symmetric-looking version was tried, and the thing that made it look
    justified was an IgBLAST truth calling 147 amplicon reads ``TRDV1*01`` + a TRAJ. arda declining
    to build that scaffold is arda being right about the biology, not a recall gap.
    """
    assert _locus("TRA").v_shared is None, (
        "TRA must not pull the TRDV stem: TRDV1/2/3 rearrange to TRDJ, and the shared TRAV/DV "
        "genes are already covered from both sides")


def test_allow_chimeras_adds_the_trdv_stem_to_tra_and_changes_nothing_else():
    """``--allow-chimeras`` is the ONLY way to get ``TRDV × TRAJ``, and it is opt-in.

    The default's refusal is a biological claim, and the external evidence disputes it: on 48,030
    TRA amplicon reads, IgBLAST calls ``TRDV1`` + a ``TRAJ`` on 530 (1.10 % of the library, median
    v_score 93.8, every one carrying a junction) and MiXCR independently agrees, while arda calls
    the same J as both and emits no ``v_call`` -- 83 % of its whole remaining ``v_gene`` gap
    (arda-benchmark ``results/round18`` §5c). Which side is right is a domain judgement, so it is a
    flag and not a silent default change.

    What this test pins is that the flag is *surgical*: TRA gains the TRDV stem, TRD keeps its
    ``TRAV/DV`` sharing, and no other locus moves.
    """
    from arda.refbuild.loci import loci_for

    assert loci_for() is LOCI, "the default must be the shipped table, not a rebuilt copy"

    chim = {loc.name: loc for loc in loci_for(allow_chimeras=True)}
    assert chim["TRA"].v_shared == ("TRDV", ""), "TRA must pull the WHOLE TRDV stem"
    assert chim["TRD"].v_shared == ("TRAV", "/DV"), "TRD's own sharing must be untouched"
    for loc in LOCI:
        if loc.name != "TRA":
            assert chim[loc.name] == loc, f"{loc.name} must be unchanged by --allow-chimeras"


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
