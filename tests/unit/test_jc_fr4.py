"""FR4 on `J + C` scaffolds — the derivation and its shipped result.

A `J + C` scaffold has no V, so it cannot be annotated by IgBLAST and every V-anchored region
stays -1. FR4 is the one exception: it lies wholly inside J, and IgBLAST's aux file states where
(columns 4 = CDR3 stop, 5 = extra bp beyond the J coding end). ``_fr4_span`` turns those into a
J-local FR4 span; ``build_jc_scaffolds`` writes it into the scaffold.

These tests cover the pure span math (offline) and then guard the shipped markup against drift.
"""

from __future__ import annotations

import pytest

from arda import paths
from arda.refbuild.constant import _fr4_span, build_jc_scaffolds
from arda.refbuild.combinations import load_j_fr4_offsets

# ---- pure span math: no reference DB, no aux file, no mmseqs -------------------------------------


def test_fr4_span_single_allele_round_trips():
    """``fwr4 = jseq[cdr3_stop+1 : len - extra_bp]``, returned 1-based closed."""
    jseq = "AAA" + "TGGGGCCAAGGGACCCTGGTCACCGTC" + "CCC"   # 3 + 27 + 3 = 33 nt
    s, e = _fr4_span(jseq, ["IGXJ9*01"], {"IGXJ9*01": (2, 3)})   # 0-based stop 2, 3 extra
    assert (s, e) == (4, 30)                                     # 1-based; drops the 3-nt tail
    assert jseq[s - 1:e] == jseq[3:30]


def test_fr4_span_extra_bp_zero_is_valid():
    """IGHJ6*02 has ``extra_bp = 0`` and a 34-nt FR4 -- IgBLAST's own convention, not an error.
    A zero tail must give a span that runs to the J terminus."""
    jseq = "GG" + "TGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG"          # stop after 2 nt, no tail
    s, e = _fr4_span(jseq, ["IGHJ6*02"], {"IGHJ6*02": (1, 0)})
    assert (s, e) == (3, len(jseq))
    assert jseq[s - 1:e] == jseq[2:]


def test_fr4_span_refuses_when_collapsed_alleles_disagree():
    """Alleles collapsed onto one identical J sequence must share an FR4 offset. If the aux file
    disagrees between them, report nothing -- a wrong FR4 silently shifts fwr4_aa's frame, which is
    worse than an absent one."""
    jseq = "ACGT" * 10
    assert _fr4_span(jseq, ["JA*01", "JB*01"], {"JA*01": (2, 1), "JB*01": (5, 1)}) == (-1, -1)
    # ...but identical offsets on both is fine.
    assert _fr4_span(jseq, ["JA*01", "JB*01"], {"JA*01": (2, 1), "JB*01": (2, 1)}) != (-1, -1)


def test_fr4_span_refuses_when_no_allele_has_an_offset():
    """Pseudogene J entries carry only three aux columns, so they never enter the offset map."""
    assert _fr4_span("ACGT" * 10, ["IGHJ1P*01"], {}) == (-1, -1)


def test_fr4_span_refuses_an_out_of_range_span():
    """A garbage offset that would run past the J must not produce a bogus span."""
    jseq = "ACGT" * 5                                            # 20 nt
    assert _fr4_span(jseq, ["J*01"], {"J*01": (25, 0)}) == (-1, -1)   # stop beyond the sequence
    assert _fr4_span(jseq, ["J*01"], {"J*01": (2, 25)}) == (-1, -1)   # tail longer than what's left


# ---- the shipped reference: does the derivation actually fire, and stay consistent? --------------

_HAS_BUNDLE = (paths.database_dir() / "c_genes" / "human.fasta").exists()
_HAS_AUX = bool(load_j_fr4_offsets("human"))
requires_jc_inputs = pytest.mark.skipif(
    not (_HAS_BUNDLE and _HAS_AUX), reason="C-gene bundle or IgBLAST aux file not present")


@pytest.mark.skipif(not _HAS_AUX, reason="IgBLAST aux file not present")
def test_load_j_fr4_offsets_parses_both_columns_from_the_real_aux():
    """Pin ``load_j_fr4_offsets`` to ground truth in bin/optional_file/human_gl.aux -- BOTH column 4
    (cdr3_stop) and column 5 (extra_bp). Every other FR4 test either reads the shipped markup (built
    with correct code, so blind to a parser regression) or feeds ``_fr4_span`` offsets directly, so
    without this a change that drops or mis-indexes column 5 ships silently. ``extra_bp`` matters:
    IGHJ6*02 has extra_bp=0 (34-nt FR4) vs IGHJ4*02's 1, and hardcoding it to 0 passes all other
    tests."""
    off = load_j_fr4_offsets("human")
    assert off["IGHJ4*02"] == (13, 1)
    assert off["IGHJ6*02"] == (28, 0)      # extra_bp genuinely 0 here -- the discriminating case
    assert off["IGKJ1*01"] == (6, 1)
    assert "IGHJ1P*01" not in off           # a 3-column pseudogene row has no FR4 offset


@requires_jc_inputs
@pytest.mark.parametrize("organism,species_dir", [("human", "Homo_sapiens"), ("mouse", "Mus_musculus")])
def test_built_jc_scaffolds_carry_an_in_bounds_fr4(organism, species_dir):
    """Most J+C scaffolds get an FR4, and where they do it is an in-bounds span of a plausible FR4
    length (~10-11 codons). The *string's* correctness is guarded byte-for-byte against the V-J
    scaffolds below; this only checks the span geometry.

    Deliberately no amino-acid assertion. ``TRAJ35*01`` opens FR4 with a real germline **C** (not the
    canonical [WF]), and ~24 mouse TRAJ/TRGJ alleles do not translate cleanly at frame 0 of the
    isolated slice -- arda translates ``fwr4_aa`` from slice-frame-0, which is not the J coding frame
    for those genes (a pre-existing property of the V-J path too; see the module task note)."""
    jc = build_jc_scaffolds(organism, species_dir)
    with_fr4 = [s for s in jc if s.fwr4_start > 0]
    assert len(with_fr4) >= 0.5 * len(jc), f"{organism}: only {len(with_fr4)}/{len(jc)} got an FR4"
    for s in with_fr4:
        assert 1 <= s.fwr4_start <= s.fwr4_end <= len(s.sequence)
        fr4 = s.sequence[s.fwr4_start - 1:s.fwr4_end]
        assert fr4 == s.sequence[s.fwr4_start - 1:s.fwr4_end]     # the stored coords, round-tripped
        assert 24 <= len(fr4) <= 45, f"{organism}/{s.scaffold_id}: implausible FR4 length {len(fr4)}"


@requires_jc_inputs
def test_fr4_is_absent_exactly_when_the_aux_has_no_agreeing_offset():
    """FR4 is emitted iff the aux file gives one -- ``build_jc_scaffolds`` must agree with the pure
    ``_fr4_span``. Some human J+C scaffolds legitimately have no FR4 (a J allele the aux does not
    cover, e.g. a pseudogene J with only three columns); those must stay -1, never borrow a frame."""
    offsets = load_j_fr4_offsets("human")
    jc = build_jc_scaffolds("human", "Homo_sapiens")
    without = [s for s in jc if s.fwr4_start == -1]
    assert without, "expected some human J+C scaffolds with no aux FR4 (pseudogene / uncovered J)"
    for s in jc:
        jseq = s.sequence[:s.j_len]
        expect = _fr4_span(jseq, s.j_call.split(","), offsets)
        assert (s.fwr4_start, s.fwr4_end) == expect, f"{s.scaffold_id}: build vs _fr4_span disagree"


@pytest.mark.skipif(not (paths.vdj_dir("human") / "markup.tsv").exists(),
                    reason="human reference DB not built")
def test_shipped_jc_fr4_matches_the_vj_scaffolds_byte_for_byte():
    """The claim behind the fix: FR4 is the same string whether read off a `J + C` scaffold or a
    V-J one. For every J allele present in both kinds, the two FR4 strings must be identical --
    otherwise a J->C read and a V-J read of one clone disagree on their only shared region."""
    import polars as pl

    df = pl.read_csv(paths.vdj_dir("human") / "markup.tsv", separator="\t", infer_schema_length=0)
    vj_fr4: dict[str, set[str]] = {}
    for r in df.filter((pl.col("v_call") != "") & (pl.col("fwr4") != "")).iter_rows(named=True):
        for al in r["j_call"].split(","):
            vj_fr4.setdefault(al, set()).add(r["fwr4"])

    checked = 0
    for r in df.filter(pl.col("c_call") != "").iter_rows(named=True):
        if not r["fwr4"]:
            continue
        for al in r["j_call"].split(","):
            if al in vj_fr4:
                assert r["fwr4"] in vj_fr4[al], \
                    f"{al}: J+C FR4 {r['fwr4']!r} not among V-J FR4s {vj_fr4[al]}"
                checked += 1
    assert checked >= 50, f"only {checked} J alleles were cross-checked — fixture may have drifted"
