"""The SHM scope: the mutation lists and `v_identity` must stop at Cys104, not at the segment end.

⛔ This is the retraction of a guarantee 2.14.0 printed in `docs/shm.rst`. Segment scoping keeps
the scaffold's N-pad out of the lists, and that is NOT the same as keeping the junction out: a
rearranged junction is *V 3' tail + N/P + J 5' head*, so both germlines' templated tails are inside
it and chew-back reads as a substitution against a germline that does not template it.

Every test here is a record-level property, so it needs no reference, no aligner and no fixture —
which is also why `arda shm` can recount a TSV that was written before the fix existed.
"""

import pytest

from arda.shm import FULL_COLUMNS, scope_record


def _rec(**over):
    """A V germline whose Cys104 codon starts at 0-based offset 100, and a J whose [FW]118 is at 8.

    So V positions 1..100 are framework and 101+ are junction; J positions 1..11 are junction and
    12+ are FR4.
    """
    rec = {
        "v_mutations": "G45A,C100T,A101G,T150C",
        "j_mutations": "G1A,C11T,A12G,T30C",
        "v_identity": 0.95,
        "v_anchor_nt": 100,
        "j_anchor_nt": 8,
        "sequence_alignment": "",
        "germline_alignment": "",
        "mmseqs2_tstart": 1,
        "mmseqs2_t_vend": 300,
    }
    rec.update(over)
    return rec


def test_framework_drops_everything_past_the_v_anchor():
    r = scope_record(_rec())
    assert r["v_mutations"] == "G45A,C100T"          # 101 and 150 are inside the junction


def test_framework_drops_everything_up_to_and_including_the_j_anchor_codon():
    """[FW]118 occupies germline positions `j_anchor_nt + 1 .. + 3`, and it is the LAST junction
    codon — so position 11 is still junction and 12 is the first FR4 base."""
    r = scope_record(_rec())
    assert r["j_mutations"] == "A12G,T30C"


def test_a_missing_anchor_leaves_the_list_ALONE_rather_than_emptying_it():
    """⛔ No anchor means arda does not know where this germline's junction starts. Emptying the
    list would be a claim; leaving it is the honest failure, and the raw anchors ship so a consumer
    can tell the two apart."""
    r = scope_record(_rec(v_anchor_nt="", j_anchor_nt=""))
    assert r["v_mutations"] == "G45A,C100T,A101G,T150C"
    assert r["j_mutations"] == "G1A,C11T,A12G,T30C"


def test_both_keeps_the_scoped_values_in_the_shipped_columns_and_adds_the_old_ones():
    """⛔ `v_identity`/`v_mutations` mean the SAME thing in every mode. `both` ADDS the legacy
    numbers under new names — a column whose meaning depends on a flag is unreadable downstream."""
    r = scope_record(_rec(), "both")
    assert r["v_mutations"] == "G45A,C100T"
    assert r["v_mutations_full"] == "G45A,C100T,A101G,T150C"
    assert r["j_mutations_full"] == "G1A,C11T,A12G,T30C"
    assert set(FULL_COLUMNS) <= set(r)


def test_off_emits_no_shm_fields():
    r = scope_record(_rec(), "off")
    assert r["v_mutations"] == "" and r["j_mutations"] == "" and r["v_identity"] == ""


def test_an_unknown_mode_raises_rather_than_falling_through_to_the_default():
    with pytest.raises(ValueError, match="shm mode"):
        scope_record(_rec(), "frameworkk")


def test_v_identity_is_remeasured_over_the_framework_only():
    """The V part of a scaffold IS the V germline verbatim, so the framework cut is the anchor in
    target coordinates. Here the read matches germline over positions 1-4 and mismatches at 5-6,
    with the anchor at 4: the whole-V identity is 4/6, the framework identity is 1.0.
    """
    calls = []

    def identity_fn(qaln, taln, tstart, t_lo, t_hi):
        # The same walk the real `_aln_identity` does, in four lines.
        calls.append((tstart, t_lo, t_hi))
        m = n = 0
        t = tstart
        for q, g in zip(qaln, taln):
            if g != "-":
                if t_lo <= t <= t_hi:
                    n += 1
                    m += q == g
                t += 1
        return m / n if n else ""

    r = scope_record(
        _rec(sequence_alignment="ACGTAA", germline_alignment="ACGTCC",
             mmseqs2_tstart=1, mmseqs2_t_vend=6, v_anchor_nt=4, v_identity=4 / 6),
        identity_fn=identity_fn)
    assert calls == [(1, 1, 4)]                      # never past the anchor
    assert r["v_identity"] == 1.0


# --- the standalone `arda shm` stage -------------------------------------------------------------

def _airr(tmp_path, cols):
    import polars as pl

    p = tmp_path / "in.airr.tsv"
    pl.DataFrame([cols]).write_csv(p, separator="\t", quote_style="never")
    return p


def test_recount_rescopes_an_existing_tsv_without_a_reference(tmp_path):
    """⛔ The whole point: `v_anchor_nt`/`j_anchor_nt` and the alignment strings are already in the
    file, so a table written before this existed can be fixed without re-mapping it."""
    import polars as pl

    from arda.shm import recount_airr

    src = _airr(tmp_path, _rec(sequence_alignment="", germline_alignment="",
                               sequence_id="r1", v_identity=0.95))
    out = tmp_path / "out.tsv"
    rep = recount_airr(src, out)
    assert rep["mutations_in"] == 8 and rep["mutations_out"] == 4 and rep["removed"] == 4
    got = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert got["v_mutations"][0] == "G45A,C100T"
    assert got["j_mutations"][0] == "A12G,T30C"


def test_recount_RAISES_on_a_file_that_predates_the_anchor_columns(tmp_path):
    """⛔ Not a pass-through with a success message. `resolve_airr` shipped exactly that failure —
    it caught the error per locus, continued with an empty germline set, and returned output
    byte-identical to its input while reporting success."""
    from arda.shm import recount_airr

    src = _airr(tmp_path, {"sequence_id": "r1", "v_mutations": "G45A", "j_mutations": ""})
    with pytest.raises(ValueError, match="v_anchor_nt"):
        recount_airr(src, tmp_path / "out.tsv")
