"""The quality-aware denoising framework: `arda.rnaseq.denoise` + `--clonotype-key`.

⛔ The invariant these tests exist to defend: **denoising MOVES reads, it never discards them.** A
read that reached a complete junction came off a real rearrangement of that locus, so deciding its
junction carries a miscall is a statement about bases, not about whether the molecule existed.
Every regime is checked for read conservation, not just for clonotype counts -- a clonotype count
falling is the point, a read count falling is a bug.
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.rnaseq.correct import CLONOTYPE_KEYS, EC_MODES, correct_airr
from arda.rnaseq.denoise import (
    REGIMES,
    DenoiseParams,
    _nearest_py,
    _read_quality_py,
    clonotype_quality,
    quality_rescue,
    read_quality,
)

try:
    from arda import _denoise as _cpp
except ImportError:                                # pragma: no cover
    _cpp = None


# --- the C++ core, against its Python reference -------------------------------------------------

def test_read_quality_matches_the_python_reference():
    jn = ["ACGT", "ACGTACGT", "ACGT", "ACGT", ""]
    q = ["IIII", "5555!!!!", "", "II", ""]
    assert read_quality(jn, q) == pytest.approx(_read_quality_py(jn, q))


def test_missing_or_mismatched_quality_is_absent_evidence_not_bad_evidence():
    """⛔ -1.0, never 0.0. A quality string of the wrong length is the one corruption nothing
    downstream can detect (right length, wrong strand or offset), so it is refused, not averaged."""
    assert read_quality(["ACGT"], [""])[0] == -1.0
    assert read_quality(["ACGT"], ["II"])[0] == -1.0
    assert read_quality(["ACGT"], ["!!!!"])[0] == 0.0        # genuinely Q0 is 0.0, not -1.0


@pytest.mark.skipif(_cpp is None, reason="extension not built")
def test_nearest_more_abundant_matches_the_python_reference():
    seqs = ["AAAAAA", "AAAAAT", "AAAATT", "TTTTTT", "AAAAAA"]
    counts = [1000, 5, 2, 1, 3]
    cand = [False, True, True, True, True]
    for k in (1, 2, 3, 6):
        got = _cpp.nearest_more_abundant(seqs, counts, cand, k, 1.0)
        assert list(got) == _nearest_py(seqs, counts, cand, k, 1.0), k


def test_a_parent_must_be_strictly_more_abundant():
    """Equal counts must not make two clonotypes each other's parent -- the root walk would hang."""
    seqs, counts, cand = ["AAAA", "AAAT"], [7, 7], [True, True]
    assert _nearest_py(seqs, counts, cand, 1, 1.0) == [-1, -1]


def test_the_abundance_ratio_gates_the_rescue():
    seqs, counts = ["AAAAAA", "AAAATT"], [100, 3]
    cand = [False, True]
    assert _nearest_py(seqs, counts, cand, 2, 10.0) == [-1, 0]     # 100 >= 3*10 -> rescued
    assert _nearest_py(seqs, counts, cand, 2, 50.0) == [-1, -1]    # 100 < 3*50 -> left alone


def test_only_equal_length_sequences_are_compared():
    """A length change is an indel -- a different event at a different rate. Comparing across
    lengths silently aligns junctions that do not correspond."""
    # One entry per sequence; the candidate at index 0 finds nothing because the only
    # more-abundant clonotype is a different length.
    assert _nearest_py(["AAAA", "AAA"], [1, 500], [True, False], 3, 1.0) == [-1, -1]


# --- the framework ------------------------------------------------------------------------------

def test_clonotype_quality_uses_the_median_and_ignores_unusable_reads():
    keys = ["a", "a", "a", "b"]
    q = [40.0, 10.0, 38.0, -1.0]
    cq = clonotype_quality(keys, q)
    assert cq["a"] == 38.0, "one catastrophic read must not drag a clean clonotype down"
    assert "b" not in cq, "a clonotype with no usable read is never a rescue candidate"


def test_rescue_is_off_unless_both_knobs_are_set():
    p = DenoiseParams(lowq_mean_q=25.0, lowq_max_subs=0)
    assert not p.enabled()
    parent, rep = quality_rescue(["AAAA"], [1], [10.0], p)
    assert parent == [None] and rep.rescued_clonotypes == 0


def test_a_low_quality_clonotype_with_no_parent_keeps_its_reads():
    """⛔ The orphan case, and the reason the framework is a rescue radius and not a filter. On a
    polyclonal repertoire a mean-Q floor at 30 would strand 3.70 % of all reads this way."""
    # 16 subs apart -- outside `amplicon`'s 12-substitution rescue radius, so there is no parent
    # to route to even though the read quality says the clonotype is junk.
    seqs = ["AAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTT"]
    parent, rep = quality_rescue(seqs, [500, 1], [40.0, 10.0], REGIMES["amplicon"])
    assert parent[1] is None
    assert rep.orphan_clonotypes == 1 and rep.orphan_reads == 1
    assert rep.rescued_clonotypes == 0


def test_a_good_quality_far_neighbour_is_never_collapsed():
    """The whole point: distance alone is not evidence. Only bad reads may be rescued."""
    seqs = ["AAAAAAAAAAAA", "AAAAAATTTTTT"]      # 6 subs
    parent, _ = quality_rescue(seqs, [5000, 1], [40.0, 40.0], REGIMES["amplicon"])
    assert parent[1] is None, "a high-quality variant is a variant, not an error"
    parent, _ = quality_rescue(seqs, [5000, 1], [40.0, 12.0], REGIMES["amplicon"])
    assert parent[1] == 0, "the same variant with bad reads IS rescued"


def test_every_shipped_regime_is_a_known_ec_mode_and_vice_versa():
    """One table, not two: `EC_MODES` names a regime and `REGIMES` must define it."""
    for mode, cfg in EC_MODES.items():
        assert cfg["regime"] in REGIMES, mode
    assert set(REGIMES) == set(EC_MODES)


# --- end to end, on a fixture that carries every class the framework distinguishes ---------------

_PARENT = "TGTGCCAGCAGCTTAGACGGGACAGGGTTC"
_SPLIT_V = "TRBV6-6*01"                     # a CALL SPLIT: same junction, different V
_NEAR = "TGTGCCAGCAGCTTAGACGGGACAGGGTTT"    # 1 sub -- the abundance model's ladder
_FAR = "TGTGCCAGCAGCTTAGACGGGTAAGCATTC"    # 5 subs -- the cliff, no ladder behind it


def _rows(junction, n, qual, *, prefix, v="TRBV20-1*01"):
    return [{"sequence_id": f"{prefix}{k}", "sequence": junction, "junction": junction,
             "junction_aa": "CASSLDGTGF", "v_call": v, "j_call": "TRBJ2-1*01", "locus": "TRB",
             "junction_quality": qual} for k in range(n)]


def _fixture(tmp_path, name="in.tsv"):
    q_good = "I" * len(_PARENT)
    q_bad = "5" * len(_PARENT)              # Q20 -- below every rescue threshold
    rows = (_rows(_PARENT, 4000, q_good, prefix="p")
            + _rows(_PARENT, 30, q_good, prefix="s", v=_SPLIT_V)   # call split, GOOD reads
            + _rows(_NEAR, 3, q_good, prefix="n")                  # near neighbour
            + _rows(_FAR, 2, q_bad, prefix="f"))                   # far + bad = the cliff class
    p = tmp_path / name
    pl.DataFrame(rows).write_csv(p, separator="\t")
    return p


def _totals(path):
    df = pl.read_csv(path, separator="\t")
    return df.height, int(df["duplicate_count"].sum())


@pytest.mark.parametrize("mode", sorted(EC_MODES))
@pytest.mark.parametrize("key", CLONOTYPE_KEYS)
def test_no_regime_or_key_ever_loses_a_read(tmp_path, mode, key):
    """⛔ THE invariant. Clonotypes may fall -- that is the job. Reads may not."""
    src = _fixture(tmp_path)
    out = tmp_path / f"{mode}_{key}.tsv"
    correct_airr(src, out, map_d=False, ec_mode=mode, clonotype_key=key, error_rate=1e-6)
    _, reads = _totals(out)
    assert reads == 4035, f"{mode}/{key} moved the read total"


def test_the_junction_key_collapses_a_call_split_that_no_error_model_can_see(tmp_path):
    """A junction byte-identical to an abundant clone's under a different V has no discriminating
    base, so `--min-junction-q` and the abundance model are both blind to it by construction."""
    src = _fixture(tmp_path)
    full = tmp_path / "full.tsv"
    jn = tmp_path / "jn.tsv"
    correct_airr(src, full, map_d=False, ec_mode="fast", clonotype_key="full", error_rate=1e-6)
    correct_airr(src, jn, map_d=False, ec_mode="fast", clonotype_key="junction", error_rate=1e-6)
    n_full, r_full = _totals(full)
    n_jn, r_jn = _totals(jn)
    assert n_jn == n_full - 1, "the call split must merge"
    assert r_jn == r_full, "...and carry its reads with it"
    top = pl.read_csv(jn, separator="\t").sort("duplicate_count", descending=True)
    assert int(top["duplicate_count"][0]) >= 4030, "the split's reads land on the dominant clone"


def test_the_rescue_reaches_the_cliff_class_the_gate_cannot(tmp_path):
    """5 subs from the parent: no ladder, no discriminating base to gate on, but bad reads."""
    src = _fixture(tmp_path)
    acc = tmp_path / "acc.tsv"
    amp = tmp_path / "amp.tsv"
    r_acc = correct_airr(src, acc, map_d=False, ec_mode="accurate", error_rate=1e-6)
    r_amp = correct_airr(src, amp, map_d=False, ec_mode="amplicon", error_rate=1e-6)
    assert r_amp.rescued_clonotypes >= 1
    assert r_acc.rescued_clonotypes == 0, "`accurate` must not silently gain the rescue"
    assert _totals(amp)[0] < _totals(acc)[0], "the cliff clonotype is absorbed"
    assert _totals(amp)[1] == _totals(acc)[1], "...without losing a read"


def test_an_unknown_clonotype_key_raises(tmp_path):
    with pytest.raises(ValueError, match="clonotype_key"):
        correct_airr(_fixture(tmp_path), tmp_path / "o.tsv", map_d=False, clonotype_key="nope")


def test_a_q1_base_in_the_quality_string_does_not_break_the_parse(tmp_path):
    """⛔ Phred+33 chr 34 is `"`, i.e. Q1 — a legitimate score any low-quality base produces.

    polars' CSV reader treats it as a quote character, so ONE such base collapsed the parse of the
    whole file (`CSV malformed: expected 1 rows, actual 155 rows`). Measured on a real Raji run:
    exactly one row of the file contained a `"` and the entire table became unreadable, which took
    out four of that sample's benchmark legs. Nothing in an AIRR TSV is ever quoted — the writer is
    `_markup.format_rows`, which emits raw fields — so quoting is disabled on read.
    """
    q_good = "I" * len(_PARENT)
    q_q1 = 'I' * 10 + '"' + "I" * (len(_PARENT) - 11)      # one Q1 base, mid-junction
    assert len(q_q1) == len(_PARENT)
    rows = _rows(_PARENT, 50, q_good, prefix="p") + _rows(_NEAR, 2, q_q1, prefix="q")
    src = tmp_path / "q1.tsv"
    pl.DataFrame(rows).write_csv(src, separator="\t", quote_style="never")
    out = tmp_path / "o.tsv"
    correct_airr(src, out, map_d=False, ec_mode="accurate", error_rate=1e-6)
    assert _totals(out)[1] == 52, "every read must survive a Q1 base in the quality string"


def test_the_rescue_never_merges_across_loci(tmp_path):
    """⛔ A rearrangement of a DIFFERENT locus is not a sequencing error of this one.

    `quality_rescue` groups candidate parents by junction LENGTH and nothing else -- there is no
    locus guard in `_nearest_py` or in the C++ `nearest_more_abundant` -- while ``--ec-mode
    amplicon`` opens the radius to 12 substitutions. Measured on SRR5233636 at full depth: of 9,025
    rescues, 3 crossed the locus, all of them 1-read TRB clonotypes absorbed into abundant TRA
    clonotypes at 11-12 substitutions. ``correct_airr`` now partitions the search by locus.

    Here a low-quality 1-read TRB clonotype sits 3 substitutions from a 5,000-read TRA clonotype and
    nothing else. It must stay put and KEEP ITS READ.
    """
    import polars as pl

    from arda.rnaseq.correct import correct_airr

    tra = "TGTGCCAGCAGTTTCTCGACCTGTTCGGCTAACTATGGCTACACCTTC"
    trb = tra[:10] + "".join("C" if c != "C" else "G" for c in tra[10:13]) + tra[13:]  # 3 subs
    assert len(tra) == len(trb) and sum(a != b for a, b in zip(tra, trb)) == 3

    def rows(jn, n, qual, locus, vc, jc, prefix):
        return [{"sequence_id": f"{prefix}{k}", "sequence": jn, "junction": jn,
                 "junction_aa": "CASSLDGTGF", "v_call": vc, "j_call": jc, "locus": locus,
                 "junction_quality": qual} for k in range(n)]

    data = (rows(tra, 5000, "I" * len(tra), "TRA", "TRAV1-2*01", "TRAJ33*01", "a")
            + rows(trb, 1, "#" * len(trb), "TRB", "TRBV12-3*01", "TRBJ1-2*01", "b"))  # '#' = Q2
    src = tmp_path / "xloc.tsv"
    pl.DataFrame(data).write_csv(src, separator="\t")
    out = tmp_path / "xloc.out"
    rep = correct_airr(src, out, map_d=False, ec_mode="amplicon")
    got = pl.read_csv(out, separator="\t").to_dicts()

    loci = {r["locus"]: r["duplicate_count"] for r in got}
    assert "TRB" in loci, "the TRB clonotype was absorbed into a TRA clonotype"
    assert loci["TRB"] == 1
    assert rep.reads_assigned == sum(r["duplicate_count"] for r in got) == 5001
