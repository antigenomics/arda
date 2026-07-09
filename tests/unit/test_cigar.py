"""Per-segment AIRR CIGAR construction (pure; no reference DB, no mmseqs).

The reference is the AIRR spec's CIGAR rules: a leading ``S`` (query 5' offset) then ``N`` (germline
5' offset) are required; the body is M/I/D; a trailing ``S`` is emitted, trailing ``N`` omitted.
"""

import re

from arda.annotate.cigar import segment_cigars, build_cigar, _classify


def _query_len_from_cigar(cig: str) -> int:
    """Query bases a CIGAR accounts for = sum of M/I/S run lengths (D and N are reference-side)."""
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDSN])", cig) if op in "MIS")


def test_spec_example_reproduced():
    """The AIRR spec's own D example: a match starting at query 419 (418S) and germline 11 (10N),
    16 nt, ending 71 nt from the query end. arda omits the optional trailing germline N (spec's 5N)."""
    assert build_cigar(q_lead=418, g_lead=10, ops=["M"] * 16, q_trail=71) == "418S10N16M71S"


def test_simple_v_j_split():
    """A gapless query over V (t1-6), the pad (7-9), then J (10-15). Both germlines start at 1, so
    no N; each segment M run with the rest of the read soft-clipped; the pad belongs to neither."""
    q = t = "ACGTACGTACGTACG"                       # 15 nt, identical -> all M
    out = segment_cigars(q, t, 1, 1, 15, t_vend=6, t_jstart=10, t_vjend=15)
    assert out == {"v_cigar": "6M9S", "j_cigar": "9S6M"}


def test_internal_deletion_and_insertion_in_v():
    out_del = segment_cigars("ACGT-CGTAC", "ACGTACGTAC", 1, 1, 9, t_vend=10, t_jstart=0, t_vjend=0)
    assert out_del == {"v_cigar": "4M1D5M"}           # target base 5 missing from the query
    out_ins = segment_cigars("ACGTXCGTAC", "ACGT-CGTAC", 1, 1, 10, t_vend=9, t_jstart=0, t_vjend=0)
    assert out_ins == {"v_cigar": "4M1I5M"}           # an extra query base vs the germline


def test_five_prime_germline_truncation_becomes_leading_N():
    """A read starting inside V (target 5) at its own base 1: the 4 missing 5' germline bases are a
    required leading ``4N`` -- NOT a soft-clip, and not silently carried only by v_germline_start."""
    out = segment_cigars("ACGTAC", "ACGTAC", 1, 5, 6, t_vend=10, t_jstart=0, t_vjend=0)
    assert out == {"v_cigar": "4N6M"}


def test_j_starting_mid_germline_gets_leading_N():
    """The BC100294.1 shape: a J alignment that begins at J germline position 56 (V-less read).
    q_start 200, so 199S; germline offset 55 -> 55N."""
    q = t = "ACGTACGT"                                # 8 nt all M, query 200..207
    out = segment_cigars(q, t, 200, 56, 1195, t_vend=0, t_jstart=1, t_vjend=100)
    assert out["j_cigar"] == "199S55N8M988S"           # 199 + 8 + 988 == 1195


def test_j_plus_c_scaffold_has_no_v_cigar():
    """A `J + C` hit (t_vend=0): J from its start, then C from its start, no V, no N."""
    q = t = "ACGTACGTACGTACGTACGT"                   # 20 nt all M
    out = segment_cigars(q, t, 1, 1, 20, t_vend=0, t_jstart=1, t_vjend=8)
    assert out == {"j_cigar": "8M12S", "c_cigar": "8S12M"}


def test_classify_boundaries():
    assert _classify(6, 6, 10, 15) == "v"             # V is inclusive to t_vend
    assert _classify(7, 6, 10, 15) is None            # pad: neither germline
    assert _classify(10, 6, 10, 15) == "j"
    assert _classify(16, 6, 10, 15) == "c"            # past the V-J end -> constant region


def test_cigar_accounts_for_the_whole_read():
    """AIRR invariant: each segment CIGAR's query-consuming ops (M/I/S) sum to the full read length,
    so a consumer can lay it over the read unambiguously (N and D are reference-side, excluded)."""
    q = "ACGTACGTAC" + "TTT" + "GGGGGCACGT"           # 23 nt, an insertion-y middle
    t = "ACGTACGTAC" + "---" + "GGGGGCACGT"           # 3-nt insertion between V and J
    out = segment_cigars(q, t, 1, 1, len(q), t_vend=10, t_jstart=11, t_vjend=20)
    assert out
    for cig in out.values():
        assert _query_len_from_cigar(cig) == len(q), f"{cig} does not cover all {len(q)} query nt"
