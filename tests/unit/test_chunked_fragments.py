"""A chunk boundary must never fall between a fragment's two mates.

`_apply_constant_rule` runs per chunk and pairs a CDR3-bearing read with its constant-region
mate to donate an isotype. Split the mates into different chunks and that donation is lost --
silently, as a slightly lower `isotype_from_mate` count. It also makes `map` output a function
of `--chunk-size`, and under sharding a function of the shard layout, which is precisely the
byte-identity the SLURM and Nextflow paths must guarantee.

No DB, no mmseqs: `chunked_fragments` is a pure generator over (id, seq) tuples.
"""

from __future__ import annotations

from arda.rnaseq.map import chunked_fragments, frag_stem


def _pairs(n: int) -> list[tuple[str, str]]:
    """`n` fragments, two mates each, in the order `read_pairs` emits them."""
    out = []
    for i in range(n):
        out += [(f"f{i}/1", "ACGT"), (f"f{i}/2", "TGCA")]
    return out


def _no_fragment_is_split(chunks: list[list[tuple[str, str]]]) -> None:
    seen: set[str] = set()
    for chunk in chunks:
        stems = {frag_stem(i) for i, _ in chunk}
        assert not (stems & seen), f"fragment(s) {stems & seen} appear in more than one chunk"
        seen |= stems


def test_a_pair_is_never_split_even_when_the_size_forces_it():
    """size=1 would put every mate in its own chunk; fragments must still stay whole."""
    chunks = list(chunked_fragments(iter(_pairs(3)), 1))
    _no_fragment_is_split(chunks)
    assert chunks == [
        [("f0/1", "ACGT"), ("f0/2", "TGCA")],
        [("f1/1", "ACGT"), ("f1/2", "TGCA")],
        [("f2/1", "ACGT"), ("f2/2", "TGCA")],
    ]


def test_odd_chunk_size_against_two_mates_per_fragment():
    """size=3 straddles the 2-records-per-fragment rhythm on every other boundary."""
    chunks = list(chunked_fragments(iter(_pairs(6)), 3))
    _no_fragment_is_split(chunks)
    assert sum(len(c) for c in chunks) == 12


def test_reconstruct_style_mixed_record_counts_stay_whole():
    """With --reconstruct a merged pair emits ONE record and an unmerged pair TWO.

    That drifting parity is what breaks plain `chunked`: boundaries stop landing between
    fragments by luck. Interleave the two shapes and require whole fragments anyway.
    """
    recs = [
        ("m0", "ACGTACGT"),                       # merged -> 1 record
        ("f1/1", "ACGT"), ("f1/2", "TGCA"),       # unmerged -> 2
        ("m2", "ACGTACGT"),
        ("f3/1", "ACGT"), ("f3/2", "TGCA"),
        ("f4/1", "ACGT"), ("f4/2", "TGCA"),
    ]
    for size in (1, 2, 3, 4, 5):
        chunks = list(chunked_fragments(iter(recs), size))
        _no_fragment_is_split(chunks)
        assert [r for c in chunks for r in c] == recs, f"records reordered or lost at size={size}"


def test_every_record_survives_in_order_for_any_size():
    recs = _pairs(17)
    for size in range(1, 40):
        chunks = list(chunked_fragments(iter(recs), size))
        _no_fragment_is_split(chunks)
        assert [r for c in chunks for r in c] == recs, f"size={size}"


def test_single_end_records_are_all_distinct_fragments():
    """Single-end input has no mate suffix; each record is its own fragment."""
    recs = [(f"s{i}", "ACGT") for i in range(10)]
    chunks = list(chunked_fragments(iter(recs), 4))
    assert [len(c) for c in chunks] == [4, 4, 2]
    assert [r for c in chunks for r in c] == recs


def test_empty_input_yields_nothing():
    assert list(chunked_fragments(iter([]), 5)) == []


def test_frag_stem_strips_mate_suffix_and_illumina_comment():
    assert frag_stem("READ1/1") == "READ1"
    assert frag_stem("READ1/2") == "READ1"
    assert frag_stem("READ1 1:N:0:ATCG") == "READ1"
    assert frag_stem("READ1") == "READ1"
