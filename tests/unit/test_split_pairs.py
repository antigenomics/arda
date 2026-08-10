"""Sharding paired FASTQ: mates together, quality intact, order preserved.

`split()` (the FASTA/amplicon path) does none of these: it drops quality and round-robins
records, so a fragment's mates land in different shards. Feeding that to the RNA-seq path
would recreate the mate-desync that produced a retracted finding in this project's own data.

`split_pairs()` cuts R1 and R2 at the same record boundaries into CONTIGUOUS blocks, which is
what makes a sharded run byte-identical to a single-node one rather than merely equivalent.

No DB, no mmseqs, no cluster.
"""

from __future__ import annotations

import gzip

import pytest

from arda.cluster import split_pairs


def _fastq(path, n: int, mate: str, *, start: int = 0, gz: bool = False):
    """`n` records with per-record distinct sequence AND quality."""
    recs = []
    for i in range(start, start + n):
        qual = "".join(chr(33 + ((i + j) % 60)) for j in range(10))  # non-trivial, varies per record
        recs.append(f"@frag{i}/{mate}\nACGTACGTAC\n+\n{qual}\n")
    data = "".join(recs)
    if gz:
        path.write_bytes(gzip.compress(data.encode()))
    else:
        path.write_text(data)
    return path


def _ids(path) -> list[str]:
    lines = path.read_text().splitlines()
    return [lines[i][1:] for i in range(0, len(lines), 4)]


def _records(path) -> list[str]:
    lines = path.read_text().splitlines()
    return ["\n".join(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def test_mates_are_never_separated(tmp_path):
    """The defect this function exists to prevent."""
    r1 = _fastq(tmp_path / "r1.fq", 10, "1")
    r2 = _fastq(tmp_path / "r2.fq", 10, "2")
    shards = split_pairs(r1, tmp_path / "sh", shards=3, r2=r2)

    for p1, p2 in shards:
        a = [i.rsplit("/", 1)[0] for i in _ids(p1)]
        b = [i.rsplit("/", 1)[0] for i in _ids(p2)]
        assert a == b, f"{p1.name} and {p2.name} hold different fragments"
        assert [i.endswith("/1") for i in _ids(p1)] == [True] * len(a)
        assert [i.endswith("/2") for i in _ids(p2)] == [True] * len(b)


def test_shards_are_contiguous_and_cover_every_pair_once(tmp_path):
    """Contiguity is what makes the merged AIRR identical, not merely a permutation."""
    r1 = _fastq(tmp_path / "r1.fq", 23, "1")
    r2 = _fastq(tmp_path / "r2.fq", 23, "2")
    shards = split_pairs(r1, tmp_path / "sh", shards=4, r2=r2)

    seen: list[int] = []
    for p1, _ in shards:
        nums = [int(i[len("frag"):-2]) for i in _ids(p1)]
        assert nums == list(range(nums[0], nums[0] + len(nums))), f"{p1.name} is not contiguous"
        seen += nums
    assert seen == list(range(23)), "concatenating shards in order must rebuild the input order"


def test_quality_strings_survive_byte_for_byte(tmp_path):
    """The anti-FASTA guard: `split()` would have thrown the quality away entirely."""
    r1 = _fastq(tmp_path / "r1.fq", 12, "1")
    r2 = _fastq(tmp_path / "r2.fq", 12, "2")
    original = _records(r1)

    shards = split_pairs(r1, tmp_path / "sh", shards=5, r2=r2)
    got = [rec for p1, _ in shards for rec in _records(p1)]
    assert got == original


def test_gzipped_input_is_accepted(tmp_path):
    r1 = _fastq(tmp_path / "r1.fq.gz", 8, "1", gz=True)
    r2 = _fastq(tmp_path / "r2.fq.gz", 8, "2", gz=True)
    shards = split_pairs(r1, tmp_path / "sh", shards=3, r2=r2)
    assert sum(len(_ids(p)) for p, _ in shards) == 8


def test_shard_names_sort_numerically(tmp_path):
    """`shard_10` must not sort before `shard_2` — the merge relies on sorted() being order."""
    r1 = _fastq(tmp_path / "r1.fq", 24, "1")
    shards = split_pairs(r1, tmp_path / "sh", shards=12)
    names = [p.name for p, _ in shards]
    assert names == sorted(names)
    assert "shard_00000_R1.fastq" in names and "shard_00011_R1.fastq" in names


def test_a_truncated_mate_is_caught_before_any_work(tmp_path):
    r1 = _fastq(tmp_path / "r1.fq", 10, "1")
    r2 = _fastq(tmp_path / "r2.fq", 9, "2")
    with pytest.raises(ValueError, match="truncated"):
        split_pairs(r1, tmp_path / "sh", shards=3, r2=r2)


def test_more_shards_than_pairs_writes_no_empty_shard(tmp_path):
    """An empty shard would make its array task die in `detect_format` and fail the afterok."""
    r1 = _fastq(tmp_path / "r1.fq", 3, "1")
    r2 = _fastq(tmp_path / "r2.fq", 3, "2")
    shards = split_pairs(r1, tmp_path / "sh", shards=7, r2=r2)
    assert len(shards) == 3
    for p1, p2 in shards:
        assert p1.stat().st_size > 0 and p2.stat().st_size > 0
    assert sorted(p.name for p in (tmp_path / "sh").glob("*_R1.fastq")) == [
        "shard_00000_R1.fastq", "shard_00001_R1.fastq", "shard_00002_R1.fastq"]


def test_single_end_yields_no_r2(tmp_path):
    r1 = _fastq(tmp_path / "r1.fq", 9, "1")
    shards = split_pairs(r1, tmp_path / "sh", shards=3)
    assert all(p2 is None for _, p2 in shards)
    assert sum(len(_ids(p1)) for p1, _ in shards) == 9


def test_fasta_input_is_refused_and_names_the_right_command(tmp_path):
    fa = tmp_path / "in.fasta"
    fa.write_text(">s0\nACGT\n>s1\nTGCA\n")
    with pytest.raises(ValueError, match="arda cluster split-fasta"):
        split_pairs(fa, tmp_path / "sh", shards=2)


def test_rejects_a_nonsense_shard_count(tmp_path):
    r1 = _fastq(tmp_path / "r1.fq", 4, "1")
    for bad in (0, -1, 10 ** 6):
        with pytest.raises(ValueError, match="shards must be"):
            split_pairs(r1, tmp_path / "sh", shards=bad)


def test_one_shard_reproduces_the_input_exactly(tmp_path):
    r1 = _fastq(tmp_path / "r1.fq", 15, "1")
    (p1, _), = split_pairs(r1, tmp_path / "sh", shards=1)
    assert p1.read_text() == r1.read_text()
