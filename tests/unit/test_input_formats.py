"""FASTA input and `--limit`, on every path a user can reach.

Both features already worked; neither was covered, which is how a working feature regresses without
anyone noticing. The paths differ in a way that matters: `read_pairs` hands the unlimited case to
**dnaio** and the limited case to a pure-Python reader (see its own comment on why), so a test that
only exercises one of them tests half the feature. Every case here asserts an exact record count.
"""
from __future__ import annotations

import gzip

import pytest

from arda.annotate import io as seqio
from arda.rnaseq.map import read_pairs

R1 = [("r1", "ACGTACGTACGTACGTAAAA"), ("r2", "TTTTGGGGCCCCAAAATTTT"),
      ("r3", "GGGGCCCCAAAATTTTGGGG"), ("r4", "CCCCAAAATTTTGGGGCCCC")]
R2 = [(i, s[::-1]) for i, s in R1]


def _write(path, body, gz):
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return path


def _fastq(path, recs, gz=False):
    return _write(path, "".join(f"@{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in recs), gz)


def _fasta(path, recs, gz=False):
    return _write(path, "".join(f">{i}\n{s}\n" for i, s in recs), gz)


# ---------------------------------------------------------------- format detection

@pytest.mark.parametrize("writer, want", [(_fasta, "fasta"), (_fastq, "fastq")])
@pytest.mark.parametrize("gz", [False, True])
def test_format_is_detected_by_content_not_extension(tmp_path, writer, want, gz):
    """The suffix is a hint, not the answer: `.fq` holding FASTA is a real thing users produce."""
    p = writer(tmp_path / ("x.dat.gz" if gz else "x.dat"), R1, gz=gz)
    assert seqio.detect_format(p) == want


# ---------------------------------------------------------------- FASTA input

@pytest.mark.parametrize("gz", [False, True])
def test_single_end_fasta(tmp_path, gz):
    p = _fasta(tmp_path / ("a.fa.gz" if gz else "a.fa"), R1, gz=gz)
    got = list(read_pairs(p))
    assert [i for i, _ in got] == ["r1", "r2", "r3", "r4"]
    assert got[0][1] == R1[0][1]


@pytest.mark.parametrize("gz", [False, True])
def test_paired_fasta_yields_both_mates_with_suffixes(tmp_path, gz):
    a = _fasta(tmp_path / ("a.fa.gz" if gz else "a.fa"), R1, gz=gz)
    b = _fasta(tmp_path / ("b.fa.gz" if gz else "b.fa"), R2, gz=gz)
    got = list(read_pairs(a, b))
    assert len(got) == 2 * len(R1)
    assert [i for i, _ in got] == ["r1/1", "r1/2", "r2/1", "r2/2",
                                   "r3/1", "r3/2", "r4/1", "r4/2"]


def test_paired_fasta_and_fastq_agree_record_for_record(tmp_path):
    """The two formats carry the same sequences, so the mapper must not be able to tell them apart
    beyond the quality string it does not use."""
    fa = list(read_pairs(_fasta(tmp_path / "a.fa", R1), _fasta(tmp_path / "b.fa", R2)))
    fq = list(read_pairs(_fastq(tmp_path / "a.fq", R1), _fastq(tmp_path / "b.fq", R2)))
    assert fa == fq


def test_a_truncated_FASTA_mate_still_raises(tmp_path):
    """The pairing check is not FASTQ-only.

    A truncated R2 makes a plain `zip` stop early and silently analyse a prefix; in this project
    that produced a published false discovery that had to be retracted. The guarantee has to hold
    whatever the input format is.
    """
    a = _fasta(tmp_path / "a.fa", R1)
    b = _fasta(tmp_path / "b.fa", R2[:2])
    with pytest.raises(ValueError):
        list(read_pairs(a, b))


# ---------------------------------------------------------------- --limit

@pytest.mark.parametrize("writer", [_fasta, _fastq])
def test_limit_is_a_head_on_single_end(tmp_path, writer):
    p = writer(tmp_path / "a.in", R1)
    assert len(list(read_pairs(p, limit=2))) == 2
    assert len(list(read_pairs(p, limit=99))) == len(R1), "a limit past the end is not an error"
    assert len(list(read_pairs(p, limit=None))) == len(R1)


@pytest.mark.parametrize("writer", [_fasta, _fastq])
def test_limit_counts_PAIRS_not_reads(tmp_path, writer):
    """`--limit 2` on a paired run means two fragments, i.e. four records out.

    Getting this backwards would silently halve or double every limited benchmark.
    """
    a = writer(tmp_path / "a.in", R1)
    b = writer(tmp_path / "b.in", R2)
    got = list(read_pairs(a, b, limit=2))
    assert len(got) == 4
    assert [i for i, _ in got] == ["r1/1", "r1/2", "r2/1", "r2/2"]


@pytest.mark.parametrize("writer", [_fasta, _fastq])
def test_limit_does_not_reach_a_divergence_beyond_it(tmp_path, writer):
    """A HEAD must not fail on something it never reads.

    This is why the limited path is pure Python rather than dnaio: dnaio validates pairing while
    filling its own buffers, so it sees — and raises on — a divergence the caller asked never to
    reach. Reading the first 2 pairs of files that disagree at record 3 must succeed.
    """
    a = writer(tmp_path / "a.in", R1)
    b = writer(tmp_path / "b.in", R2[:2] + [("SHUFFLED", "ACGT"), ("ALSO", "TTTT")])
    assert len(list(read_pairs(a, b, limit=2))) == 4
    with pytest.raises(ValueError):
        list(read_pairs(a, b))          # unlimited, so the divergence IS reached


def test_limit_zero_and_negative(tmp_path):
    """The CLI passes `limit or None`, so 0 means "no limit" there; the API takes 0 literally."""
    p = _fastq(tmp_path / "a.fq", R1)
    assert list(read_pairs(p, limit=0)) == []
