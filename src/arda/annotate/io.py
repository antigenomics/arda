"""Sequence I/O: streaming FASTA/FASTQ readers and chunking.

Native parsing (no BioPython). Transparently handles gzip by ``.gz`` extension.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator

__all__ = ["open_text", "read_sequences", "detect_format", "write_fasta", "chunked"]


def open_text(path: str | Path):
    """Open a (possibly gzipped) text file for reading."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_format(path: str | Path) -> str:
    """Return ``"fasta"`` or ``"fastq"`` by peeking at the first non-empty char."""
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                return "fasta"
            if line[0] == "@":
                return "fastq"
            raise ValueError(f"Cannot determine format of {path}: starts with {line[0]!r}")
    raise ValueError(f"Empty input: {path}")


def _read_fasta(fh) -> Iterator[tuple[str, str]]:
    sid: str | None = None
    seq: list[str] = []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if sid is not None:
                yield sid, "".join(seq)
            sid = line[1:].split()[0]
            seq = []
        elif line:
            seq.append(line)
    if sid is not None:
        yield sid, "".join(seq)


def _read_fastq(fh, with_qual: bool = False) -> Iterator[tuple]:
    while True:
        header = fh.readline()
        if not header:
            break
        seq = fh.readline().rstrip("\n")
        fh.readline()  # '+'
        qual = fh.readline()  # quality (read either way to advance the record)
        if header.startswith("@"):
            hid = header[1:].split()[0]
            yield (hid, seq, qual.rstrip("\n")) if with_qual else (hid, seq)


def read_sequences(path: str | Path, *, with_qual: bool = False) -> Iterator[tuple]:
    """Yield ``(id, sequence)`` from a FASTA or FASTQ file (auto-detected).

    ``with_qual=True`` yields ``(id, sequence, qual)`` instead -- the FASTQ Phred string, or
    ``None`` for FASTA (which has no quality). The default is unchanged and pays nothing: the
    quality line is consumed either way, only kept when asked (needed solely by the paired
    overlap-merge, :func:`arda.rnaseq.map.merge_pair`)."""
    try:
        fmt = detect_format(path)
        with open_text(path) as fh:
            if fmt == "fasta":
                for sid, seq in _read_fasta(fh):
                    yield (sid, seq, None) if with_qual else (sid, seq)
            else:
                yield from _read_fastq(fh, with_qual=with_qual)
    except (EOFError, gzip.BadGzipFile) as exc:
        # A truncated/corrupt .gz raises deep in the gzip layer mid-stream (EOFError:
        # "Compressed file ended before the end-of-stream marker"). Surface it as a clear
        # input error: a bare EOFError reaching Typer/Click is mistaken for a Ctrl-D and
        # printed as "Aborted." with no cause.
        raise ValueError(f"truncated or corrupt gzip input: {path}") from exc


def write_fasta(records: Iterator[tuple[str, str]], path: str | Path) -> Path:
    """Write ``(id, sequence)`` records to a FASTA file."""
    path = Path(path)
    with open(path, "w") as fh:
        for sid, seq in records:
            fh.write(f">{sid}\n{seq}\n")
    return path


def chunked(it: Iterator, size: int) -> Iterator[list]:
    """Yield lists of up to ``size`` items from an iterator."""
    batch: list = []
    for item in it:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch
