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


def _read_fastq(fh) -> Iterator[tuple[str, str]]:
    while True:
        header = fh.readline()
        if not header:
            break
        seq = fh.readline().rstrip("\n")
        fh.readline()  # '+'
        fh.readline()  # quality
        if header.startswith("@"):
            yield header[1:].split()[0], seq


def read_sequences(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(id, sequence)`` from a FASTA or FASTQ file (auto-detected)."""
    fmt = detect_format(path)
    with open_text(path) as fh:
        if fmt == "fasta":
            yield from _read_fasta(fh)
        else:
            yield from _read_fastq(fh)


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
