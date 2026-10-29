"""Sequence I/O: streaming FASTA/FASTQ readers and chunking.

Native parsing (no BioPython). Transparently handles gzip by ``.gz`` extension.
"""

from __future__ import annotations

import gzip
import queue
import re
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

__all__ = ["open_text", "read_sequences", "detect_format", "write_fasta", "chunked", "normalize_nt", "normalize_nt_frame"]


# Sequencing observations only: reference/template readers deliberately do not call this.
_IUPAC_TO_A = str.maketrans("RYSWKMBDHVN", "A" * 11)
_INVALID_NT = re.compile("[^ACGT]")
_INVALID_ALIGNED_NT = re.compile("[^ACGT.-]")


def normalize_nt(sequence: str, *, aligned: bool = False) -> str:
    """Uppercase sequencing DNA and replace every ambiguous DNA IUPAC base with A.

    This mandatory rule preserves positions and rejects invalid non-DNA symbols.
    ``aligned=True`` additionally permits alignment gaps (- and .), never raw
    read gaps. Never apply this sequencing-input rule to germlines.
    """
    normalized = sequence.upper().translate(_IUPAC_TO_A)
    invalid = (_INVALID_ALIGNED_NT if aligned else _INVALID_NT).search(normalized)
    if invalid:
        raise ValueError("invalid symbol in sequencing DNA input")
    return normalized


def normalize_nt_frame(frame):
    """Normalize observed AIRR nucleotide columns, preserving references and metadata.

    Re-translate existing region/junction amino acids only where their observed
    nucleotide sequence changed. Nulls, coordinates and quality strings survive.
    """
    import polars as pl
    from .._markup import translate

    regions = ("junction", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
    columns = ("sequence", "sequence_alignment", "v_sequence_alignment",
               "d_sequence_alignment", "d2_sequence_alignment", "j_sequence_alignment", "c_sequence_alignment",
               "np1", "np2", "np3", *regions)
    expressions = []
    for name in columns:
        if name not in frame.columns:
            continue
        normalized = pl.col(name).cast(pl.String).str.to_uppercase().str.replace_all("[RYSWKMBDHVN]", "A")
        invalid = "[^ACGT.-]" if name.endswith("_alignment") else "[^ACGT]"
        if frame.select(normalized.str.contains(invalid).any()).item():
            raise ValueError(f"invalid symbol in sequencing DNA column {name}")
        expressions.append(normalized.alias(name))
        aa = name + "_aa"
        if name in regions and aa in frame.columns:
            # Only changed distinct observations cross the native translation boundary.
            changed = frame.filter(pl.col(name) != normalized).select(normalized.alias(name))[name].unique().to_list()
            if changed:
                mapping = {nt: translate(nt) for nt in changed}
                expressions.append(pl.when(pl.col(name) != normalized)
                    .then(normalized.replace_strict(mapping, default=pl.col(aa)))
                    .otherwise(pl.col(aa)).alias(aa))
    return frame.with_columns(expressions)


@contextmanager
def prefetched(chunks, *, enabled=True):
    """Read at most two chunks ahead; stop the reader on consumer failure too."""
    source = iter(chunks)
    if not enabled:
        try:
            yield source
        finally:
            if hasattr(source, "close"):
                source.close()
        return
    pending = queue.Queue(maxsize=2)
    stop = threading.Event()
    errors = []
    sentinel = object()

    def put(item):
        while not stop.is_set():
            try:
                pending.put(item, timeout=0.1)
                return True
            except queue.Full:
                pass
        return False

    def read():
        try:
            for chunk in source:
                if not put(chunk):
                    break
        except BaseException as exc:
            errors.append(exc)
        finally:
            try:
                if hasattr(source, "close"):
                    source.close()
            except BaseException as exc:
                errors.append(exc)
            finally:
                put(sentinel)

    def consume():
        while (item := pending.get()) is not sentinel:
            yield item
        if errors:
            raise errors[0]

    reader = threading.Thread(target=read, name="arda-reader", daemon=True)
    reader.start()
    try:
        yield consume()
    finally:
        stop.set()
        reader.join(timeout=5)
        if reader.is_alive():
            raise RuntimeError("ARDA input reader did not stop within 5 seconds")


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
    # `zip(*[iter(fh)] * 4)` walks the file object's own buffered iterator four lines at a time
    # instead of making four `readline()` calls per record. Measured 1.43x on a 1 M-read FASTQ
    # (3.02 vs 2.11 M reads/s) -- worth having now that reading is 65 % of a prefiltered bulk run,
    # where before the prefilter it was 3 % and not worth touching.
    #
    # A trailing partial record is dropped by `zip`, the same outcome the old loop reached through
    # its `if header.startswith("@")` guard. That guard is kept: it is what stops a stray non-@
    # line from being read as a record header.
    if with_qual:
        for header, seq, _plus, qual in zip(*[iter(fh)] * 4):
            if header[:1] == "@":
                yield header[1:].split(None, 1)[0], seq.rstrip("\n"), qual.rstrip("\n")
    else:
        for header, seq, _plus, _qual in zip(*[iter(fh)] * 4):
            if header[:1] == "@":
                yield header[1:].split(None, 1)[0], seq.rstrip("\n")


def read_sequences(path: str | Path, *, with_qual: bool = False, seqtype: str = "nt") -> Iterator[tuple]:
    """Yield ``(id, sequence)`` from a FASTA or FASTQ file (auto-detected).

    Sequencing DNA IUPAC bases are always replaced with A. ``seqtype="aa"``
    denotes protein input, which is preserved.

    ``with_qual=True`` yields ``(id, sequence, qual)`` instead -- the FASTQ Phred string, or
    ``None`` for FASTA (which has no quality). The quality line is consumed either way,
    only kept when asked (needed solely by the paired
    overlap-merge, :func:`arda.rnaseq.map.merge_pair`)."""
    if seqtype not in ("nt", "aa"):
        raise ValueError("seqtype must be nt or aa")
    try:
        fmt = detect_format(path)
        with open_text(path) as fh:
            if fmt == "fasta":
                for sid, seq in _read_fasta(fh):
                    seq = normalize_nt(seq) if seqtype == "nt" else seq
                    yield (sid, seq, None) if with_qual else (sid, seq)
            else:
                for record in _read_fastq(fh, with_qual=with_qual):
                    sid, seq, *qual = record
                    yield (sid, normalize_nt(seq) if seqtype == "nt" else seq, *qual)
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
