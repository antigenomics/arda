"""Lift a cell barcode out of a sequence identifier that already carries one.

arda has no barcode or UMI concept and does not need one: the tools that produce single-cell input
put the barcode in the record NAME, and the name is what arda reads -- ``map`` yields
``rec.id``, the first whitespace-delimited token. (``dnaio`` itself keeps the FASTQ comment on
``.comment``; arda is what drops it. See ``project/design-singlecell.md``.) This module is the
whole of the single-cell front end -- a pure parser, no I/O.

Three dialects ship, and a regex escape hatch for the rest:

===============  ==============================================================
``cellranger``   ``AAACCTGAGAAACCAT-1_contig_2`` -- 10x, with a gem-group suffix
``migec``        ``PBMC.AAACCTGCAGCCTGTT.CGTTTTTATC`` and its ``.c<k>``/``.<m>``
                 contig and split variants
``prefix``       ``<barcode>_<anything>`` / ``<barcode>.<anything>``
===============  ==============================================================

Named dialects rather than regex alone, because each carries a validator a regex cannot express --
see ``parse_migec`` below, which is a right-to-left walk and not a pattern at all.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

__all__ = ["CellKey", "DIALECTS", "parse", "sniff", "make_parser"]

_DNA = re.compile(r"^[ACGTNacgtn]+$")
_CELLRANGER = re.compile(r"^(?P<cell>[ACGTNacgtn]+)(?:-\d+)?_contig_\d+$")
_PREFIX = re.compile(r"^(?P<cell>[ACGTNacgtn]+)[_.]")

DIALECTS = ("cellranger", "migec", "prefix")

#: What :func:`sniff` will consider, and it is NOT all of :data:`DIALECTS`.
#:
#: Never: ``prefix`` is opt-in only. Its rule is "leading run of ACGTN before a separator", which any
#: bulk sample id drawn from those four letters satisfies -- ``TCGA.CGTTTTTATC`` parses to a cell
#: ``TCGA``, consistently, across every record in the file, so neither the match rate nor the
#: length check catches it. Naming a dialect is a claim the caller is entitled to make; guessing
#: this one is a claim arda is not.
_SNIFFABLE = ("cellranger", "migec")

#: Shortest barcode :func:`sniff` will believe. The narrowest in use is 10 nt (Drop-seq's cell
#: barcode); 10x is 16. Below this a "barcode" is far likelier to be a sample id that looks like
#: DNA than a real one.
_MIN_CELL_LEN = 10


@dataclass(frozen=True, slots=True)
class CellKey:
    """What a single-cell identifier resolves to.

    ``contig`` numbers the fragments of ONE molecule; ``molecule`` numbers DIFFERENT molecules that
    shared a barcode. Both default to 1, which is what an identifier carrying neither means.
    """

    cell: str | None = None
    umi: str | None = None
    sample: str | None = None
    contig: int = 1
    molecule: int = 1


def parse_cellranger(sequence_id: str) -> CellKey | None:
    m = _CELLRANGER.match(sequence_id)
    return CellKey(cell=m.group("cell").upper()) if m else None


def parse_prefix(sequence_id: str) -> CellKey | None:
    m = _PREFIX.match(sequence_id)
    return CellKey(cell=m.group("cell").upper()) if m else None


def parse_migec(sequence_id: str) -> CellKey | None:
    """``<sample>[.<cell>].<umi>[.c<k>][.<m>]``, parsed RIGHT TO LEFT.

    Never: **not a left-to-right split.** migec's ``validate_sample_id`` forbids ``/``, ``\\``,
    ``..``, a leading ``-`` and control characters -- it does NOT forbid a ``.``, so ``PBMC.rep2``
    is a legal sample id and splitting from the left reads ``rep2`` as the barcode.

    Never: **the optional suffixes are emitted conditionally**, so the field count is 2, 3, 4 or 5
    and a fixed index is wrong for some of them (``migec/src/assemble.cpp:503-505``).
    """
    fields = sequence_id.split(".")
    if len(fields) < 2:
        return None
    molecule = contig = 1
    if len(fields) > 2 and fields[-1].isdigit():
        molecule = int(fields.pop())
    if len(fields) > 2 and len(fields[-1]) > 1 and fields[-1][0] in "cC" \
            and fields[-1][1:].isdigit():
        contig = int(fields.pop()[1:])
    if not _DNA.match(fields[-1]):
        return None
    umi = fields.pop().upper()
    cell = fields.pop().upper() if len(fields) > 1 and _DNA.match(fields[-1]) else None
    return CellKey(cell=cell, umi=umi, sample=".".join(fields) or None,
                   contig=contig, molecule=molecule)


_PARSERS = {"cellranger": parse_cellranger, "migec": parse_migec, "prefix": parse_prefix}


def parse(sequence_id: str, dialect: str) -> CellKey | None:
    try:
        fn = _PARSERS[dialect]
    except KeyError:
        raise ValueError(f"unknown cell-id dialect {dialect!r}; pick from {', '.join(DIALECTS)} "
                         f"or pass a regex") from None
    return fn(sequence_id)


def sniff(sequence_ids, *, min_match: float = 0.95) -> str | None:
    """The dialect that parses at least ``min_match`` of ``sequence_ids`` to a consistent cell.

    Never: **an alphabet check alone fabricates cells.** A bulk sample named ``TCGA``, ``CAG`` or
    any donor code over those four letters parses as a one-cell library, and nothing downstream
    flags it. So a dialect is only accepted when every parsed cell has the SAME LENGTH -- a real
    10x cell field is exactly 16 nt in every record, and a sample id that happens to look like DNA
    is not.

    The caller is responsible for handing this a REPRESENTATIVE sample. migec writes in barcode
    order, so the first N identifiers of a file share their leading bases and are not one.
    """
    ids = list(sequence_ids)
    if not ids:
        return None
    for dialect in _SNIFFABLE:
        keys = [parse(i, dialect) for i in ids]
        cells = [k.cell for k in keys if k is not None and k.cell]
        if len(cells) < min_match * len(ids):
            continue
        lengths = {len(c) for c in cells}
        if len(lengths) != 1 or min(lengths) < _MIN_CELL_LEN:
            continue
        return dialect
    return None


def make_parser(dialect: str = "auto", regex: str | None = None, *, sample=(),
                with_sample: bool = False):
    """A ``sequence_id -> cell or None`` callable.

    ``regex`` wins when given and must define a ``cell`` group. ``auto`` runs :func:`sniff` over
    ``sample`` and refuses rather than guessing when nothing fits.

    ``with_sample`` returns ``(sample, cell)`` instead, so a caller can check that one file holds
    one sample. A 10x barcode comes from a fixed 737,280-entry whitelist, so the same ``cell``
    across two samples is guaranteed, not unlucky -- keying on the barcode alone merges two cells
    and nothing downstream can tell. Only the ``migec`` dialect carries a sample; the others
    return ``None`` for it.
    """
    if regex is not None:
        rx = re.compile(regex)
        if "cell" not in (rx.groupindex or {}):
            raise ValueError(f"--cell-regex {regex!r} has no (?P<cell>...) group")
        if with_sample:
            return lambda sid: (None, m.group("cell").upper() if (m := rx.match(sid)) else None)
        return lambda sid: (m.group("cell").upper() if (m := rx.match(sid)) else None)
    if dialect == "auto":
        found = sniff(sample)
        if found is None:
            raise ValueError(
                "--cell-from auto: no dialect parses a consistent cell barcode out of these "
                f"identifiers (tried {', '.join(DIALECTS)}). Name the dialect, or pass "
                "--cell-regex")
        dialect = found
    if dialect not in _PARSERS:
        raise ValueError(f"unknown cell-id dialect {dialect!r}; pick from {', '.join(DIALECTS)}")
    fn = _PARSERS[dialect]
    if with_sample:
        return lambda sid: ((k.sample, k.cell) if (k := fn(sid)) else (None, None))
    return lambda sid: (k.cell if (k := fn(sid)) else None)


def modal_length(sequence_ids, dialect: str) -> int | None:
    """The most common parsed cell length, for a caller that wants to report the disagreement."""
    lengths = Counter(len(k.cell) for i in sequence_ids
                      if (k := parse(i, dialect)) is not None and k.cell)
    return lengths.most_common(1)[0][0] if lengths else None
