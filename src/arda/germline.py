"""Per-allele V and J germline nucleotide sequences, derived from the COMMITTED reference.

arda's runtime reference is a V·J **scaffold** cross product: ``alleles.fasta`` holds assembled
``V + N·pad + J`` sequences under positional ids (``TRA_0``), and ``markup.tsv`` says which
``v_call``/``j_call`` each one carries and where the segments start and end. So the per-allele
germlines are already in there -- ``scaffold[:v_sequence_end]`` for V, ``scaffold[j_sequence_start
- 1:vj_end]`` for J -- and this module is the ten lines that take them back out.

Why not :func:`arda.refbuild.imgt.load_functional_alleles`
---------------------------------------------------------
That is where :func:`arda.annotate.ties.resolve_airr` used to read them, and it is wrong on three
counts:

**Never: it reads the IMGT SOURCE directory, which a ``pip install`` does not have.** The source
tree lives under ``data/imgt`` and is populated by ``arda build-db``; the shipped reference tarball
carries neither. So ``arda resolve-ties`` raised on every plain install -- the same shape as the
``segments.fasta`` deploy trap, which also looked fine in a source checkout.

**Never: the coordinate spaces differ.** ``build_locus_scaffolds`` trims each V to its coding frame
(``combinations.py``: ``s[detect_coding_frame(s):]``) before assembling, so the scaffold's V -- and
therefore every ``v_germline_start``/``v_germline_end`` arda emits -- is shifted 1 or 2 nt against
the raw IMGT sequence for **10 of 884** human V alleles. Slicing a raw germline with scaffold
coordinates silently tests a window that is off by up to two bases.

**Never: it returns alleles the reference cannot call.** 884 human V alleles are functional in
IMGT; 775 survive ``drop_unanchorable`` / ``drop_truncated`` and complete markup into the built
reference. Expanding a tie list into one of the other 109 offers an allele that could never have
been the answer.

All three go away by deriving from the artifact that is actually shipped and actually searched.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache

import polars as pl

from ._log import logger
from .paths import vdj_dir
from .refbuild.imgt import read_fasta

__all__ = ["SEGMENTS", "segment_germlines"]

#: Segments this module can recover. D is absent on purpose: D germlines are not part of a
#: scaffold (the junction interior is N-padded), and they already ship whole as
#: ``d_germlines.fasta``.
SEGMENTS = ("v", "j")

#: ``segment -> (call column, start column or None for "the scaffold start", end column)``.
#: Coordinates in ``markup.tsv`` are 1-based closed, as everywhere else in arda.
_SPAN = {
    "v": ("v_call", None, "v_sequence_end"),
    "j": ("j_call", "j_sequence_start", "vj_end"),
}


def _int(v) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


@lru_cache(maxsize=16)
def segment_germlines(organism: str, segment: str = "v") -> dict[str, str]:
    """``{allele: germline_nt}`` for one segment, across every locus of ``organism``.

    The sequences are in **scaffold space** -- V trimmed to its coding frame -- which is the space
    arda's ``*_germline_start``/``*_germline_end`` columns are already in, so a span may be sliced
    out of these directly.

    Raises ``FileNotFoundError`` when the reference is not present.
    """
    if segment not in _SPAN:
        raise ValueError(f"unknown segment {segment!r}; one of {', '.join(SEGMENTS)}")
    base = vdj_dir(organism)
    markup, fasta = base / "markup.tsv", base / "alleles.fasta"
    if not markup.exists() or not fasta.exists():
        raise FileNotFoundError(
            f"no reference for {organism!r} at {base} (looked for markup.tsv + alleles.fasta)")

    seqs = dict(read_fasta(fasta))
    call_col, start_col, end_col = _SPAN[segment]
    # Same reader as `annotate.reference._load_markup`: default quoting, so a blank cell written
    # as the two characters `""` parses back to the empty string rather than staying literal.
    df = pl.read_csv(markup, separator="\t", infer_schema_length=0)
    if call_col not in df.columns or end_col not in df.columns:
        raise FileNotFoundError(f"{markup} has no {call_col}/{end_col} column")

    observed: dict[str, Counter[str]] = {}
    for row in df.iter_rows(named=True):
        allele = (row.get(call_col) or "").strip()
        end = _int(row.get(end_col))
        scaffold = seqs.get(row["scaffold_id"])
        if not allele or end <= 0 or scaffold is None:
            continue
        start = _int(row.get(start_col)) - 1 if start_col else 0
        if start < 0 or end <= start:
            continue
        germline = scaffold[start:end]
        # Never: `markup.tsv`'s `v_call` is NOT always one allele. 23 of the human reference's 775
        # V entries are comma-joined groups of alleles whose scaffolds came out byte-identical and
        # were deduplicated (`IGKV1-33*01,IGKV1D-33*01`, ...), hiding 49 alleles. Keying on the
        # group STRING silently breaks every consumer: `TieResolver.expand` takes
        # `call.split(",")[0]` and looks that up, misses the group key, and returns the call
        # untouched -- a tie list that never fires, reported as success. Split to 801 real keys.
        for member in (a.strip() for a in allele.split(",")):
            if member:
                observed.setdefault(member, Counter())[germline] += 1

    # Never: an allele sits on every scaffold that pairs it with a J, and those copies are NOT
    # guaranteed identical -- `v_sequence_end` is IgBLAST's markup of one assembled scaffold, not
    # a property of the allele, and it wobbles. Measured across all five shipped organisms: human,
    # rat, rabbit and rhesus agree everywhere; mouse disagrees on 4 of 897 V and 10 of 109 J
    # alleles, always by 1-2 nt at one edge and always with a landslide majority (`TRAV16*02` is
    # 288 nt on 58 scaffolds and 290 nt on 1; `TRAJ42*01` 64 nt on 218 and 63 nt on 1).
    # So: the MODE, not the first row -- last-wins would have handed the whole reference's tie
    # resolution to whichever scaffold happened to sort last. Ties break by longer, then by
    # sequence, which is a TOTAL order: two runs over one reference cannot disagree.
    out = {a: min(seen, key=lambda g: (-seen[g], -len(g), g)) for a, seen in observed.items()}
    wobbled = sum(1 for seen in observed.values() if len(seen) > 1)
    if wobbled:
        logger.debug("%s %s: %d/%d alleles had inconsistent markup across scaffolds; took the mode",
                     organism, segment.upper(), wobbled, len(observed))
    return out
