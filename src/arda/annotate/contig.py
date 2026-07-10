"""Two ways to give an assembled contig its AIRR cigars — and they agree.

A Stage-3 contig is a consensus of reads that Stage 1 already aligned to a scaffold.
Its ``v_cigar``/``j_cigar``/``c_cigar`` + alignment strings can be produced two ways:

* :func:`reannotate_contigs` — treat the contig as one long query and run it back
  through :func:`~arda.annotate.mapper.annotate_records` (one mmseqs alignment, then
  ``segment_cigars``). Simple, exact, no new code; the cost is a second alignment pass.
* :func:`merge_contigs` — stitch the reads' existing alignments into the contig's
  (the C++ ``_markup.merge_alignment`` per-column consensus over N reads), skipping the
  alignment pass. Wins when a sample has ~10^5 contigs (scRNA-seq).

Both converge on the same synthetic ``hit`` and reuse :func:`~arda.annotate.transfer.transfer_hit`,
so their output is field-for-field comparable. Which is optimal is a measured question —
see ``tests/unit/test_contig_merge.py`` and the arda-benchmark Phase-D benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .. import _markup
from .mapper import annotate_records
from .reference import Reference, load_reference
from .transfer import transfer_hit

__all__ = ["ReadPlacement", "Contig", "reannotate_contigs", "merge_contig", "merge_contigs"]


@dataclass(slots=True)
class ReadPlacement:
    """One read's placement in a contig: its scaffold alignment + contig offset.

    ``qaln``/``taln`` are the read's coding-strand aligned strings vs the scaffold
    (``-`` for gaps, as Stage 1 emits). ``qstart``/``tstart`` are 1-based starts in the
    read / scaffold. ``offset`` is the 0-based position of the read's first base within
    the contig (the assembly layout), in contig orientation.
    """

    qaln: str
    taln: str
    qstart: int
    tstart: int
    offset: int


@dataclass(slots=True)
class Contig:
    """An assembled contig + the reads it was built from, all hitting one scaffold."""

    sequence_id: str
    sequence: str
    target: str                                   # scaffold id the reads agree on
    reads: list[ReadPlacement] = field(default_factory=list)


def reannotate_contigs(
    records: list[tuple[str, str]],
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    sensitivity: float | None = None,
    strand: str = "both",
    map_d: bool = True,
) -> list[dict]:
    """Annotate assembled contigs by re-aligning them (baseline path).

    ``records`` are ``(contig_id, contig_seq)``. A thin wrapper over
    :func:`~arda.annotate.mapper.annotate_records`: a contig is just a long query.
    """
    return annotate_records(records, organism, seqtype, threads=threads,
                            sensitivity=sensitivity, strand=strand, map_d=map_d)


def merge_contig(contig: Contig, reference: Reference, *, map_d: bool = True) -> dict:
    """Annotate one contig by stitching its reads' alignments (merge path).

    ``reference`` is a preloaded :class:`~arda.annotate.reference.Reference` (load it
    once for a whole sample; see :func:`merge_contigs`). Raises ``KeyError`` if the
    contig's ``target`` scaffold is absent from the reference.
    """
    entry = reference.get(contig.target)
    if entry is None:
        raise KeyError(f"scaffold {contig.target!r} not in reference {reference.organism!r}")
    reads = contig.reads
    qaln, taln, qstart, qend, tstart, tend = _markup.merge_alignment(
        [r.qaln for r in reads], [r.taln for r in reads],
        [r.qstart for r in reads], [r.tstart for r in reads],
        [r.offset for r in reads], contig.sequence)
    # A synthetic hit identical in shape to a fresh mmseqs hit, so transfer_hit does the
    # rest (cigars, region markup, germline coords). Alignment-quality scalars are blank:
    # a merged alignment has no single bit score / E-value.
    hit = {"qaln": qaln, "taln": taln, "qstart": qstart, "qend": qend,
           "tstart": tstart, "tend": tend, "qlen": len(contig.sequence),
           "tlen": "", "bits": "", "evalue": "", "pident": "", "target": contig.target}
    dg = reference.d_germlines.get(entry.locus) if map_d else None
    return transfer_hit(contig.sequence_id, contig.sequence, hit, entry, reference.seqtype,
                        d_germlines=dg, submitted_seq=contig.sequence,
                        anchors=reference.anchors)


def merge_contigs(
    contigs: list[Contig],
    organism: str = "human",
    seqtype: str = "nt",
    *,
    reference: Reference | None = None,
    map_d: bool = True,
) -> list[dict]:
    """Annotate contigs by the merge path; loads the reference once for all of them."""
    ref = reference or load_reference(organism, seqtype)
    return [merge_contig(c, ref, map_d=map_d) for c in contigs]
