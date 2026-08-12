"""Stage 1 — RNA-seq filter + map.

Reuses the streaming, memory-bounded annotator (``annotate.mapper._prep`` +
``_annotate_chunk`` + the background-reader/bounded-queue loop of ``annotate_file``):
MMseqs2 is the parallel layer and its k-mer prefilter rejects non-receptor reads
before alignment, so mostly-non-receptor RNA-seq is cheap. The difference from
``arda annotate`` is that we write **only the reads that map** (keyed by read id, so
the AIRR TSV *is* the read-id → junction map), plus an optional candidate FASTA and a
run report.

Paired FASTQ mates are streamed independently, tagged ``<id>/1`` / ``<id>/2`` so query
ids stay unique; a pair is kept if either mate maps (recall-first — the base id
recovers the pair).
"""

from __future__ import annotations

import gzip
import json
import logging
import queue
import threading
from itertools import islice, zip_longest
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

try:                      # C-accelerated FASTQ/FASTA parsing; see _read_pairs_dnaio
    import dnaio as _dnaio
except ImportError:       # pure-Python fallback below keeps a source checkout working
    _dnaio = None

from ._res import Stage, peak_rss_mb
from .._log import Throttle
from ..annotate import io as seqio
from ..annotate import mapper
from ..annotate.airr_out import airr_header, format_rows
from ..refbuild.translate import reverse_complement
from ..shm import FULL_COLUMNS, SHM_MODES

logger = logging.getLogger(__name__)

__all__ = ["map_rnaseq", "read_pairs", "merge_pair", "junction_quality", "mutation_quality",
           "RnaseqReport"]

_MERGE_ANCHOR = 12  # exact k-mer used to locate the R1/rc(R2) overlap in O(len)

#: Non-schema AIRR column written by ``--junction-quality``: the read's Phred+33 string over
#: exactly the bases of ``junction``, in the same orientation, so position *i* of one indexes
#: position *i* of the other. Stage 2 reads it for ``correct --min-junction-q``.
JUNCTION_QUALITY = "junction_quality"

#: Non-schema AIRR columns written by ``--mutation-quality``: the Phred score of the READ BASE
#: behind each entry of ``v_mutations`` / ``j_mutations``, comma-joined, one-for-one and in the
#: same order. This is what separates a NOVEL ALLELE from somatic hypermutation from a miscall:
#: all three look identical in the mutation list, and only a recurrent, high-quality mutation is
#: evidence of a germline the reference does not carry. ``arda stats`` reads them.
MUTATION_QUALITY = ("v_mutation_quality", "j_mutation_quality")


def mutation_quality(rec: dict, qual: str) -> dict[str, str]:
    """Per-mutation Phred scores for ``rec["v_mutations"]`` / ``rec["j_mutations"]``.

    ⛔ **Driven by the mutation list that was EMITTED, not by re-deriving one.** Walking the
    alignment and scoring every mismatch reproduces what ``_markup.segment_cigars`` found -- which
    since 2.16.0 is a SUPERSET of what the columns carry, because ``arda.shm`` then drops the
    junction-internal entries (measured on this repo's own fixture: 25 of 242 V rows had more
    mismatches than mutations). The result lines up in length only by accident and pairs entry *i*
    with a different base's score. So the walk builds *germline position -> query position* and
    each emitted entry looks its own position up; an entry whose position the alignment does not
    cover yields ``""`` for the whole segment rather than a short, misaligned list.

    ⛔ Quality is oriented like the READ AS SUBMITTED; the alignment and every coordinate here are
    on the coding strand. Same reversal rule as :func:`junction_quality`.

    ``ponytail:`` a Python pass over the alignment, not a second output from the C++ walk that
    already visits these columns. It rides its own flag, so it never lands on a mode run; move it
    into ``segment_cigars`` if it ever shows up in a profile.
    """
    out = {c: "" for c in MUTATION_QUALITY}
    qaln, taln = rec.get("sequence_alignment") or "", rec.get("germline_alignment") or ""
    lists = [str(rec.get(c) or "") for c in ("v_mutations", "j_mutations")]
    if not qaln or not taln or not qual or not any(lists):
        return out
    seq = rec.get("sequence") or ""
    if str(rec.get("rev_comp") or "").upper() in ("T", "TRUE", "1"):
        qual, seq = qual[::-1], reverse_complement(seq)
    if len(qual) != len(seq):
        return out
    try:
        # `_num` writes these as floats, a bare producer as ints; both parse through float().
        q = int(float(rec["mmseqs2_qstart"]))
        t = int(float(rec["mmseqs2_tstart"]))
        t_vend = int(float(rec.get("mmseqs2_t_vend") or 0))
        t_jstart = int(float(rec.get("mmseqs2_t_jstart") or 0))
        t_vjend = int(float(rec.get("mmseqs2_t_vjend") or 0))
    except (KeyError, TypeError, ValueError):
        return out

    # germline position -> 1-based query position, for V (0) and J (1). Same segment test as
    # `seg_key`: the N-pad between them and the C region carry no SHM evidence.
    where: tuple[dict[int, int], dict[int, int]] = ({}, {})
    for qc, tc in zip(qaln, taln):
        cq, ct = qc != "-", tc != "-"
        if cq and ct:
            if t_vend and t <= t_vend:
                where[0][t] = q
            elif t_jstart and t_vjend and t_jstart <= t <= t_vjend:
                where[1][t - t_jstart + 1] = q
        q += cq
        t += ct

    for col, entries, table in zip(MUTATION_QUALITY, lists, where):
        if not entries:
            continue
        scores: list[str] = []
        for entry in entries.split(","):
            try:                              # `G191C` -> germline position 191
                pos = table[int(entry[1:-1])]
            except (KeyError, ValueError):
                scores = []
                break
            if not 1 <= pos <= len(qual):
                scores = []
                break
            scores.append(str(ord(qual[pos - 1]) - 33))
        out[col] = ",".join(scores)
    return out


def junction_quality(rec: dict, qual: str) -> str:
    """The read's Phred+33 substring covering ``rec["junction"]``, or ``""``.

    ⛔ The quality string belongs to the read AS SUBMITTED, while the junction and every coordinate
    in the record are on the CODING strand. For ``rev_comp == "T"`` the two run in opposite
    directions, so the quality is reversed (not complemented -- a Phred char has no complement)
    before slicing. Getting that backwards yields a string of the right LENGTH holding the wrong
    bases' qualities, which no length or format check downstream can catch. So the slice is
    verified against the junction it claims to describe before it is returned.

    The junction is ``coding[cdr3_start - 3 : cdr3_end + 3]`` -- the CDR3 flanked by the Cys104 and
    [FW]118 anchor codons, see :func:`arda.annotate.transfer._junction_nt` -- so the coordinates
    place it in O(1). They can be absent (a producer that did not fill the region columns), hence
    the ``find`` fallback; if neither reproduces the junction, return ``""`` rather than a
    misaligned string.
    """
    jn = rec.get("junction") or ""
    if not jn or not qual:
        return ""
    seq = rec.get("sequence") or ""
    if str(rec.get("rev_comp") or "").upper() in ("T", "TRUE", "1"):
        qual, seq = qual[::-1], reverse_complement(seq)
    if len(qual) != len(seq):
        return ""
    cs, ce = rec.get("cdr3_start"), rec.get("cdr3_end")
    if cs and ce:
        s = int(cs) - 4                                   # 0-based start of the Cys104 codon
        if s >= 0 and seq[s:s + len(jn)] == jn:
            return qual[s:s + len(jn)]
    p = seq.find(jn)
    return qual[p:p + len(jn)] if p >= 0 else ""


def _constant_only(rec: dict) -> bool:
    """Does this alignment lie WHOLLY inside the constant region?

    `mmseqs2_t_vjend` is the end of the V-J part of the winning scaffold: the scaffold length for a
    V-J scaffold (so this is never true), the J length for a `J + C` scaffold. An alignment starting
    at or after it never touched a J, so the read is receptor mRNA carrying no V(D)J -- a real
    transcript, but no clonotype.
    """
    vjend = rec.get("mmseqs2_t_vjend")
    tstart = rec.get("mmseqs2_tstart")
    if not vjend or not tstart:
        return False
    return float(tstart) >= float(vjend)


def _fragment_id(sequence_id: str) -> str:
    return sequence_id[:-2] if sequence_id.endswith(("/1", "/2")) else sequence_id


def _apply_constant_rule(records: list[dict]) -> tuple[list[dict], int, int]:
    """Drop constant-only *fragments*; keep a constant-only mate's isotype for its V(D)J partner.

    This has to be decided per FRAGMENT, not per read. Insert size exceeds 2x read length for a large
    share of libraries -- 36.4 % of pairs in PRJNA371303, median insert 145 nt against 2x100 bp reads
    -- so the mates do not overlap and each samples a different part of the transcript. The commonest
    informative layout is exactly: **R1 across V/CDR3, R2 deep in the constant region with no J.**

    Judged read by read, that R2 is "constant-only" and gets discarded -- along with the only isotype
    evidence the fragment carries. So: a fragment survives if ANY of its reads touched V(D)J, and a
    constant-only mate then donates its `c_call`/`c_class` to the fragment instead of being thrown
    away. A fragment whose every read lies inside C is receptor mRNA with no rearrangement: dropped.

    Mates are adjacent in the stream (`read_pairs` interleaves `/1`,`/2`), so they share a chunk except
    for at most one pair per chunk boundary; that pair simply loses the isotype donation.

    Returns:
        (kept records, fragments dropped as constant-only, isotypes donated by a mate).
    """
    by_frag: dict[str, list[dict]] = {}
    for r in records:
        by_frag.setdefault(_fragment_id(r["sequence_id"]), []).append(r)

    drop_ids: set[int] = set()
    dropped_frags = donated = 0
    for recs in by_frag.values():
        vdj = [r for r in recs if not _constant_only(r)]
        if not vdj:
            dropped_frags += 1
            drop_ids.update(id(r) for r in recs)
            continue
        donor = next((r for r in recs if _constant_only(r) and r.get("c_call")), None)
        if donor is not None:
            for r in vdj:
                if not r.get("c_call"):
                    r["c_call"], r["c_class"] = donor["c_call"], donor["c_class"]
            donated += 1
        drop_ids.update(id(r) for r in recs if _constant_only(r))

    return [r for r in records if id(r) not in drop_ids], dropped_frags, donated

# Minimum MMseqs2 bit score for a mapped read to be reported.
#
# The knob stays; what changed is that its exact value no longer matters. Once the reference carries
# `J + C` scaffolds, a read crossing the J->C splice has somewhere legal to end, and the recall/score
# curve goes flat (arda-benchmark OPTIMIZATION.md §6.2, §6.5):
#
#     min_score        0      40      55      75      85
#     V+pad+J only   .9556   .9550   .9273   .8761   .8422     <- 8 points of recall live in the knob
#     + J+C, + P1    1.000   1.000   1.000   .9994   .9803     <- flat; the knob does nothing
#
# and precision rises monotonically with it (.9327 -> .9586 across 40..75) at no recall cost. 75 sits
# in the flat interior. The old calibration comment claimed a "cliff at 78, collapse to 0.74 at 80":
# that cliff was the missing constant region, not a property of the score. Re-derive with
# `scripts/reference_variants.py` if the reference changes again.
_MIN_SCORE = 75.0


def _overlap_consensus(a: str, qa: str, b: str, qb: str) -> str:
    """Per-base consensus of two aligned overlap strings: keep the higher-Phred base on a mismatch
    (R1 wins a tie). Phred+33 is monotonic in the quality char, so compare the chars directly."""
    return "".join(a[i] if (a[i] == b[i] or qa[i] >= qb[i]) else b[i] for i in range(len(a)))


def merge_pair(s1: str, s2: str, *, q1: str | None = None, q2: str | None = None,
               min_overlap: int = 12, max_mismatch_rate: float = 0.1) -> str | None:
    """Overlap-merge a read pair into one fragment, or ``None`` if they don't overlap.

    Aligns ``s1`` (R1) with ``reverse_complement(s2)`` (R2 flipped to the same strand)
    by finding an exact ``_MERGE_ANCHOR``-mer from the flipped R2's 5' end inside R1
    (C-level ``str.find`` → O(len), so non-overlapping pairs — the common RNA-seq case —
    cost almost nothing), then verifying the implied overlap; the mate provides the V/J
    context a short read lacks.

    In the overlap the two mates may disagree. Given Phred qualities ``q1``/``q2`` the base
    with the higher quality wins per position (``rc(R2)``'s quality is ``q2`` reversed, not
    complemented); without them R2 wins the whole overlap (the historical behaviour). Outside
    the overlap, R1 supplies its 5' part and R2 its 3' tail.
    """
    r2 = reverse_complement(s2)
    if len(s1) < min_overlap or len(r2) < min_overlap:
        return None
    pos = s1.find(r2[:_MERGE_ANCHOR])
    if pos < 0:
        return None
    ov = min(len(s1) - pos, len(r2))
    if ov < min_overlap:
        return None
    a, b = s1[pos:pos + ov], r2[:ov]
    mism = sum(1 for x, y in zip(a, b) if x != y)
    if mism > max_mismatch_rate * ov:
        return None
    if q1 is None or q2 is None:
        overlap = b                                   # no quality: R2 wins the overlap (legacy)
    else:
        overlap = _overlap_consensus(a, q1[pos:pos + ov], b, q2[::-1][:ov])
    return s1[:pos] + overlap + r2[ov:]

# Records per mmseqs invocation. Larger chunks amortise the fixed per-call cost; memory stays
# bounded because only one chunk is resident.
#
# 400k, not 200k, and not "as large as possible" -- both ends were measured on real data
# (10 full-depth cluster runs, 754.7M reads, plus a chunk-size sweep on 1.2M records, 8 threads):
#
#   chunks:    12        6        3        1
#   wall:   51.3s    47.6s    46.3s    47.7s      <- minimum at 3 chunks (~400k records)
#
# Going UP from 400k is slower, not faster. The reader is a daemon thread behind a
# `queue.Queue(maxsize=2)`, so with one giant chunk there is nothing left to overlap and the
# +3.0% regression is the pipeline serialising. The measured per-invocation intercept is
# **0.56 s** -- not the ~1.6 s an earlier 200k-read profile implied -- which is why this knob is
# worth only ~2.6% and why "batch every mmseqs call into one" is NOT the big lever it looked
# like. The real cost model is:
#
#     wall_map ~= total_reads / 44,470  +  mapped_reads / 681
#
# i.e. a read that HITS costs ~65x one that does not, so on anything but a cold library the
# alignment term dominates and chunking cannot touch it. See `scripts/bench_cost_model.py`.
#
# Safe to change because `chunked_fragments` made the output invariant to it -- verified:
# identical AIRR checksum and identical `isotype_from_mate` at 50k / 200k / 400k.
_RNASEQ_CHUNK = 400_000


def frag_stem(i: str) -> str:
    """Fragment id for a read id: the mate suffix stripped.

    Tolerates the ``/1``,``/2`` and `` 1:N:0:`` conventions already present in the wild.
    """
    i = i.split()[0]
    return i[:-2] if i.endswith(("/1", "/2")) else i


def chunked_fragments(records: Iterator[tuple], size: int) -> Iterator[list]:
    """Chunk records without ever splitting a FRAGMENT across two chunks.

    ``_apply_constant_rule`` decides per fragment *within a chunk*, so a boundary falling
    between two mates costs that fragment its isotype donation (`isotype_from_mate`). Plain
    ``chunked`` gets away with it by accident: in the default paired path records arrive
    strictly two-per-fragment, so a boundary at a multiple of the chunk size never lands
    mid-fragment. ``--reconstruct`` destroys that accident -- a merged pair emits one record
    and an unmerged pair two, so the parity drifts and boundaries do land mid-fragment.

    Without this, `map` output depends on ``--chunk-size``, and under sharding it would
    depend on the shard layout too -- which is exactly the byte-identity the SLURM and
    Nextflow paths are supposed to guarantee.

    A chunk may exceed ``size`` by at most one fragment.

    ``frag_stem`` runs ONLY at a possible cut point, i.e. once the batch is already full. It used
    to run per record -- 1.2 M calls on a 1.2 M-read run, each doing a ``split()`` that allocates a
    list, to make about ten decisions. Profiled at 1.35 s, which is real now that the prefilter has
    removed the search it used to hide behind. The cut decisions are unchanged: the previous
    record's stem is read back off ``batch[-1]`` at the moment it is actually needed.
    """
    batch: list = []
    for rec in records:
        if len(batch) >= size and frag_stem(rec[0]) != frag_stem(batch[-1][0]):
            yield batch
            batch = []
        batch.append(rec)
    if batch:
        yield batch


def read_pairs(r1: str | Path, r2: str | Path | None = None,
               *, reconstruct: bool = False, limit: int | None = None,
               with_qual: bool = False) -> Iterator[tuple]:
    """Stream ``(id, sequence)`` reads for single-end (``r1`` only) or paired input.

    ``with_qual`` yields ``(id, sequence, quality)`` instead — the Phred+33 string, ``""`` for
    FASTA input. It is incompatible with ``reconstruct``: a merged fragment's bases come from two
    different reads, so no single input quality string describes it (:func:`map_rnaseq` refuses
    the combination rather than emitting a quality that does not belong to the sequence).

    For paired input the two mates carry the same id, so they are tagged ``<id>/1`` and
    ``<id>/2`` to keep query ids unique (strip the suffix to recover the pair). With
    ``reconstruct``, overlapping mates are merged into one fragment (:func:`merge_pair`)
    keyed by the bare id — giving a short read the mate's V/J context; non-overlapping
    mates fall back to the tagged-independent form.

    Args:
        limit: analyse only the first ``limit`` input records — reads (single-end) or read
            pairs (paired) — then stop, without decompressing the rest of the file. ``None``
            reads everything. The mate-order / truncation checks below still run on every
            record actually read; a truncation *beyond* ``limit`` is simply never reached —
            that is the intent of a head-style limit, not a hole in the check.

    Raises:
        ValueError: if the two files disagree on read names or record count. This is not paranoia:
            a truncated R2 makes ``zip`` stop early and silently analyse a prefix, and a shuffled R2
            pairs mate 1 of one fragment with mate 2 of another. Both were observed in this project's
            own data and produced a *published* false discovery (a spurious R2-only blind spot) that
            had to be retracted. A pair of FASTQs is an assertion; check it.
    """
    # dnaio only on the unlimited path. `limit` is a HEAD: a truncation beyond it must never be
    # reached, which is why `read_pairs(r1, r2, limit=2)` succeeds on a file pair that diverges at
    # record 3. dnaio validates pairing as it fills its own buffers, so it sees -- and raises on --
    # a divergence the caller asked never to reach. A limited run is small by construction, so it
    # gives up nothing that matters to take the pure-Python path here.
    if _dnaio is not None and limit is None:
        yield from _read_pairs_dnaio(r1, r2, reconstruct=reconstruct, with_qual=with_qual)
        return

    if r2 is None:
        it = seqio.read_sequences(r1, with_qual=with_qual)
        if with_qual:                       # FASTA has no quality: `read_sequences` yields None
            it = ((i, s, q or "") for i, s, q in it)
        yield from islice(it, limit) if limit is not None else it
        return

    _stem = frag_stem
    # zip_longest, not zip: plain zip stops at the shorter file AND consumes one extra record from the
    # longer one, so a truncated mate file is invisible both during and after the loop.
    # Quality is read only when reconstructing (merge_pair's tie-break needs it) -- the default path
    # keeps discarding it, so it costs nothing.
    wq = reconstruct or with_qual
    n = 0
    for n, (a, b) in enumerate(zip_longest(seqio.read_sequences(r1, with_qual=wq),
                                           seqio.read_sequences(r2, with_qual=wq))):
        if limit is not None and n >= limit:
            break
        if a is None or b is None:
            raise ValueError(
                f"R1 and R2 differ in length (diverge at record {n}); one file is truncated.")
        if wq:
            (i1, s1, q1), (i2, s2, q2) = a, b
        else:
            (i1, s1), (i2, s2) = a, b
        if _stem(i1) != _stem(i2):
            raise ValueError(
                f"R1/R2 mate mismatch at record {n}: {i1!r} vs {i2!r}. "
                f"The FASTQs are not in the same order.")
        if reconstruct:
            merged = merge_pair(s1, s2, q1=q1, q2=q2)
            if merged is not None:
                yield i1, merged
                continue
        if with_qual:
            yield f"{i1}/1", s1, q1 or ""
            yield f"{i2}/2", s2, q2 or ""
        else:
            yield f"{i1}/1", s1
            yield f"{i2}/2", s2


def _read_pairs_dnaio(r1, r2, *, reconstruct: bool,
                      with_qual: bool = False) -> Iterator[tuple]:
    """The same stream, parsed in C.

    Once ``--prefilter`` removed the search, reading became the largest single cost of a bulk run
    -- 65 % of a 0.024 %-receptor library, where before it was 3 % and explicitly not worth
    touching. Reworking the pure-Python loop bought 1.43x; dnaio is 2x on top of that including
    the mate tagging (3.25 vs 1.63 M records/s on a 1 M-pair fixture).

    It is used only because it makes the SAME two assertions this function has always made, and
    both are load-bearing: a truncated mate file and a shuffled mate file each produced a
    *published* false discovery in this project (a spurious R2-only blind spot) that had to be
    retracted. dnaio raises on both -- "There are more reads in file 1 than in file 2" and "Read
    name 'b' in file 1 does not match 'zzz'" -- and they are re-raised as ``ValueError`` here so
    callers and tests see the error type they always have.
    """
    try:
        if r2 is None:
            with _dnaio.open(str(r1)) as fh:
                for rec in fh:
                    yield ((rec.id, rec.sequence, rec.qualities or "") if with_qual
                           else (rec.id, rec.sequence))
            return
        with _dnaio.open(str(r1), str(r2)) as fh:
            for a, b in fh:
                if reconstruct:
                    merged = merge_pair(a.sequence, b.sequence,
                                        q1=a.qualities, q2=b.qualities)
                    if merged is not None:
                        yield a.id, merged
                        continue
                if with_qual:
                    yield f"{a.id}/1", a.sequence, a.qualities or ""
                    yield f"{b.id}/2", b.sequence, b.qualities or ""
                else:
                    yield f"{a.id}/1", a.sequence
                    yield f"{b.id}/2", b.sequence
    except _dnaio.FileFormatError as exc:
        # Re-word to the two phrases callers and tests match on. dnaio says "Reads are improperly
        # paired. There are more reads in file 1 than in file 2" and "Read name 'b' in file 1 does
        # not match 'zzz'"; this pipeline has always distinguished a TRUNCATED mate file from a
        # SHUFFLED one, because they are different data-integrity failures with different fixes.
        msg = str(exc)
        if "does not match" in msg:
            raise ValueError(f"R1/R2 mate mismatch: {msg}. "
                             f"The FASTQs are not in the same order.") from exc
        # ⛔ Only claim a truncation when dnaio actually reported one. It raises the same exception
        # type for a malformed RECORD, and a real example from this project's own data is a `+`
        # line that kept the original SRA description after the `@` line was renamed to carry a
        # mate suffix (`@SRR5233635.1/2` against `+SRR5233635.1 1 length=151`). Calling that "one
        # file is truncated" sends the reader hunting for a truncation that does not exist -- and a
        # fabricated data-integrity finding has already cost this project a retraction once.
        if "improperly paired" in msg or "more reads in file" in msg:
            raise ValueError(f"R1 and R2 differ in length; one file is truncated. {msg}") from exc
        raise ValueError(f"malformed FASTQ in {r1} / {r2}: {msg}") from exc
    except (EOFError, gzip.BadGzipFile) as exc:
        # Same surfacing as `seqio.read_sequences`: a bare EOFError from the gzip layer reaches
        # Typer as "Aborted." with no cause, which reads like a Ctrl-D rather than a bad input.
        raise ValueError(f"truncated or corrupt gzip input: {r1} / {r2}") from exc


@dataclass
class RnaseqReport:
    """Counts + timing for one ``map`` run (written as JSON with ``--report``)."""

    input: str
    organism: str
    #: Library shape, recorded here because nothing downstream can recover it: the AIRR carries
    #: only the reads that mapped, so its row count, its `sequence` lengths and its mate suffixes
    #: all describe the receptor subset rather than the library. `arda stats` reports them.
    paired: bool = False
    input_bytes: int = 0                # sum of the R1 (+ R2) file sizes AS SUBMITTED (gzip: compressed)
    read_length_min: int = 0
    read_length_max: int = 0
    read_length_mean: float = 0.0
    total_reads: int = 0
    mapped_reads: int = 0
    per_locus: dict[str, int] = field(default_factory=dict)
    # fragments whose every read lay wholly inside the constant region: receptor mRNA, no V(D)J
    constant_only_fragments: int = 0
    # fragments whose isotype was donated by a constant-region mate (free `c_call` on a gapped pair)
    isotype_from_mate: int = 0
    min_score: float = 0.0  # a mapped_reads count is uninterpretable without its cutoff
    threads: int = 0
    wall_seconds: float = 0.0
    reads_per_second: float = 0.0
    peak_rss_mb: float = 0.0          # whole-process high-water mark at stage end
    rss_gain_mb: float = 0.0          # how much THIS stage raised it
    # Two-pass segment search accounting, empty when it is off. `rescued` reads cost a full-
    # reference realignment; they are the price of the fast path never dropping one.
    segment_search: dict = field(default_factory=dict)
    # k-mer prefilter accounting, empty when it is off. `prefilter_passed / prefilter_seen` is the
    # only number that says whether it earned its keep on this library.
    prefilter_stats: dict = field(default_factory=dict)

    @property
    def mapped_fraction(self) -> float:
        return self.mapped_reads / self.total_reads if self.total_reads else 0.0

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["mapped_fraction"] = self.mapped_fraction
        return d


def map_rnaseq(
    r1: str | Path,
    output: str | Path,
    *,
    r2: str | Path | None = None,
    organism: str = "human",
    seqtype: str = "nt",
    threads: int = 0,
    sensitivity: float | None = None,
    strand: str = "both",
    chunk_size: int = _RNASEQ_CHUNK,
    map_d: bool = True,
    d_max_evalue: float | None = None,
    reconstruct: bool = False,
    min_score: float = _MIN_SCORE,
    max_seqs: int = mapper._MAX_SEQS,
    kmer: int | None = -1,
    drop_constant_only: bool = True,
    limit: int | None = None,
    emit_reads: str | Path | None = None,
    report_path: str | Path | None = None,
    two_pass: bool = False,
    adaptive: bool = False,
    fast_segments: bool = False,
    indel_rescue: bool = False,
    segment_only_v: bool = False,
    prefilter: bool = False,
    with_junction_quality: bool = False,
    with_mutation_quality: bool = False,
    shm: str = "framework",
    complete_junction_nt: int = 0,
) -> RnaseqReport:
    """Filter + map an RNA-seq FASTQ (single or paired); write mapped reads as AIRR.

    Args:
        r1: FASTA/FASTQ (gzip by ``.gz``). Single-end, or R1 of a pair.
        output: AIRR TSV of the mapped reads only (keyed by ``sequence_id``).
        r2: R2 FASTQ for paired input; ``None`` for single-end.
        min_score: drop mapped reads below this MMseqs2 bit score. ``0`` disables the
            filter (recall-max). See :data:`_MIN_SCORE` for the calibration.
        kmer: MMseqs2 ``-k``. The memory knob: the nucleotide prefilter allocates 4**k index
            entries, so the tool default k=15 costs ~8.4 GB peak RSS whatever else you set.
            arda defaults to 13 (~0.7 GB, and never slower). ``None`` = MMseqs2's default.
        max_seqs: MMseqs2 target hits per read. Does not change which reads are kept, only
            which V/J scaffold wins. See :data:`arda.annotate.mapper._MAX_SEQS`.
        limit: analyse only the first ``limit`` reads (single-end) / read pairs (paired), then
            stop — a native head, so a subsample no longer needs an external ``zcat | head |
            gzip`` round-trip. ``None`` maps the whole file.
        with_junction_quality: also emit a ``junction_quality`` column — the read's Phred+33
            string over exactly the bases of ``junction``, same orientation (see
            :func:`junction_quality`). OFF by default: it appends a non-schema column, so the
            default output stays byte-identical. This is the only place the FASTQ quality is
            still in hand — Stage 1 otherwise discards it — and it is what
            ``correct --min-junction-q`` gates on. Refused with ``reconstruct`` (a merged
            fragment has no single input quality string).
        emit_reads: optional path — write the mapped reads' sequences as FASTA
            (coding-strand oriented) for downstream handoff.
        report_path: optional path — write the :class:`RnaseqReport` as JSON.
        fast_segments: with ``two_pass``, answer the segment pass structurally instead of with
            `mmseqs search` -- 37x on that step, agreeing with it on .9997 of V alleles and .9998
            of J. It only NOMINATES candidates; the winner is still aligned against the full
            scaffold by MMseqs2, so the contract is that the AIRR output does not move. Ignored
            without ``two_pass``, since there is no segment pass to replace.
        indel_rescue: with ``fast_segments``, send reads carrying the two-diagonal signature of an
            indel to the GAPPED rescue path instead of resolving them on the fast path. One
            ungapped extension scores such a read only up to the indel, so its segment score is
            systematically low. Measured on 341,294 real IGH mates: 3.18 % of reads carry a V
            indel, rising to 8.00 % below 90 % V identity. Reroutes, never drops.
        adaptive: cap alignments per read and re-search only the reads whose capped score is
            low (:func:`arda.annotate.mapper._extend_uncertain`). Measured 2.17x on 1 M bulk reads
            with **zero reads lost** — but read preservation is not the whole guarantee.
            **OFF by default**: on the real-read fixture it also changes `junction_aa` on 3 of 453
            reads, and two of them scored 128 and 131, far above the 90-bit trigger. So a high
            score does NOT certify that the best alignment was found, and the trigger cannot be
            calibrated on score alone. Opt in only where a junction-level difference is acceptable.
        two_pass: use the segment reference to shortlist a single V×J scaffold per read before
            aligning (:func:`arda.annotate.mapper._segment_best_hits`). Reads it cannot resolve
            are realigned against the full reference, so nothing is dropped — see
            :mod:`arda.annotate.shortlist`. Requires ``segments.fasta`` (written by
            ``arda build-index``); silently falls back to the one-pass search when it is absent.

            ⛔ **The win is set by whether reads SPAN V INTO J, not by the library type**, and the
            predictor is ``fast_fraction`` in the report. Measured: 3.51× on a TCR amplicon (fast
            path 85 %), 2.96× on a 100 %-receptor human TRB set (95.6 %), 2.64× on mouse TRA
            (89 %) — but **1.03× slower** on the human IGH leg of that *same* 100 %-receptor
            dataset (16.3 %: those reads cover V and stop short of the short IGHJ target), and
            0.762× on 2.74 %-receptor bulk (5 %). Off by default because no library type predicts
            it; run a sample and read ``fast_fraction``.

    Returns:
        The run :class:`RnaseqReport` (also printed by the CLI).
    """
    want_qual = with_junction_quality or with_mutation_quality
    if with_mutation_quality and reconstruct:
        raise ValueError("--mutation-quality cannot be combined with --reconstruct: a merged "
                         "fragment has no single input quality string")
    if with_junction_quality and reconstruct:
        # A merged fragment's bases come from two reads, so no input quality string describes it.
        # Emitting one of the two would put a quality beside bases it does not belong to -- the
        # exact silent-misalignment failure `junction_quality` verifies against. Refuse instead.
        raise ValueError("--junction-quality cannot be combined with --reconstruct: a merged "
                         "fragment has no single input quality string")
    if shm not in SHM_MODES:
        # Reject rather than fall through to the default: a mode that is accepted and silently
        # does nothing is the failure this project keeps hitting.
        raise ValueError(f"--shm must be one of {SHM_MODES}, got {shm!r}")
    output = Path(output)
    ref, target_db, threads, sens, mm_strand = mapper._prep(
        organism, seqtype, threads, sensitivity, strand)

    # Built and parsed ONCE for the whole run, not per chunk: the segment DB is an mmseqs build
    # and `combinations.tsv` is 550 KB. Both are read-only afterwards.
    segment_db = combos = None
    if two_pass:
        segment_db = mapper._cached_segment_db(ref, organism)
        if segment_db is None:
            logger.warning("--two-pass: no segments.fasta for %s (run `arda build-index`); "
                           "falling back to the one-pass search", organism)
        else:
            from ..annotate.shortlist import load_combinations
            combos = load_combinations(ref.target_fasta.parent / "combinations.tsv")

    chunks: queue.Queue = queue.Queue(maxsize=2)
    # The reader runs in a daemon thread, so anything it raises -- a missing FASTQ, a shuffled mate
    # file, a truncated gzip -- dies there unheard while `finally` still posts the sentinel. The main
    # loop then sees a clean end-of-stream and reports "0/0 reads mapped", exit 0. A pipeline reads
    # that as "this sample has no receptor reads". Capture it and re-raise on the consumer side.
    reader_exc: list[BaseException] = []

    def reader():
        try:
            pairs = read_pairs(r1, r2, reconstruct=reconstruct, limit=limit,
                               with_qual=want_qual)
            for chunk in chunked_fragments(pairs, chunk_size):
                chunks.put(chunk)
        except BaseException as exc:  # noqa: BLE001 — re-raised in the consumer
            reader_exc.append(exc)
        finally:
            chunks.put(None)  # sentinel

    t = threading.Thread(target=reader, daemon=True)
    t.start()

    report = RnaseqReport(input=str(r1), organism=organism, threads=threads,
                          min_score=min_score, paired=r2 is not None,
                          input_bytes=sum(Path(p).stat().st_size
                                          for p in (r1, r2) if p and Path(p).exists()))
    reads_fh = open(emit_reads, "w") if emit_reads else None
    stage = Stage()
    tick = Throttle()
    _len_sum = 0
    # Every extra goes at the END, in a fixed order, so a consumer reading the shipped set by
    # position is unaffected whichever combination is on.
    extra_cols: tuple[str, ...] = (
        ((JUNCTION_QUALITY,) if with_junction_quality else ())
        + (MUTATION_QUALITY if with_mutation_quality else ())
        + (FULL_COLUMNS if shm == "both" else ()))
    logger.info("map: %s%s -> %s | %d threads, chunk %d, min_score %g",
                r1, f" + {r2}" if r2 else "", output, threads, chunk_size, min_score)
    try:
        with open(output, "w") as fh:
            fh.write(airr_header(extra_cols) + "\n")

            def flush(batch: list) -> None:
                """Search one batch and write its mapped reads."""
                quals = {}
                if want_qual:
                    # Split the quality off HERE, not in the reader: `_annotate_chunk`,
                    # `write_fasta` and `prefilter.keep_records` all take `(id, sequence)` pairs,
                    # and a third element would reach every one of them.
                    quals = {r[0]: r[2] for r in batch}
                    batch = [(r[0], r[1]) for r in batch]
                keep = mapper._annotate_chunk(
                    batch, ref, target_db, seqtype, threads=threads,
                    sensitivity=sens, mm_strand=mm_strand, map_d=map_d,
                    d_max_evalue=d_max_evalue,
                    mapped_only=True, max_seqs=max_seqs, kmer=kmer,
                    segment_db=segment_db, combos=combos, adaptive=adaptive,
                    fast_segments=fast_segments,
                    indel_rescue=indel_rescue,
                    segment_only_v=segment_only_v, shm=shm,
                    complete_junction_nt=complete_junction_nt,
                    report=report.segment_search if segment_db else None)
                if drop_constant_only:
                    keep, n_drop, n_iso = _apply_constant_rule(keep)
                    report.constant_only_fragments += n_drop
                    report.isotype_from_mate += n_iso
                if min_score > 0:
                    keep = [r for r in keep
                            if float(r.get("mmseqs2_score") or 0) >= min_score]
                if not keep:
                    return
                if with_junction_quality:
                    for r in keep:
                        r[JUNCTION_QUALITY] = junction_quality(r, quals.get(r["sequence_id"], ""))
                if with_mutation_quality:
                    for r in keep:
                        r.update(mutation_quality(r, quals.get(r["sequence_id"], "")))
                fh.write(format_rows(keep, extra_cols))
                report.mapped_reads += len(keep)
                for r in keep:
                    loc = r.get("locus") or "?"
                    report.per_locus[loc] = report.per_locus.get(loc, 0) + 1
                    if reads_fh is not None:
                        reads_fh.write(f">{r['sequence_id']}\n{r['sequence']}\n")

            # The prefilter runs HERE, in the consumer, not in the reader thread. Moving it into
            # the reader to overlap it with `mmseqs search` was tried and measured WORSE
            # (SRR10611239 7.79 -> 9.12 s), for two reasons worth writing down:
            #
            #   * the overlap already exists, just split differently -- parsing is in the reader
            #     thread and the prefilter here, so the two run concurrently. Putting both in one
            #     thread serialises them, and on a cold library the prefilter (1.3 s) dwarfs the
            #     search (1.35 s of a 7.79 s run), so that loss is the bigger term.
            #   * at a 2 % pass rate there is only ever ONE batch to overlap. Filling 400 k
            #     survivors on SRR8363894 takes 19.4 M reads scanned -- more than the whole file --
            #     so the reader must consume everything before the first search can start.
            #
            # READ chunk size and SEARCH batch size are not the same thing once the prefilter is
            # on. Reading stays chunked to bound memory, but a prefiltered chunk is tiny -- 0.47 %
            # of reads survive on a 0.024 %-receptor library -- and every `mmseqs search` call
            # costs ~0.7 s of fixed setup whatever it is given. Ten near-empty searches were 25.8 s
            # against 13.4 s for one (SRR10611239). So survivors accumulate until they amount to a
            # full chunk's worth of work, and only then is MMseqs2 invoked.
            #
            # Flushing only ever happens on a chunk boundary, never inside one: `chunked_fragments`
            # guarantees a fragment's mates land in the same chunk, and splitting them is a bug
            # this pipeline has already shipped once under --reconstruct.
            pending: list = []
            while True:
                chunk = chunks.get()
                if chunk is None:
                    if reader_exc:
                        raise reader_exc[0]
                    break
                report.total_reads += len(chunk)
                # Read length, over the WHOLE library rather than the mapped subset -- the AIRR
                # only holds reads that mapped, so nothing downstream can recover this. One pass
                # over a list already in memory; `_len_sum` divides at the end, so a 300 M-read
                # library costs one int per chunk rather than a length histogram.
                lens = [len(r[1]) for r in chunk]
                if lens:
                    _len_sum += sum(lens)
                    lo, hi = min(lens), max(lens)
                    report.read_length_min = lo if not report.read_length_min else min(
                        report.read_length_min, lo)
                    report.read_length_max = max(report.read_length_max, hi)
                if tick.ready():
                    logger.info("map: %s reads read, %s mapped (%.2f %%), %.0f reads/s, %.0f MB",
                                f"{report.total_reads:,}", f"{report.mapped_reads:,}",
                                100.0 * report.mapped_reads / max(1, report.total_reads),
                                report.total_reads / max(1e-9, stage.wall_seconds), peak_rss_mb())
                if prefilter and seqtype == "nt":
                    from ..prefilter import keep_records  # noqa: PLC0415 — optional native ext
                    if want_qual:
                        qmap = {r[0]: r[2] for r in chunk}
                        survivors = [(i, s, qmap[i]) for i, s in
                                     keep_records([(r[0], r[1]) for r in chunk],
                                                  ref.target_fasta, threads=threads)]
                    else:
                        survivors = keep_records(chunk, ref.target_fasta, threads=threads)
                    report.prefilter_stats["seen"] = (
                        report.prefilter_stats.get("seen", 0) + len(chunk))
                    report.prefilter_stats["passed"] = (
                        report.prefilter_stats.get("passed", 0) + len(survivors))
                else:
                    survivors = chunk
                pending.extend(survivors)
                if len(pending) >= chunk_size:
                    flush(pending)
                    pending = []
            if pending:
                flush(pending)
    finally:
        if reads_fh is not None:
            reads_fh.close()
    t.join()

    stage.finish(report)   # wall_seconds / peak_rss_mb / rss_gain_mb, same definition as Stages 2-3
    report.reads_per_second = round(
        report.total_reads / report.wall_seconds if report.wall_seconds else 0.0, 1)
    report.read_length_mean = round(_len_sum / report.total_reads, 1) if report.total_reads else 0.0
    logger.info("map: %s/%s reads mapped (%.2f %%) in %.1f s (%.0f reads/s), peak %.0f MB; loci=%s",
                f"{report.mapped_reads:,}", f"{report.total_reads:,}",
                report.mapped_fraction * 100, report.wall_seconds, report.reads_per_second,
                report.peak_rss_mb, report.per_locus)
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report
