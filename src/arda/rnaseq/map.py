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

import json
import logging
import queue
import threading
from itertools import islice, zip_longest
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from ._res import Stage
from ..annotate import io as seqio
from ..annotate import mapper
from ..annotate.airr_out import airr_header, format_rows
from ..refbuild.translate import reverse_complement

logger = logging.getLogger(__name__)

__all__ = ["map_rnaseq", "read_pairs", "merge_pair", "RnaseqReport"]

_MERGE_ANCHOR = 12  # exact k-mer used to locate the R1/rc(R2) overlap in O(len)


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
    """
    batch: list = []
    prev: str | None = None
    for rec in records:
        stem = frag_stem(rec[0])
        if batch and len(batch) >= size and stem != prev:
            yield batch
            batch = []
        batch.append(rec)
        prev = stem
    if batch:
        yield batch


def read_pairs(r1: str | Path, r2: str | Path | None = None,
               *, reconstruct: bool = False, limit: int | None = None) -> Iterator[tuple[str, str]]:
    """Stream ``(id, sequence)`` reads for single-end (``r1`` only) or paired input.

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
    if r2 is None:
        it = seqio.read_sequences(r1)
        yield from islice(it, limit) if limit is not None else it
        return

    _stem = frag_stem
    # zip_longest, not zip: plain zip stops at the shorter file AND consumes one extra record from the
    # longer one, so a truncated mate file is invisible both during and after the loop.
    # Quality is read only when reconstructing (merge_pair's tie-break needs it) -- the default path
    # keeps discarding it, so it costs nothing.
    wq = reconstruct
    n = 0
    for n, (a, b) in enumerate(zip_longest(seqio.read_sequences(r1, with_qual=wq),
                                           seqio.read_sequences(r2, with_qual=wq))):
        if limit is not None and n >= limit:
            break
        if a is None or b is None:
            raise ValueError(
                f"R1 and R2 differ in length (diverge at record {n}); one file is truncated.")
        if reconstruct:
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
        yield f"{i1}/1", s1
        yield f"{i2}/2", s2


@dataclass
class RnaseqReport:
    """Counts + timing for one ``map`` run (written as JSON with ``--report``)."""

    input: str
    organism: str
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
    prefilter: bool = False,
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
        emit_reads: optional path — write the mapped reads' sequences as FASTA
            (coding-strand oriented) for downstream handoff.
        report_path: optional path — write the :class:`RnaseqReport` as JSON.
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
            :mod:`arda.annotate.shortlist`. Off by default: the win scales with the library's
            receptor fraction, so it pays on amplicon and barely moves a 0.0003 %-receptor
            negative. Requires ``segments.fasta`` (written by ``arda build-index``); silently
            falls back to the one-pass search when it is absent.

    Returns:
        The run :class:`RnaseqReport` (also printed by the CLI).
    """
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
            pairs = read_pairs(r1, r2, reconstruct=reconstruct, limit=limit)
            for chunk in chunked_fragments(pairs, chunk_size):
                chunks.put(chunk)
        except BaseException as exc:  # noqa: BLE001 — re-raised in the consumer
            reader_exc.append(exc)
        finally:
            chunks.put(None)  # sentinel

    t = threading.Thread(target=reader, daemon=True)
    t.start()

    report = RnaseqReport(input=str(r1), organism=organism, threads=threads,
                          min_score=min_score)
    reads_fh = open(emit_reads, "w") if emit_reads else None
    stage = Stage()
    try:
        with open(output, "w") as fh:
            fh.write(airr_header() + "\n")

            def flush(batch: list) -> None:
                """Search one batch and write its mapped reads."""
                keep = mapper._annotate_chunk(
                    batch, ref, target_db, seqtype, threads=threads,
                    sensitivity=sens, mm_strand=mm_strand, map_d=map_d,
                    mapped_only=True, max_seqs=max_seqs, kmer=kmer,
                    segment_db=segment_db, combos=combos, adaptive=adaptive,
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
                fh.write(format_rows(keep))
                report.mapped_reads += len(keep)
                for r in keep:
                    loc = r.get("locus") or "?"
                    report.per_locus[loc] = report.per_locus.get(loc, 0) + 1
                    if reads_fh is not None:
                        reads_fh.write(f">{r['sequence_id']}\n{r['sequence']}\n")

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
                if prefilter and seqtype == "nt":
                    from ..prefilter import keep_mask  # noqa: PLC0415 — optional native ext
                    mask = keep_mask(chunk, ref.target_fasta, threads=threads)
                    survivors = [r for r, m in zip(chunk, mask) if m]
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
    if report_path is not None:
        Path(report_path).write_text(json.dumps(report.as_dict(), indent=2) + "\n")
    return report
