"""Runtime annotation: map input sequences to the reference and transfer markup.

Pipeline: read input (FASTA/FASTQ) -> MMseqs2 search against the curated scaffold
DB -> best hit per query -> project reference region markup onto the query (C++
hot path) -> AIRR TSV.
"""

from __future__ import annotations

import logging
import os
import queue
import shutil
import tempfile
import threading
from pathlib import Path

import polars as pl

from .. import mmseqs
from .._locking import build_lock
from ..paths import data_dir, vdj_dir
from ..refbuild.translate import reverse_complement
from . import io as seqio
from .reference import load_reference, Reference
from .transfer import transfer_hit, AIRR_COLUMNS
from .airr_out import airr_header, format_rows

logger = logging.getLogger(__name__)

__all__ = ["annotate_file", "annotate_records", "build_index"]

# Streaming defaults: process the input in bounded chunks so memory stays flat
# regardless of input size (a 30M-read FASTQ never fully loads).
_CHUNK_SIZE = 50_000

_SEARCH_TYPE = {"nt": mmseqs.SEARCH_TYPE_NUCLEOTIDE, "aa": mmseqs.SEARCH_TYPE_PROTEIN}

# Tuned defaults (see memory/mmseqs-params.md). Short germline-similar queries:
# moderately high sensitivity, keep only a few best hits, no coverage filter so
# partial RNA-seq reads still map. The whole reference (all loci) is one DB, so a
# single search annotates mixed bulk RNA-seq across all loci at once.
_SENSITIVITY = {"nt": 7.0, "aa": 7.0}
# Max MMseqs2 target hits kept per read, i.e. how many V-J scaffolds may compete before
# best-hit selection. This is NOT a cosmetic knob: the V(D)J loci are paralog-dense, and at 50
# the true scaffold was being truncated out of the candidate list *before* `_best_hits` ever saw
# it. Measured on PRJNA371303 TRA amplicon (20k reads, arda-benchmark RESULTS.md §16):
#
#   max_seqs   V-gene concord.   J-gene   junction_aa exact   wall(amplicon)  wall(RNA-seq)
#       50         83.4 %        76.0 %       84.77 %            7.6 s           48 s
#      300         96.1 %        88.4 %       87.37 %           13.4 s           58 s
#
# The mapped read SET is identical either way (bit scores can only rise with more candidates —
# verified: 68 reads scored higher at 300, 0 scored lower), so the filter and the `--min-score 75`
# calibration are untouched; only the V/J/junction CALLS improve. This is why the V+J re-ranking
# experiment (RESULTS §5) was a wash: no re-ranking can recover a scaffold that was never a
# candidate. Cost is +21 % wall and 2.3x peak RSS (820 -> 1908 MB at chunk 200k); lower it via
# `--max-seqs` if memory-bound.
_MAX_SEQS = 300

# How many exactly-tied V alleles may contribute a scaffold to the two-pass alignment. V alleles of
# one gene differ by a nucleotide or two, so a 100 nt read commonly cannot separate them; aligning
# all of the tied ones lets `_best_hits` break the tie by the same rule the one-pass search uses.
# 8 is a degeneracy guard, not a tuning knob -- a read tying against more than that is not
# resolvable at allele level by any means.
_MAX_TIED_V = 8

# `--max-seqs` for the SEGMENT pass, which is a structurally different database from the one 300
# was calibrated on and must not inherit its value.
#
# 300 exists because the full reference is a CROSS-PRODUCT: 15,414 V×J scaffolds from 1,244
# distinct segments, so a read covering its V matches ~277 scaffolds at near-identical scores and
# truncating that list truncates the true one. Measured on the TRA amplicon, `--max-seqs` against
# the full reference: at 150 `junction_aa` agreement is already 98.20 %, at 75 it is 85.30 %, and
# at 25 it is 56.17 % -- 44 % of junctions moved while the mapped read count falls 0.013 %, which
# hides the whole thing.
#
# The segment reference has no cross-product. A read's V competes only against the alleles of its
# own gene and its close paralogs, and the pass needs exactly three things: the best V, the best J,
# and up to `_MAX_TIED_V` exactly-tied V alleles. 300 of 1,244 targets is 24 % of the whole
# database per query to answer that -- against 1.9 % for 300 of 15,414 on the full reference. The
# two knobs share a name and nothing else.
#
# Anything wrong here is bounded rather than silent: a mis-implied scaffold either loses the
# `(V, J)` lookup and goes to `_full_rescue`, or is caught by the tied-V expansion. Set to a value
# that comfortably clears `_MAX_TIED_V` on both sides.
_SEGMENT_MAX_SEQS = 50

# Adaptive search. `--max-accept` bounds how many alignments mmseqs performs per query before it
# stops; it is UNBOUNDED by default, so arda aligns every hitting read against all ~300 of its
# prefilter candidates and then keeps exactly one. Capping it is the single largest lever on the
# align term, which is 75 % of bulk search wall (`wall = reads/46,353 + hits/350`).
#
# The cap is not free on its own: mmseqs orders candidates by prefilter score, which predicts the
# true best alignment only 55.8 % of the time, so a capped search can stop before reaching a read's
# real scaffold. But the reads that suffer are identifiable from the output -- every read lost at
# `--max-accept 40` scored 75-83, i.e. just above the `--min-score 75` cutoff. A read returning 300
# bits has plainly found its scaffold; one returning 76 may not have.
#
# So: cap everything, then re-search UNCAPPED only the reads whose capped best score is below
# `_ADAPTIVE_TRIGGER`. Measured on 1 M real bulk reads: 2.17x end to end with **zero reads lost**,
# where the uncertain set was 4,816 reads (0.5 % of the library) costing 3.9 s against a 50 s
# saving. No single `--max-accept` value achieves this -- its own lossless point is 1.25x.
#
# ⚠ OFF BY DEFAULT, because preserving the read SET is not the whole guarantee. On the real-read
# fixture the adaptive search also changes `junction_aa` on 3 of 453 reads -- and two of them
# scored **128 and 131**, far above the 90-bit trigger. A high score therefore does NOT certify
# that the best alignment was found: a read can be comfortably above threshold on a scaffold that
# is not its best, and the junction moves even though the read is kept. The bulk measurement did
# not catch this because it compared read sets and winning targets, not junctions.
#
# Any future calibration of `_ADAPTIVE_TRIGGER` has to be judged on junction identity, not on
# read survival -- and since the counter-examples sit at 128-131 bits, a score-only trigger may
# not be calibratable at all.
_MAX_ACCEPT = 40
_ADAPTIVE_TRIGGER = 90.0

# MMseqs2's nucleotide prefilter allocates a k-mer index table of 4**k entries, so peak RSS tracks
# 4**k * 8 B almost exactly and is INDEPENDENT of database size, --max-seqs, --chunk-size and
# --threads. Measured on the V+J | J+C reference, 100 k reads, 8 threads:
#
#     k       11      12      13      14      15
#     RSS   202 MB  298 MB  697 MB  2313 MB  ~8.4 GB
#     wall  3.76 s  3.60 s  4.63 s  4.15 s      --
#
# and recall/precision are INVARIANT over k=11..14 (recall 1.0000, precision .9463-.9469, FP 269-270).
# So k is a pure memory/speed knob in this range: 12 is both the fastest measured and 2.3x smaller
# than 13. arda-benchmark OPTIMIZATION.md §6.3c. Nucleotide only: the aa prefilter is a different index.
_KMER = {"nt": 12, "aa": None}


# Columns of DEFAULT_FORMAT_OUTPUT that are genuinely text; everything else is numeric.
_STR_COLS = frozenset({"query", "target", "cigar", "qaln", "taln"})


def _dbtype(seqtype: str) -> int:
    return 2 if seqtype == "nt" else 1


def _committed_index(organism: str, seqtype: str, target_fasta: Path | None = None) -> Path | None:
    """A precompiled mmseqs DB shipped in ``database/`` — if it is version-matched **and current**.

    The staleness check is not optional. ``build-db`` rewrites ``alleles.fasta`` in place, and without
    this every later run silently searched the *previous* reference; the only symptom was a result that
    quietly refused to change. Found exactly that way: a rebuild adding 345 constant-region scaffolds
    produced zero constant-region hits. The ``data/mmseqs_db`` fallback below has always compared
    mtimes -- this path did not.
    """
    d = vdj_dir(organism) / "mmseqs" / seqtype
    db, ver = d / "db", d / "VERSION"
    if not (db.exists() and ver.exists()):
        return None
    if target_fasta is not None and db.stat().st_mtime < Path(target_fasta).stat().st_mtime:
        return None  # the reference was rebuilt after this index was compiled
    try:
        if mmseqs.versions_compatible(ver.read_text(), mmseqs.version()):
            return db
    except Exception:  # noqa: BLE001 — mmseqs missing; fall through to rebuild
        return None
    return None


def _createdb_atomic(target_fasta: Path, db: Path, dbtype: int) -> None:
    """Build an mmseqs DB at ``db``, safely when several arda processes start at once.

    `mmseqs createdb` writes ~6 sibling files (``db``, ``db.index``, ``db_h``, ``db_h.index``,
    ``db.dbtype``, ``db.lookup``). Building them in place is not safe, and arda is routinely run
    concurrently against the SAME reference (see :mod:`arda._locking`): every process finds no index
    and starts building it into the same paths.

    They then interleave writes -- and `build_index` additionally unlinks the files another process
    is mid-read on. The observed failure is silent and total: a 0-byte ``db``, after which every read
    fails to map (``0/200000 reads mapped, loci={}``) with no error and a clean exit code.

    So: hold the build lock, build into a private temp dir, and move the finished files into place
    with ``db`` **last** -- readers test ``db.exists()``, so it must not appear until its siblings are
    all there. A killed builder leaves a temp dir, never a half-built DB that looks complete.

    "Done" means a **current** db, not merely a present one: an existing db OLDER than
    ``target_fasta`` is stale and must be rebuilt. Gating on bare existence made a stale cache
    permanent -- ``_cached_target_db`` decides to rebuild on mtime, then this guard sees the old
    file, calls it done, and skips. The reference was rebuilt (e.g. `build-db`, or a pulled
    `alleles.fasta`) but every run kept searching the previous scaffolds, silently shifting markup.
    """
    def _current() -> bool:
        return db.exists() and db.stat().st_mtime >= target_fasta.stat().st_mtime

    with build_lock(db.parent / f".{db.name}.lock", done=_current) as ours:
        if not ours:
            return
        tmp = db.parent / f".{db.name}.tmp.{os.getpid()}"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)
        try:
            mmseqs.createdb(target_fasta, tmp / db.name, dbtype=dbtype)
            built = sorted(tmp.glob(f"{db.name}*"))
            # `db` itself last: it is the existence check every reader uses.
            for f in sorted(built, key=lambda p: p.name == db.name):
                os.replace(f, db.parent / f.name)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def _cached_target_db(target_fasta: Path, organism: str, seqtype: str) -> Path:
    """Resolve the mmseqs target DB for the reference scaffolds.

    Prefers the precompiled DB shipped in ``database/vdj/<org>/mmseqs/<seqtype>``
    (used out of the box when its mmseqs version matches). Otherwise builds once
    into a ``data/mmseqs_db`` cache and reuses it (rebuilt if the FASTA is newer).
    """
    committed = _committed_index(organism, seqtype, target_fasta)
    if committed is not None:
        return committed
    cache = data_dir() / "mmseqs_db" / f"{organism}_{seqtype}"
    cache.mkdir(parents=True, exist_ok=True)
    db = cache / "db"
    if not db.exists() or db.stat().st_mtime < Path(target_fasta).stat().st_mtime:
        _createdb_atomic(Path(target_fasta), db, _dbtype(seqtype))
    return db


def _cached_segment_db(ref: Reference, organism: str) -> Path | None:
    """The mmseqs DB for the segment reference, or ``None`` if it has not been built.

    Nothing precompiled ships for this one -- `segments.fasta` is generated by ``build-index``
    (see :mod:`arda.refbuild.segments`), so this always goes through the local cache. Cached
    per *chunk-independent* key so a 400 k-read run builds it once, not once per chunk;
    `_createdb_atomic` holds the build lock, which is what makes that safe for the Nextflow
    process-per-sample and SLURM task-per-shard layouts.

    Returns ``None`` rather than raising when `segments.fasta` is absent: the two-pass is then
    simply unavailable and the caller falls back to the one-pass search.
    """
    fasta = ref.target_fasta.parent / "segments.fasta"
    if not fasta.exists():
        return None
    cache = data_dir() / "mmseqs_db" / f"{organism}_segments"
    cache.mkdir(parents=True, exist_ok=True)
    db = cache / "db"
    if not db.exists() or db.stat().st_mtime < fasta.stat().st_mtime:
        _createdb_atomic(fasta, db, _dbtype("nt"))
    return db


def build_index(organism: str = "all", *, force: bool = False) -> None:
    """(Re)build the precompiled mmseqs DBs shipped under ``database/``.

    Writes ``database/vdj/<org>/mmseqs/<seqtype>/db*`` + a ``VERSION`` marker so
    the runtime can use them out of the box (and detect a mmseqs-version mismatch).
    Skips up-to-date DBs unless ``force``.
    """
    from ..igblast import SUPPORTED_ORGANISMS
    orgs = SUPPORTED_ORGANISMS if organism == "all" else (organism,)
    ver = mmseqs.version()
    for org in orgs:
        for seqtype in ("nt", "aa"):
            fasta = vdj_dir(org) / ("alleles.aa.fasta" if seqtype == "aa" else "alleles.fasta")
            if not fasta.exists():
                continue
            out = vdj_dir(org) / "mmseqs" / seqtype
            db, vfile = out / "db", out / "VERSION"
            # Skip only if the DB is version-matched AND newer than the FASTA. Matching on the mmseqs
            # version alone made `build-index` a silent no-op after `build-db` rewrote alleles.fasta:
            # the rebuilt reference was never compiled, and the runtime kept searching the old one.
            fresh = db.exists() and db.stat().st_mtime >= fasta.stat().st_mtime
            if fresh and not force and vfile.exists() and \
                    mmseqs.versions_compatible(vfile.read_text(), ver):
                continue
            out.mkdir(parents=True, exist_ok=True)
            # No unlink-then-rebuild: that deleted the files a concurrently-running arda was reading.
            # `_createdb_atomic` builds in a temp dir under a lock and moves the finished files in.
            vfile.unlink(missing_ok=True)     # drop the marker first: a DB without one is "unusable",
            if force:                         # so a crash mid-rebuild degrades to a rebuild, not a lie
                for stale in out.glob("db*"):
                    stale.unlink()
            _createdb_atomic(fasta, db, _dbtype(seqtype))
            vfile.write_text(ver + "\n")      # marker written LAST: it is what certifies the DB
            print(f"[arda] built mmseqs index {org}/{seqtype} ({ver})")

        # The segment reference (V / J / J+C as separate targets) is derived from `alleles.fasta` +
        # `markup.tsv` in well under a second, so it is generated here rather than committed or
        # shipped in the release tarball -- the same argument that already excludes the mmseqs
        # indexes. Its own mmseqs DB is built lazily by `_cached_target_db` on first use.
        alleles = vdj_dir(org) / "alleles.fasta"
        seg = vdj_dir(org) / "segments.fasta"
        if alleles.exists() and (force or not seg.exists()
                                 or seg.stat().st_mtime < alleles.stat().st_mtime):
            from ..refbuild.segments import build_segment_reference
            stats = build_segment_reference(org)
            print(f"[arda] built segment reference {org} ({stats.total} targets)")


def _best_hits(tsv: Path) -> dict[str, dict]:
    """Parse the mmseqs TSV and return the top-scoring hit per query.

    Best = highest whole-scaffold bit score. (A V-segment-restricted re-ranking was
    tried to fix V-paralog mis-assignment on short 3'-V-anchored reads; it only traded
    V accuracy for J and cost ~17 % throughput — the disagreement is largely irreducible
    paralog ambiguity — so max-bits is kept. See arda-benchmark RESULTS.md.)
    """
    cols = mmseqs.DEFAULT_FORMAT_OUTPUT.split(",")
    if tsv.stat().st_size == 0:  # no hits at all
        return {}
    # Typed, not `infer_schema_length=0`. Reading 17 columns as Utf8 and then casting one of them was
    # a large share of peak RSS before `mmseqs.top_hit` shrank this file from 194 MB to 1 MB; keep the
    # schema explicit so it stays cheap if `top_hit` is ever bypassed.
    schema = {c: (pl.Utf8 if c in _STR_COLS else pl.Float64) for c in cols}
    df = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=cols,
                     schema_overrides=schema)
    if df.height == 0:
        return {}
    df = df.with_columns(pl.col("bits").cast(pl.Float64, strict=False))
    # `target` is a tie-break key, and `maintain_order` is not decoration. Paralogous scaffolds
    # tie exactly on whole-scaffold bits (IGHV3-30*01 vs IGHV3-30-3*01, IGKV1-x vs IGKV3-x), and
    # polars' sort and `unique` are both unordered by default -- so the winner, and therefore the
    # emitted v_call, depended on the order mmseqs happened to list the alignments in. Measured:
    # every one of 25 tied queries flipped its call when the input rows were reversed.
    # That makes a run irreproducible, and it makes byte-identity between the single-node,
    # sharded and Nextflow paths unprovable -- a shard boundary changes row order, nothing else.
    # Lexicographically smallest target among exact ties: arbitrary, but the same everywhere.
    # An alignment row whose query NAME is empty cannot be attributed to a read. It happens on the
    # two-pass path: `_align_implied` searches a hand-built prefilter sub-DB, and `createsubdb`
    # does not carry the header/`.lookup` entries for everything in it, so `convertalis` emits a
    # blank query field for the uncovered ones. Left in, the null became a dict key and
    # `seqs[None]` raised `KeyError: None` -- `--two-pass` died ~90 % of the way through a 1 M-pair
    # amplicon run, after writing a partial output that looked like a completed one.
    #
    # Dropping them is the CORRECT repair, not a silent workaround: a read absent from this
    # mapping falls into `set(want) - set(hits)` in the caller and is realigned against the full
    # reference, which is the exactness guarantee the two-pass path is built around. It is warned
    # about rather than passed over, because a large count means the sub-DB is malformed.
    #
    # `bits` is cast with `strict=False`, so an unparseable score becomes null rather than raising
    # at parse time -- and a null then reaches `float(row["bits"])` in the J+C contest and dies
    # there instead. Both symptoms are the same malformed row, so both are filtered here, at the
    # one place that knows the row is unusable.
    n_before = df.height
    df = df.filter(pl.col("query").is_not_null() & (pl.col("query") != "")
                   & pl.col("bits").is_not_null())
    if df.height != n_before:
        logger.warning("%d alignment rows were unusable (no query id or no score) and were routed "
                       "to the full-reference rescue; a large count means the hand-built sub-DB "
                       "is missing lookup entries", n_before - df.height)
    if df.height == 0:
        return {}
    df = (df.sort(["bits", "target"], descending=[True, False], maintain_order=True)
            .unique(subset="query", keep="first", maintain_order=True))
    return {row["query"]: row for row in df.iter_rows(named=True)}


# The segment pass needs only these five fields. NOT `DEFAULT_FORMAT_OUTPUT`: that carries
# `cigar`/`qaln`/`taln` for up to `--max-seqs` rows per read, which is the 194 MB / 877 MB
# regression `top_hit` exists to prevent (see its docstring). Here we cannot reduce to one row
# per query first -- a junction-spanning read must contribute its best V AND its best J -- so the
# only lever is asking for fewer columns.
_SEGMENT_FORMAT = "query,target,bits,qstart,qend,tstart"


def _segment_rows(tsv: Path) -> list[dict]:
    """The alignment rows from a segment-pass convertalis TSV that can affect the answer.

    The segment search runs at the same ``--max-seqs`` as the full search, so an amplicon read
    covering its V emits a row for **every** V allele it clears threshold against -- measured at
    tens of rows per read. :func:`_segment_best_hits` consumes at most the top row per
    ``(query, side)`` plus up to ``_MAX_TIED_V`` exactly-tied V rows, and ``continue``s on all the
    rest, so materialising the whole file as Python dicts builds tens of millions of objects to
    throw nearly all of them away. That is why ``--two-pass`` peaked at **3.3 GB** on a 1 M-pair
    amplicon while the one-pass baseline it beats used 2.5 GB -- a mode doing strictly *less*
    alignment work should not cost more memory, and that was the tell.

    The reduction happens in polars instead, so only the rows the loop can act on are converted.
    **The loop's semantics are unchanged**: everything dropped here is a row it would have skipped.
    """
    if not tsv.exists() or tsv.stat().st_size == 0:
        return []
    cols = _SEGMENT_FORMAT.split(",")
    schema = {c: (pl.Utf8 if c in _STR_COLS else pl.Float64) for c in cols}
    df = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=cols,
                     schema_overrides=schema)
    # Deterministic: exact bit-score ties between paralogous targets are broken on `target`,
    # the same rule `_best_hits` uses. Resolving them by TSV row order (as an earlier draft did)
    # silently undid the determinism fix -- 25/25 tied queries flipped when rows were reversed.
    df = df.sort(["bits", "target"], descending=[True, False], maintain_order=True)

    # `V|`/`J|`/`JC|` prefix -> which side of the scaffold this row is evidence for. An
    # unrecognised prefix is dropped, matching the loop, which never treats one as a J.
    kind = pl.col("target").str.split("|").list.first()
    df = df.filter(kind.is_in(["V", "J", "JC"])).with_columns(
        _side=pl.when(kind == "V").then(pl.lit("V")).otherwise(pl.lit("J")))

    # ⛔ `over`, NOT `group_by`. A window function maps its result back to the ORIGINAL row
    # positions, so on an already-sorted frame it is deterministic; `group_by` is a multithreaded
    # hash aggregation, and using it unordered is precisely what made `correct` nondeterministic
    # across runs while the row count stayed stable.
    df = df.with_columns(
        _rank=pl.int_range(pl.len()).over(["query", "_side"]),
        _top=pl.col("bits").first().over(["query", "_side"]))
    keep = (pl.col("_rank") == 0) | (
        (pl.col("_side") == "V") & (pl.col("bits") == pl.col("_top"))
        & (pl.col("_rank") < _MAX_TIED_V))
    df = df.filter(keep).drop("_side", "_rank", "_top")
    return list(df.iter_rows(named=True))


def _subset_db(src_db: Path, ids: list[str], dst: Path,
               keys: dict[str, str] | None = None) -> Path:
    """A query sub-DB holding exactly ``ids``, keys and headers intact.

    Two things an earlier draft got wrong, both silent:

    * it read a ``.lookup`` beside the *output* DB, which ``createsubdb`` never writes, so every
      id missed and the entire fast path fell through to rescue while still producing correct
      output -- fast to miss, because nothing failed;
    * it subset only the sequence DB. ``createsubdb`` acts on ONE database, so without the
      matching ``_h`` call ``convertalis`` cannot print query *names* and emits numeric keys,
      which would key the result dict by something the caller never asked about.

    ``--id-mode 1`` is not usable for the header DB (``<db>_h`` has no ``.lookup`` of its own), so
    both calls go through numeric keys taken from the source DB's lookup. ``createsubdb``
    preserves those keys, which is why the caller can keep using the source mapping afterwards.

    It also writes ``<dst>.lookup`` itself -- verified against mmseqs 18-8cc5c, which emits
    ``dst``, ``dst.index``, ``dst.dbtype``, ``dst.lookup``, ``dst.source`` and the ``_h`` pair.
    An earlier version of this note claimed otherwise; writing the file by hand here on that
    assumption OVERWROTE mmseqs' own and broke the rescue path (490 ids "absent from lookup").
    So: the blank query names seen on the two-pass path do **not** come from a missing lookup,
    and their cause is still open.

    Raises:
        MMseqsError: if any id is absent from the source lookup. Dropping them silently is what
            would lose reads, and this function sits on the no-read-lost path.
    """
    keys = keys if keys is not None else _db_keys(src_db)
    missing = [i for i in ids if i not in keys]
    if missing:
        raise mmseqs.MMseqsError(
            f"{len(missing)} read id(s) absent from {src_db}.lookup, e.g. {missing[:3]}")
    listing = dst.parent / f"{dst.name}.keys"
    listing.write_text("\n".join(keys[i] for i in ids) + "\n")
    mmseqs.run(["createsubdb", str(listing), str(src_db), str(dst)])
    mmseqs.run(["createsubdb", str(listing), f"{src_db}_h", f"{dst}_h"])

    return dst


def _db_keys(db: Path) -> dict[str, str]:
    """``{fasta_identifier: numeric_db_key}`` from an mmseqs ``.lookup``.

    Needed because a hand-built prefilter DB addresses entries by numeric key, not by name.
    """
    out: dict[str, str] = {}
    lk = Path(f"{db}.lookup")
    if not lk.exists():
        return out
    with open(lk) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[1]] = parts[0]
    return out


def _v_end(ref: Reference | None, scaffold_id: str) -> int:
    """Scaffold nt position where the V germline ends, or 0 if unknown.

    This is what makes two scaffolds built from different-length alleles of one gene comparable:
    the read's offset into the scaffold differs by exactly this much.
    """
    if ref is None:
        return 0
    entry = ref.get(scaffold_id)
    return getattr(entry, "v_sequence_end", 0) if entry is not None else 0


def _align_implied(query_db: Path, full_db: Path, implied: dict[str, str],
                   seg_rows: dict[tuple[str, str], dict], seqs: dict[str, str], tmp: Path, *,
                   threads: int, side: str = "V", ref: Reference | None = None,
                   tag_prefix: str = "") -> tuple[dict[str, dict], set[str]]:
    """Align each read against ONLY the scaffold its (best V, best J) implies.

    Measured: 10,787 alignments instead of 3.04 M -- 5.36 s -> 0.044 s (122x) on a 20 k-read
    amplicon.

    **Strand is taken per read from what the segment pass observed, never assumed.** Library
    strandedness is a property of the prep, not a constant: measured across this project's own
    panel, TruSeq-stranded dUTP libraries put R1 antisense (`fr-firststrand`, Illumina's documented
    behaviour), Raji's DNBSEQ prep runs the *opposite* convention with R1 sense, TCR amplicons have
    **both** mates antisense (primer-defined, so not a transcriptome convention at all), and K562
    is unstranded at ~25/25/25/25. Four regimes in one panel.

    `mmseqs search` copes by running `extractframes` internally (doubling the query DB) and
    mapping back with `offsetalignment`. A hand-built prefilter DB has no such mechanism, so
    instead we reverse-complement the reads whose segment hit was on the minus strand, align them
    forward, and flip the reported coordinates back. The result is indistinguishable from the
    one-pass output: mmseqs also reports reverse hits with ``qstart > qend`` and aligned strings
    already on the coding strand, which is exactly what `_annotate_chunk` expects.

    ``side`` names which segment hit supplies the alignment diagonal — ``"V"`` for a V×J scaffold
    (the read enters it through V), ``"J"`` for a J+C scaffold (there is no V to enter through).

    Returns ``(hits, failed)``. Anything in ``failed`` is rescued by the caller, so a miss here
    costs time, never a read.
    """
    want = sorted(implied)
    qkey = _db_keys(query_db)
    tkey = _db_keys(full_db)

    fwd, rc, skipped = [], [], set()
    targets: dict[str, list[str]] = {}
    base: dict[str, str] = {}
    for qid in want:
        seg = seg_rows.get((qid, side))
        want_t = implied[qid]
        want_t = [want_t] if isinstance(want_t, str) else list(want_t)
        # A scaffold missing from the target DB is dropped from the candidate list rather than
        # failing the read -- only a read left with NO candidate is handed to the rescue pass.
        want_t = [x for x in want_t if x in tkey]
        if seg is None or qid not in qkey or not want_t:
            skipped.add(qid)
            continue
        targets[qid] = want_t
        base[qid] = want_t[0]        # the diagonal below is measured against this one
        qs, qe, ts = int(seg["qstart"]), int(seg["qend"]), int(seg["tstart"])
        (rc if qs > qe else fwd).append((qid, qs, qe, ts, float(seg["bits"])))

    hits: dict[str, dict] = {}
    # Forward-strand reads reuse the original query DB; reverse-strand reads need a private DB of
    # reverse complements, because `align` has no way to be told to flip.
    for tag, group, flip in ((f"{tag_prefix}f", fwd, False), (f"{tag_prefix}r", rc, True)):
        if not group:
            continue
        ids = [g[0] for g in group]
        if flip:
            rc_fa = tmp / f"rc_{tag}.fasta"
            with open(rc_fa, "w") as fh:
                for qid, *_ in group:
                    fh.write(f">{qid}\n{reverse_complement(seqs[qid])}\n")
            sub = tmp / f"impQ_{tag}"
            mmseqs.createdb(rc_fa, sub, dbtype=2)
            keys = _db_keys(sub)
        else:
            sub = _subset_db(query_db, ids, tmp / f"impQ_{tag}", keys=qkey)
            keys = qkey

        rows = []
        for qid, qs, qe, ts, bits in group:
            if flip:
                # 1-based position p on a length-L read maps to L - p + 1 on its reverse
                # complement, so a minus-strand hit (qs > qe) becomes a plus-strand one.
                qs = len(seqs[qid]) - qs + 1
            # One prefilter line per candidate, each with ITS OWN diagonal.
            #
            # A V×J scaffold is `V + pad + J` with the V left-aligned, so a read sits at the same
            # offset only in scaffolds whose V allele has the same length. Alleles of one gene do
            # NOT: 11 human V genes carry alleles of differing length, 8 of them differing by
            # >= 35 nt and one by 72. `mmseqs align` returns NOTHING once the diagonal is off by
            # more than ~35 nt, so a single shared diagonal silently drops the sibling scaffold --
            # including, on a real read, the very scaffold the one-pass search wins on
            # (TRBV29-1*03 at bits 121 lost, leaving *01 at 113).
            #
            # Shifting by the V-end difference corrects it. The correction is not always exact --
            # alleles can differ by more than pure 5' truncation, measured 2 nt of residual on the
            # TRBV29-1 case -- but the tolerance band is ~35 nt, so landing within a few nt is
            # what matters. A candidate whose geometry is unknown keeps the unshifted diagonal;
            # it is no worse off than before, and a miss costs a rescue, never a read.
            v_end0 = _v_end(ref, base[qid])
            for tgt in targets[qid]:
                shift = 0
                v_end = _v_end(ref, tgt)
                if v_end and v_end0:
                    shift = v_end - v_end0
                rows.append(f"{keys[qid]}\t{tkey[tgt]}\t{int(bits)}\t{qs - ts + shift}")
        pref_tsv = tmp / f"implied_{tag}.pref.tsv"
        pref_tsv.write_text("\n".join(rows) + "\n")
        pref_db = tmp / f"impPref_{tag}"
        mmseqs.run(["tsv2db", str(pref_tsv), str(pref_db), "--output-dbtype", "7"])
        aln = tmp / f"impAln_{tag}"
        # No `--search-type`: `mmseqs align` does not accept it (it infers nt from the DB type).
        mmseqs.run(["align", str(sub), str(full_db), str(pref_db), str(aln),
                    "-a", "--alignment-mode", "3", "--threads", str(threads)])
        out_tsv = tmp / f"implied_{tag}.tsv"
        # Reduce to one alignment per query BEFORE materialising it as text, as the one-pass path
        # does: `convertalis` writes `cigar`/`qaln`/`taln` per row, and offering every tied V
        # allele emits a measured mean of 3.45 rows per read (2.88x after the `_MAX_TIED_V` cap),
        # which rebuilds a fraction of the 194 MB -> 877 MB regression `top_hit` exists to prevent.
        #
        # `mmseqs align` DOES emit each query's results score-descending even from a hand-built
        # prefilter DB -- verified by removing this call: agreement fell 0.9956 -> 0.9735, i.e.
        # `_best_hits`' lexicographic tie-break disagrees with mmseqs' own ordering more often
        # than positional selection does. Keeping `top_hit` on both paths is what makes them
        # break exact ties by the SAME rule.
        mmseqs.convertalis(sub, full_db, mmseqs.top_hit(aln, tmp / f"impBest_{tag}"),
                           out_tsv, threads=threads, search_type=_SEARCH_TYPE["nt"])
        got = _best_hits(out_tsv)
        if flip:
            for qid, row in got.items():
                L = len(seqs[qid])
                row = dict(row)
                # Back to original-read coordinates, restoring the qstart > qend convention that
                # signals a minus-strand hit downstream.
                row["qstart"], row["qend"] = L - int(row["qstart"]) + 1, L - int(row["qend"]) + 1
                hits[qid] = row
        else:
            hits.update(got)

    return hits, skipped | (set(want) - set(hits))


def _full_rescue(query_db: Path, full_db: Path, ids: list[str], tmp: Path, *,
                 threads: int, sensitivity: float, mm_strand: int | None,
                 max_seqs: int, kmer: int | None, search_type: int) -> dict[str, dict]:
    """Realign unresolved reads against the FULL reference -- the exactness guarantee.

    This is what makes the fast path safe to enable: whatever the segment pass could not resolve
    gets precisely the treatment it would have received with no optimisation at all. Measured
    cost: 11 % of the new total on amplicon, 20 % on bulk.

    Raises:
        MMseqsError: if the sub-DB cannot be built. Returning ``{}`` here would silently drop the
            entire rescue set -- the exact failure this function exists to prevent -- so it is
            deliberately loud.
    """
    if not ids:
        return {}
    sub = _subset_db(query_db, sorted(ids), tmp / "rescQ")
    res, out_tsv = tmp / "rescRes", tmp / "rescue.tsv"
    mmseqs.search(sub, full_db, res, tmp / "resc_tmp", search_type=search_type,
                  sensitivity=sensitivity, max_seqs=max_seqs, threads=threads, kmer=kmer,
                  extra=(["--strand", str(mm_strand)] if mm_strand is not None else None))
    mmseqs.convertalis(sub, full_db, mmseqs.top_hit(res, tmp / "rescBest"),
                       out_tsv, threads=threads, search_type=search_type)
    return _best_hits(out_tsv)


def _segment_best_hits(
    query_db: Path, seg_db: Path, full_db: Path, tmp: Path, ref: Reference, *,
    threads: int, sensitivity: float, mm_strand: int | None, max_seqs: int, kmer: int | None,
    search_type: int, seqs: dict[str, str],
    combos: dict[tuple[str, str], str] | None = None,
) -> tuple[dict[str, dict], dict]:
    """Two-pass best hits: cheap segment pass, then align only the implied scaffold.

    Produces the same ``{query: hit_row}`` shape as the one-pass path, so everything downstream
    (strand handling, markup transfer, D mapping) is untouched.

    **No read is lost.** Whatever the segment pass cannot resolve to a single V×J scaffold --
    V-only, J-only, a pair absent from the reference, a reverse-strand hit, or a second-pass
    miss -- is realigned against the FULL reference exactly as today. The partition is asserted,
    not assumed, and `_full_rescue` raises rather than returning empty.
    """
    from .shortlist import load_combinations, shortlist

    seg_res, seg_tsv = tmp / "segRes", tmp / "seg.tsv"
    mmseqs.search(query_db, seg_db, seg_res, tmp / "seg_tmp",
                  search_type=search_type, sensitivity=sensitivity,
                  max_seqs=min(max_seqs, _SEGMENT_MAX_SEQS),
                  threads=threads, kmer=kmer,
                  extra=(["--strand", str(mm_strand)] if mm_strand is not None else None))
    # NOT `top_hit` here. `filterdb --extract-lines 1` keeps ONE row per query, which destroys
    # exactly the pairing this pass exists to find: a junction-spanning read must contribute its
    # best V AND its best J. Reducing first left `implied` at 0 -- correct output, zero speedup.
    mmseqs.convertalis(query_db, seg_db, seg_res, seg_tsv, threads=threads,
                       search_type=search_type, format_output=_SEGMENT_FORMAT)

    # Best V and best J per read. `JC|` targets are named by SCAFFOLD id, not by allele, so they
    # are resolved through the segment markup -- feeding the raw name into the combination lookup
    # silently fails for every J->C read and collapsed the fast path from 85.3 % to 0.1 % once.
    best_v: dict[str, str] = {}
    best_j: dict[str, str] = {}
    top: dict[tuple[str, str], float] = {}
    seg_rows: dict[tuple[str, str], dict] = {}
    # Reads whose best J-side evidence came from a `JC|` target, mapped to that J+C scaffold.
    # `JC|` targets are named by scaffold id, so the name IS the competing target.
    jc_scaffold: dict[str, str] = {}
    # Every V allele tied at the read's best segment score. V alleles of one gene differ by a
    # nucleotide or two, so a short read routinely cannot separate them and several tie exactly --
    # and picking one by segment score alone diverged from the one-pass V CALL on 38 % of amplicon
    # reads. Since the clonotype key is `(locus, v_call, j_call, junction)` at ALLELE level
    # (`rnaseq.correct`), that silently splits and merges clonotypes. All tied alleles' scaffolds
    # are aligned instead, so `_best_hits` applies the same rule the one-pass search does.
    tied_v: dict[str, list[str]] = {}
    for row in _segment_rows(seg_tsv):
        q, t, bits = row["query"], row["target"], float(row["bits"])
        kind, sep, name = t.partition("|")
        if not sep or kind not in ("V", "J", "JC"):
            continue                      # unrecognised target: never silently treat it as a J
        side = "V" if kind == "V" else "J"
        if (q, side) in top:              # rows arrive pre-sorted by (bits desc, target asc)
            # Keep the ties, drop everything below them. Capped: a read that ties against dozens
            # of alleles is degenerate, and the cap bounds the alignment work without affecting
            # the answer in any non-degenerate case.
            if side == "V" and bits == top[(q, side)] and len(tied_v[q]) < _MAX_TIED_V:
                tied_v[q].append(name)
            continue
        top[(q, side)] = bits
        seg_rows[(q, side)] = row
        allele = ref.segment_j_call(name) if kind == "JC" else name
        (best_v if side == "V" else best_j)[q] = allele
        if kind == "V":
            tied_v[q] = [name]
        if kind == "JC":
            jc_scaffold[q] = name

    if combos is None:
        combos = load_combinations(ref.target_fasta.parent / "combinations.tsv")
    sl = shortlist(best_v, best_j, combos)

    best, failed = ({}, set())
    if sl.implied:
        # Expand each read's single implied scaffold to every scaffold its exactly-tied V alleles
        # imply against the same best J. `dict.fromkeys` keeps order and de-duplicates; the
        # shortlist's own choice stays first so a read with no ties is unchanged.
        candidates = {
            q: list(dict.fromkeys(
                [sid] + [combos[(v, best_j[q])] for v in tied_v.get(q, ())
                         if (v, best_j[q]) in combos]))
            for q, sid in sl.implied.items()}
        best, failed = _align_implied(query_db, full_db, candidates, seg_rows, seqs, tmp,
                                      threads=threads, ref=ref)

    # A read that reaches its J and keeps going into C has TWO plausible homes: the V×J scaffold
    # its (V, J) pair names, and the J+C scaffold the segment pass actually hit. The one-pass rule
    # is "highest whole-scaffold bit score wins", so both must compete -- forcing the V×J choice
    # took a J->C read scoring 141 on a J+C scaffold, re-seated it on a V×J scaffold at 99, and
    # **fabricated a junction** from a weak spurious V hit (a 100 nt read always has *some* V above
    # threshold among 775). It also destroyed the `c_call`, i.e. the isotype. Measured on real
    # reads: 4 of 453 mapped reads, 3 of them gaining an invented junction_aa.
    # Two targets per read is still ~138x fewer alignments than the full reference.
    contest = {q: s for q, s in jc_scaffold.items() if q in sl.implied}
    if contest:
        jc_hits, _ = _align_implied(query_db, full_db, contest, seg_rows, seqs, tmp,
                                    threads=threads, side="J", ref=ref, tag_prefix="jc")
        for q, row in jc_hits.items():
            cur = best.get(q)
            if cur is None or float(row["bits"]) > float(cur["bits"]):
                best[q] = row
                failed.discard(q)         # the J+C alignment succeeded where V×J may not have

    rescue = sorted(set(sl.rescue) | failed)
    if rescue:
        best.update(_full_rescue(query_db, full_db, rescue, tmp, threads=threads,
                                 sensitivity=sensitivity, mm_strand=mm_strand,
                                 max_seqs=max_seqs, kmer=kmer, search_type=search_type))

    seen = set(best_v) | set(best_j)
    assert set(sl.implied) | set(sl.rescue) == seen, "a read was lost during shortlisting"
    n_fast = len(sl.implied) - len(failed & set(sl.implied))
    # Reads the segment pass never hit at all. Nearly all are genuinely non-receptor -- on bulk
    # that is ~97% of the library, which is the entire point -- but a few are real: a
    # junction-spanning read with weak homology can clear threshold on the concatenated
    # `V+pad+J` scaffold while neither half clears it alone, and raising sensitivity does not
    # recover them (measured to `-s 7.5 -e 1.0`). Counted here so the exposure is auditable
    # rather than implicit. Measured: 1 of 5,278 on amplicon at bits 56, 0 of 1,956 on bulk --
    # and **none at or above the default `--min-score 75`**, so nothing arda would report is lost.
    n_unseen = max(0, len(_db_keys(query_db)) - len(seen))
    report = {"implied": n_fast, "rescued": len(rescue), "no_segment_hit": n_unseen,
              "fast_fraction": round(n_fast / len(seen), 4) if seen else 0.0,
              "reasons": {**sl.reasons, **({"second_pass_failed": len(failed)} if failed else {})}}
    return best, report


def _extend_uncertain(best: dict[str, dict], records: list[tuple[str, str]], target_db: Path,
                      tmp: Path, *, threads: int, sensitivity: float, max_seqs: int,
                      kmer: int | None, seqtype: str,
                      strand: list[str]) -> tuple[dict[str, dict], int]:
    """Re-search, uncapped, the reads whose capped best score is low. Returns (hits, n_rechecked).

    A capped search can stop before reaching a read's true scaffold, and the reads that suffer are
    exactly the low-scoring ones (see `_ADAPTIVE_TRIGGER`). Re-running those without a cap is
    authoritative -- it is precisely the search arda does today -- so its answer replaces the
    capped one unconditionally rather than being compared to it.

    The uncertain set is small (measured 0.5 % of a bulk library), so this costs a few percent of
    the saving. If it is ever NOT small the adaptive path degrades to two full searches, which is
    slower than one; that is the failure mode to watch if the trigger is ever raised.
    """
    uncertain = {q for q, v in best.items() if float(v["bits"]) < _ADAPTIVE_TRIGGER}
    if not uncertain:
        return best, 0
    sub = [(q, s) for q, s in records if q in uncertain]
    if not sub:
        return best, 0
    d = tmp / "adaptive"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    fa = seqio.write_fasta(iter(sub), d / "q.fasta")
    qdb = d / "qDB"
    mmseqs.createdb(fa, qdb, dbtype=2 if seqtype == "nt" else 1)
    mmseqs.search(qdb, target_db, d / "res", d / "mt", search_type=_SEARCH_TYPE[seqtype],
                  sensitivity=sensitivity, max_seqs=max_seqs, threads=threads, kmer=kmer,
                  extra=strand or None)
    out_tsv = d / "hits.tsv"
    mmseqs.convertalis(qdb, target_db, mmseqs.top_hit(d / "res", d / "best"), out_tsv,
                       threads=threads, search_type=_SEARCH_TYPE[seqtype])
    best = dict(best)
    best.update(_best_hits(out_tsv))
    shutil.rmtree(d, ignore_errors=True)
    return best, len(sub)


def _merge_segment_report(acc: dict, chunk: dict) -> None:
    """Accumulate per-chunk shortlist counters into the run report.

    ``fast_fraction`` is recomputed from the totals rather than averaged -- averaging per-chunk
    fractions weights a 12-read tail chunk the same as a 400 k one.
    """
    for k in ("implied", "rescued", "no_segment_hit"):
        acc[k] = acc.get(k, 0) + chunk.get(k, 0)
    reasons = acc.setdefault("reasons", {})
    for k, v in chunk.get("reasons", {}).items():
        reasons[k] = reasons.get(k, 0) + v
    seen = acc["implied"] + acc["rescued"]
    acc["fast_fraction"] = round(acc["implied"] / seen, 4) if seen else 0.0


def _annotate_chunk(
    records: list[tuple[str, str]],
    ref: Reference,
    target_db: Path,
    seqtype: str,
    *,
    threads: int,
    sensitivity: float,
    mm_strand: int | None,
    map_d: bool = True,
    mapped_only: bool = False,
    max_seqs: int = _MAX_SEQS,
    kmer: int | None = -1,
    segment_db: Path | None = None,
    combos: dict[tuple[str, str], str] | None = None,
    report: dict | None = None,
    adaptive: bool = False,
) -> list[dict]:
    """Annotate one batch against a preloaded reference + cached target DB.

    ``mapped_only`` skips the empty record for non-hits (the RNA-seq filter path,
    where 95-99 % of reads have no hit and building throwaway records dominates).

    ``segment_db`` opts into the two-pass segment search (:func:`_segment_best_hits`): search the
    1,244-target segment reference, then align only the one V×J scaffold each read's V+J pair
    implies. Same ``{query: hit}`` shape, so nothing below this line changes. ``None`` (the
    default) keeps the one-pass search. ``combos`` is ``combinations.tsv`` preloaded by the
    caller -- it is 550 KB and re-parsing it per chunk is pure waste. Counters are accumulated
    into ``report`` when given.
    """
    if not records:
        return []
    if kmer == -1:  # sentinel: caller did not override, use the seqtype default
        kmer = _KMER[seqtype]
    with tempfile.TemporaryDirectory(prefix="arda_") as td:
        tmp = Path(td)
        query_fa = seqio.write_fasta(iter(records), tmp / "query.fasta")
        query_db = tmp / "queryDB"
        res_db = tmp / "resDB"
        out_tsv = tmp / "hits.tsv"
        mmseqs.createdb(query_fa, query_db, dbtype=2 if seqtype == "nt" else 1)
        if segment_db is not None and seqtype == "nt":
            best, seg_report = _segment_best_hits(
                query_db, segment_db, target_db, tmp, ref, threads=threads,
                sensitivity=sensitivity, mm_strand=mm_strand, max_seqs=max_seqs, kmer=kmer,
                search_type=_SEARCH_TYPE[seqtype], seqs=dict(records), combos=combos)
            if report is not None:
                _merge_segment_report(report, seg_report)
        else:
            strand = ["--strand", str(mm_strand)] if mm_strand is not None else []
            capped = ["--max-accept", str(_MAX_ACCEPT)] if adaptive else []
            mmseqs.search(
                query_db, target_db, res_db, tmp / "mmseqs_tmp",
                search_type=_SEARCH_TYPE[seqtype], sensitivity=sensitivity,
                max_seqs=max_seqs, threads=threads, kmer=kmer,
                extra=(strand + capped) or None,
            )
            # Reduce to one alignment per query BEFORE materialising it as text. With --max-seqs 300
            # a 100 k-read chunk emits ~804 k alignment rows -- 194 MB of cigar/qaln/taln -- of which
            # we keep 4 k. Parsing the other 800 k was arda's single largest memory consumer (877 MB
            # peak, vs 284 MB for mmseqs itself). Bit-identical: same target and score on all 4,101
            # queries.
            mmseqs.convertalis(query_db, target_db, mmseqs.top_hit(res_db, tmp / "bestDB"),
                               out_tsv, threads=threads, search_type=_SEARCH_TYPE[seqtype])
            best = _best_hits(out_tsv)
            if adaptive:
                best, n_re = _extend_uncertain(
                    best, records, target_db, tmp, threads=threads, sensitivity=sensitivity,
                    max_seqs=max_seqs, kmer=kmer, seqtype=seqtype, strand=strand)
                if report is not None:
                    report["adaptive_rechecked"] = report.get("adaptive_rechecked", 0) + n_re

    out: list[dict] = []
    for qid, qseq in records:
        hit = best.get(qid)
        entry = ref.get(hit["target"]) if hit else None
        if hit is None or entry is None:
            if mapped_only:
                continue
            rec = {c: "" for c in AIRR_COLUMNS}
            rec["sequence_id"], rec["sequence"] = qid, qseq
            out.append(rec)
            continue
        # mmseqs reports reverse-strand nt hits with qstart > qend and aligned
        # strings already on the coding strand. Re-orient: work on the revcomp,
        # remap the alignment start to forward coords on it. All markup, coords and
        # CIGARs are computed on ``work`` (the coding strand); the AIRR ``sequence``
        # field keeps the read AS SUBMITTED (``qseq``), with rev_comp=T signalling that
        # the output data are on its reverse complement -- per the AIRR spec.
        qs, qe = int(hit["qstart"]), int(hit["qend"])
        rev = qs > qe
        work = qseq
        if rev:
            work = reverse_complement(qseq)
            qlen = len(qseq)
            hit = dict(hit)
            hit["qstart"], hit["qend"] = qlen - qs + 1, qlen - qe + 1
        dg = ref.d_germlines.get(entry.locus) if map_d else None
        out.append(transfer_hit(qid, work, hit, entry, seqtype, rev_comp=rev,
                                d_germlines=dg, submitted_seq=qseq, anchors=ref.anchors))
    return out


def _prep(organism, seqtype, threads, sensitivity, strand):
    if seqtype not in _SEARCH_TYPE:
        raise ValueError(f"seqtype must be 'nt' or 'aa', got {seqtype!r}")
    # Reject an unknown strand rather than silently searching forward-only. On stranded paired
    # libraries R2 is antisense, so a typo here quietly discards ~40 % of the recoverable
    # repertoire (the R2-only fragments; arda-benchmark RESULTS.md §13d).
    if strand not in ("both", "forward"):
        raise ValueError(f"strand must be 'both' or 'forward', got {strand!r}")
    ref = load_reference(organism, seqtype)
    threads = threads or (os.cpu_count() or 1)
    sensitivity = _SENSITIVITY[seqtype] if sensitivity is None else sensitivity
    mm_strand = (2 if strand == "both" else 1) if seqtype == "nt" else None
    target_db = _cached_target_db(ref.target_fasta, organism, seqtype)
    return ref, target_db, threads, sensitivity, mm_strand


def annotate_records(
    records: list[tuple[str, str]],
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    sensitivity: float | None = None,
    strand: str = "both",
    map_d: bool = True,
) -> list[dict]:
    """Annotate in-memory ``(id, sequence)`` records; return AIRR record dicts.

    Args:
        strand: ``"both"`` (default, nt only) searches both strands and re-orients
            reverse-complement hits; ``"forward"`` searches the plus strand only.
            Ignored for protein input.
        map_d: ``True`` (default) maps D segments into the junction of VDJ-locus
            hits (``d_call``/``d2_call``/``np*``); ``False`` skips D mapping (nt
            input only — D mapping never runs for protein input).
    """
    ref, target_db, threads, sensitivity, mm_strand = _prep(
        organism, seqtype, threads, sensitivity, strand)
    return _annotate_chunk(records, ref, target_db, seqtype,
                           threads=threads, sensitivity=sensitivity,
                           mm_strand=mm_strand, map_d=map_d)


def annotate_file(
    input: str | Path,
    output: str | Path,
    organism: str = "human",
    seqtype: str = "nt",
    *,
    threads: int = 0,
    sensitivity: float | None = None,
    strand: str = "both",
    chunk_size: int = _CHUNK_SIZE,
    map_d: bool = True,
) -> Path:
    """Annotate a FASTA/FASTQ file and stream an AIRR TSV.

    The input is processed in bounded chunks with a background reader thread that
    prefetches the next chunk while the current one is annotated (mmseqs releases
    the GIL during its subprocess), so memory stays flat for arbitrarily large
    FASTQ and read parsing overlaps compute. The reference + target DB are loaded
    once and reused across all chunks.
    """
    output = Path(output)
    ref, target_db, threads, sensitivity, mm_strand = _prep(
        organism, seqtype, threads, sensitivity, strand)

    chunks: queue.Queue = queue.Queue(maxsize=2)
    # The reader runs in a daemon thread; anything it raises -- a missing file, a truncated gzip --
    # would otherwise die there unheard while `finally` still posts the sentinel, so the main loop
    # sees a clean end-of-stream and writes partial output with exit 0 (silent truncation). Capture
    # it and re-raise on the consumer side. Same guard as `rnaseq.map.map_rnaseq`.
    reader_exc: list[BaseException] = []

    def reader():
        try:
            for chunk in seqio.chunked(seqio.read_sequences(input), chunk_size):
                chunks.put(chunk)
        except BaseException as exc:  # noqa: BLE001 — re-raised in the consumer
            reader_exc.append(exc)
        finally:
            chunks.put(None)  # sentinel

    t = threading.Thread(target=reader, daemon=True)
    t.start()
    with open(output, "w") as fh:
        fh.write(airr_header() + "\n")
        while True:
            chunk = chunks.get()
            if chunk is None:
                if reader_exc:
                    raise reader_exc[0]
                break
            recs = _annotate_chunk(chunk, ref, target_db, seqtype,
                                   threads=threads, sensitivity=sensitivity,
                                   mm_strand=mm_strand, map_d=map_d)
            fh.write(format_rows(recs))
    t.join()
    return output
