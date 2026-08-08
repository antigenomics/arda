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
from .reference import load_reference, Reference, REGIONS
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

# The same problem exists on the J side and had no handling at all: J alleles of one gene are short
# and differ by a base or two, so a read routinely ties EXACTLY between two `J|` targets, and the
# tie was then broken lexicographically -- an arbitrary rule that decides which V×J scaffold the
# read is aligned against. Measured on the real-read fixture: SRR5233639.3589/2 ties
# `J|IGLJ2*01,IGLJ3*01` and `J|IGLJ2A*01` at **54 bits each**; the comma sorts before `A`, so the
# composite won and the read was seated on a scaffold scoring **93** while its true home scores
# **96**. Offering both lets `_best_hits` decide on whole-scaffold bit score, exactly as the
# one-pass search does and exactly as `_MAX_TIED_V` already does for V.
#
# Smaller than the V cap because J targets are ~40-60 nt: past a handful of exact ties the read
# carries no information that could separate them.
_MAX_TIED_J = 4

#: side -> how many exactly-tied rows may be kept for it. ⛔ ONE mapping, read by both
#: `_segment_rows` (polars) and `_segment_best_hits` (Python) -- see `_SEGMENT_SIDE` for what
#: happens when a rule about segment rows is spelled out twice. `C` is absent: a constant-region
#: hit nominates nothing that could be tied.
_MAX_TIED = {"V": _MAX_TIED_V, "J": _MAX_TIED_J}

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


def _has_jc_targets(fasta: Path) -> bool:
    """Was this ``segments.fasta`` written before the constant region became its own target?

    Scans headers only, and stops at the first ``>JC|``. The whole file is ~250 KB, so this costs
    nothing next to the mmseqs build it guards.
    """
    with open(fasta) as fh:
        return any(line.startswith(">JC|") for line in fh)


def _regenerate_segments(fasta: Path, organism: str) -> None:
    """Rewrite a pre-2.8.0 ``segments.fasta`` (+ markup), safely under concurrency.

    ``build_segment_reference`` truncates both artifacts **in place**. That was fine while it only
    ran from ``build-index``; 2.8.0 put it on the map path, where arda is concurrent by design
    (Nextflow process-per-sample, SLURM task-per-shard). N array tasks would all see ``JC|``, all
    truncate the same file, and ``_createdb_atomic``'s **mtime** gate would then happily compile an
    mmseqs DB from whatever bytes were on disk -- a short or interleaved reference, exit 0, targets
    quietly missing, every affected read pushed to rescue forever. Exactly the class of
    ``mmseqs createdb`` (0-byte ``db``) and ``fetch_database`` (cross-filesystem ``shutil.move``).

    So: take the same build lock the DB build takes, and re-check the format after acquiring it --
    the common case is not contention over the work, it is that the process we queued behind
    already did it.

    A read-only reference tree (container image, HPC module) is not an error either: the two-pass
    is simply unavailable, which is the contract ``_cached_segment_db`` already documents.
    """
    from .._locking import build_lock
    from ..refbuild.segments import build_segment_reference

    lock = fasta.parent / ".segments.build.lock"
    try:
        with build_lock(lock, done=lambda: not _has_jc_targets(fasta)) as mine:
            if not mine:
                return
            logger.info("segments.fasta predates the J+C collapse; regenerating for %s", organism)
            build_segment_reference(organism, out_dir=fasta.parent)
    except OSError as exc:
        logger.warning("could not regenerate segments.fasta for %s (%s); --two-pass will use the "
                       "pre-2.8.0 reference, which is correct but ~1.9x slower", organism, exc)


def _cached_segment_db(ref: Reference, organism: str) -> Path | None:
    """The mmseqs DB for the segment reference, or ``None`` if it has not been built.

    Nothing precompiled ships for this one -- `segments.fasta` is generated by ``build-index``
    (see :mod:`arda.refbuild.segments`), so this always goes through the local cache. Cached
    per *chunk-independent* key so a 400 k-read run builds it once, not once per chunk;
    `_createdb_atomic` holds the build lock, which is what makes that safe for the Nextflow
    process-per-sample and SLURM task-per-shard layouts.

    Returns ``None`` rather than raising when `segments.fasta` is absent: the two-pass is then
    simply unavailable and the caller falls back to the one-pass search.

    **A pre-2.8.0 `segments.fasta` is regenerated, not used.** Upgrading arda does not rewrite a
    generated artifact, and this one changed shape in 2.8.0: the 345 `JC|` scaffolds became 25 `C|`
    targets. The mapper still reads `JC|` (so a mixed-vintage install is correct, not broken),
    which is exactly what makes the stale case invisible -- an upgraded user passing ``--two-pass``
    would get correct output, no error, and none of the 1.89x, forever. Detecting it by FORMAT
    rather than by mtime or version is what makes that self-healing: the marker is the thing that
    actually changed.
    """
    fasta = ref.target_fasta.parent / "segments.fasta"
    if not fasta.exists():
        return None
    if _has_jc_targets(fasta):
        _regenerate_segments(fasta, organism)
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
    #
    # A TARGET-INVERTED row (`tstart > tend`) is the third shape of unusable, and the most
    # dangerous, because it does not raise -- it silently produces a well-formed WRONG junction.
    # arda detects a reverse-strand nt hit only from the QUERY side (`rev = qs > qe`, ~:1280), so
    # when mmseqs expresses the minus strand on the TARGET instead, the row looks forward.
    # `_markup.transfer_regions` then walks the target strictly forward from `tstart`
    # (`_markup/markup.cpp:185-201`), sliding the whole scaffold markup by
    # `(tlen + 1 - tstart) - tstart` nt and moving the junction window off Cys104 onto whatever
    # codon lands there. Measured on a delivered Jurkat run (arda 2.5.6, ERR3003543): tlen 349,
    # true tstart 170, reported tstart 180 -> the window slid exactly 10 nt into V framework 3,
    # started on a spurious TGT, ended on TGG, passed `assemble._CANON`'s `^C...[FW]$` and became
    # a 7,408-read phantom clonotype in a MONOCLONAL cell line -- 48 % of the true clone, from
    # which it had stolen 5,758 reads (ablation: 15,380 -> 21,138).
    #
    # These are not recoverable minus-strand hits to be reflected into forward coordinates. They
    # are internally inconsistent: on that read `germline_alignment` is the reverse complement of
    # the scaffold while the query matches the PLUS strand at 91.4 % identity (the row's own
    # reported `pident`), and `identity(qaln, taln)` is 0.232 -- which is why `v_identity` came out
    # 0.216 against ~0.98 for every normal record. Reflecting the coordinates would keep a garbage
    # alignment; dropping routes the read to the full-reference rescue, or leaves it unmapped, and
    # that is the only option that cannot ship a well-formed junction that is wrong.
    #
    # Emission is mmseqs-BUILD dependent -- 120 such rows across six delivered samples, 0 on the
    # build in the local env at the same arda version -- so this is a robustness gate, not a
    # regression fix, and comparing two arda versions on one machine cannot catch it.
    n_before = df.height
    df = df.filter(pl.col("query").is_not_null() & (pl.col("query") != "")
                   & pl.col("bits").is_not_null()
                   & (pl.col("tstart") <= pl.col("tend")))
    if df.height != n_before:
        logger.warning("%d alignment rows were unusable (no query id, no score, or an inverted "
                       "target span) and were routed to the full-reference rescue; a large count "
                       "means the hand-built sub-DB is missing lookup entries",
                       n_before - df.height)
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

#: Segment target prefix -> which side of a scaffold that row is evidence for. Anything not listed
#: is dropped and is NEVER silently treated as a J.
#:
#: ⛔ ONE mapping, used by both `_segment_rows` (which reduces the alignment TSV in polars) and the
#: loop in `_segment_best_hits` (which consumes it). They are the same rule, and when they were
#: written out separately, adding `C|` targets to the reference and updating only the loop made the
#: reduction discard every C row: `best_c` stayed empty, no constant-only read was ever rescued,
#: and 15 J->C reads vanished **without `no_segment_hit` moving**, because the rows were gone
#: before any counter saw them. Sharing the mapping makes that divergence unrepresentable.
#:
#: `C` is its OWN side, not the J side: a constant-region hit says what the isotype is and nothing
#: at all about which J the read carries. `JC` is the pre-split kind, kept on the J side so this
#: mapper still works against a reference built before the constant region became its own target.
_SEGMENT_SIDE = {"V": "V", "J": "J", "JC": "J", "C": "C"}


def _alleles(call: str) -> list[str]:
    """Split a possibly comma-joined call into individual alleles.

    A `v_call`/`j_call`/`c_call` names a GROUP of alleles arda could not tell apart, and different
    parts of the reference build group them by different rules. Anything matching one call against
    another must compare alleles, not the joined strings.
    """
    return [a for a in (s.strip() for s in call.split(",")) if a]


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
    # ⛔ The SAME unusable-row filter `_best_hits` applies, for a sharper reason: polars sorts
    # nulls FIRST under `descending=True`, so a row with an empty `bits` field becomes `_rank == 0`
    # for its `(query, side)` and **evicts the read's real best hit**. Reproduced: with rows
    # `V|A*01 <empty bits>` and `V|B*01 120` for one read, the reduction returns `V|A*01` and
    # `V|B*01` is gone. The null then reaches `float(row["bits"])` in `_segment_best_hits` and
    # raises mid-chunk -- after earlier chunks were already written, i.e. a partial AIRR file that
    # looks complete. An empty query id is the same malformed row and would key every downstream
    # dict on `None`. Both route the read to the full-reference rescue instead, which is the
    # guarantee `--two-pass` is built on.
    n_before = df.height
    df = df.filter(pl.col("query").is_not_null() & (pl.col("query") != "")
                   & pl.col("bits").is_not_null())
    if df.height != n_before:
        logger.warning("%d segment alignment rows were unusable (no query id or no score) and were "
                       "dropped before best-hit selection", n_before - df.height)
    # Deterministic: exact bit-score ties between paralogous targets are broken on `target`,
    # the same rule `_best_hits` uses. Resolving them by TSV row order (as an earlier draft did)
    # silently undid the determinism fix -- 25/25 tied queries flipped when rows were reversed.
    df = df.sort(["bits", "target"], descending=[True, False], maintain_order=True)

    # Which side of a scaffold each row is evidence for, from the SHARED `_SEGMENT_SIDE` mapping --
    # see its comment for why this must not be spelled out separately here.
    kind = pl.col("target").str.split("|").list.first()
    df = df.filter(kind.is_in(list(_SEGMENT_SIDE))).with_columns(
        _side=kind.replace_strict(_SEGMENT_SIDE))

    # ⛔ `over`, NOT `group_by`. A window function maps its result back to the ORIGINAL row
    # positions, so on an already-sorted frame it is deterministic; `group_by` is a multithreaded
    # hash aggregation, and using it unordered is precisely what made `correct` nondeterministic
    # across runs while the row count stayed stable.
    df = df.with_columns(
        _rank=pl.int_range(pl.len()).over(["query", "_side"]),
        _top=pl.col("bits").first().over(["query", "_side"]))
    # Top row per (query, side), plus exactly-tied rows up to that side's cap -- from the shared
    # `_MAX_TIED`, so the reduction and the consuming loop cannot disagree about it.
    keep = (pl.col("_rank") == 0) | (
        (pl.col("bits") == pl.col("_top"))
        & (pl.col("_rank") < pl.col("_side").replace_strict(_MAX_TIED, default=1)))
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


#: Header of the ``ARDA_VONLY_DUMP`` calibration table. Order is the contract the offline fit reads.
_VONLY_COLS = ("read_id", "v_allele", "seg_bits", "seg_qstart", "seg_qend", "seg_tstart",
               "read_len", "scaffold", "scaffold_bits")

#: Header of the ``ARDA_PROJECT_DUMP`` validation table.
_PROJECT_COLS = ("read_id", "locus", "v_call", "j_call", "refusal",
                 "proj_junction", "proj_start", "proj_end", "rev_comp")


def _dump_projection(best_v: dict[str, str], best_j: dict[str, str],
                     seg_rows: dict[tuple[str, str], dict], seqs: dict[str, str],
                     ref: Reference, *, split_checked: bool) -> None:
    """Append the ARITHMETIC junction for every read carrying both anchors, or why it was refused.

    Validation only, and deliberately **out of the output path**: the projection is not yet what
    arda reports, so this writes it beside the scaffold-derived answer and lets an offline script
    join the two on ``sequence_id``. Wiring it in first and comparing second would make a
    disagreement look like a regression instead of a measurement.

    The refusal counts are as interesting as the junctions. Fast-path yield is 87.0 % of hit
    TRA-amplicon reads against 7.2 % of bulk reads, because bulk reads mostly do not span a junction
    at all -- and a fast path that silently covered 7 % of a library while looking like a speedup is
    exactly the "correct output, zero speedup" failure this project has shipped twice.

    No-op unless ``ARDA_PROJECT_DUMP`` names a path. Appends, because ``map`` runs per chunk.
    """
    path = os.environ.get("ARDA_PROJECT_DUMP")
    if not path:
        return
    from ..refbuild.translate import reverse_complement
    from .project import _anchor, project_junction

    p = Path(path)
    rows = []
    for q in sorted(set(best_v) & set(best_j)):
        v_row, j_row = seg_rows.get((q, "V")), seg_rows.get((q, "J"))
        seq = seqs.get(q)
        if not v_row or not j_row or not seq:
            continue
        # `project_junction` works in the frame the hits were measured in, so a minus-strand read is
        # handed its reverse complement. Reflecting coordinates afterwards instead is how sign
        # errors get in -- see the module docstring.
        rc = v_row["qstart"] > v_row["qend"]
        strand_seq = reverse_complement(seq) if rc else seq
        proj, why = project_junction(strand_seq, len(seq), v_row=v_row, j_row=j_row,
                                     v_call=best_v[q], j_call=best_j[q], anchors=ref.anchors,
                                     split_checked=split_checked)
        # Locus from the J anchor, not the V: TRAV/DV alleles pair with either TRAJ or TRDJ and
        # **the J decides the locus**. Taking it from the V would mislabel every TRD read.
        ja = _anchor(ref.anchors, "J", best_j[q])
        rows.append("\t".join(str(x) for x in (
            q, ja.locus if ja else "", best_v[q], best_j[q],
            why, proj.junction if proj else "", proj.start if proj else "",
            proj.end if proj else "", int(rc))))
    if not rows:
        return
    new = not p.exists()
    with p.open("a") as fh:
        if new:
            fh.write("\t".join(_PROJECT_COLS) + "\n")
        fh.write("\n".join(rows) + "\n")


def _dump_vonly(rescue: list[str], best_v: dict[str, str], best_j: dict[str, str],
                seg_rows: dict[tuple[str, str], dict], best: dict[str, dict],
                seqs: dict[str, str]) -> None:
    """Append one row per ``v_only`` read: its SEGMENT score beside its SCAFFOLD score.

    A ``v_only`` read hit a V and never reached a J, because on a 5'RACE amplicon **there is no J
    in the read**. It is 43 % of amplicon mates and 93.4 % of the rescue set, and the only reason
    it goes to the full 15,414-scaffold reference at all is to obtain a score on the scale
    ``--min-score`` is defined in. Skipping that search is worth ~2.5x (measured, round 14 step 1),
    but only if a segment-scale threshold can reproduce the same kept set.

    This writes the raw pair of scores and nothing else. It deliberately does **not** fit or apply
    a threshold: the fit belongs offline, where it can be cross-validated per locus against the
    round-12/13 truth, and where getting it wrong cannot silently change a shipped call. Round 6
    measured ``--min-score 60`` taking precision 94.3 % -> 65.5 %, which is the cost of guessing.

    ⚠ The two scores are **not** on one scale and are not meant to be: ``seg_bits`` is an ungapped
    ``MATCH 2 / MISMATCH -3`` score (``src/_segmap/segmap.cpp``) that grows with aligned length,
    while ``scaffold_bits`` is an MMseqs2 bit score. Finding the map between them is the experiment.
    ``seg_qstart``/``seg_qend``/``read_len`` are emitted so the fit can normalise by aligned length
    rather than assuming a single global cut -- V genes differ in length across loci, so a raw
    ungapped score is not comparable between them.

    No-op unless ``ARDA_VONLY_DUMP`` names a path. Appends, because ``map`` runs per chunk.
    """
    path = os.environ.get("ARDA_VONLY_DUMP")
    if not path:
        return
    p = Path(path)
    rows = []
    for q in rescue:
        # The `v_only` predicate, restated from `shortlist()`: a V but no J. `Shortlist.reasons`
        # only carries counts, so recomputing here is what keeps this out of the shipped partition.
        if not best_v.get(q) or best_j.get(q):
            continue
        seg = seg_rows.get((q, "V")) or {}
        hit = best.get(q) or {}
        rows.append("\t".join(str(x) for x in (
            q, best_v[q], seg.get("bits", ""), seg.get("qstart", ""), seg.get("qend", ""),
            seg.get("tstart", ""), len(seqs.get(q, "")),
            hit.get("target", ""), hit.get("bits", ""))))
    if not rows:
        return
    new = not p.exists()
    with p.open("a") as fh:
        if new:
            fh.write("\t".join(_VONLY_COLS) + "\n")
        fh.write("\n".join(rows) + "\n")


def _segment_best_hits(
    query_db: Path, seg_db: Path, full_db: Path, tmp: Path, ref: Reference, *,
    threads: int, sensitivity: float, mm_strand: int | None, max_seqs: int, kmer: int | None,
    search_type: int, seqs: dict[str, str],
    combos: dict[tuple[str, str], str] | None = None,
    fast_segments: bool = False,
    indel_rescue: bool = False,
    segment_only_v: bool = False,
) -> tuple[dict[str, dict], dict]:
    """Two-pass best hits: cheap segment pass, then align only the implied scaffold.

    Produces the same ``{query: hit_row}`` shape as the one-pass path, so everything downstream
    (strand handling, markup transfer, D mapping) is untouched.

    **No read is lost.** Whatever the segment pass cannot resolve to a single V×J scaffold --
    V-only, J-only, a pair absent from the reference, a reverse-strand hit, or a second-pass
    miss -- is realigned against the FULL reference exactly as today. The partition is asserted,
    not assumed, and `_full_rescue` raises rather than returning empty.
    """
    from .shortlist import _lookup, load_combinations, shortlist

    if fast_segments:
        # Structure-aware path: seed, vote by diagonal, extend ungapped. 37x the search below on
        # the same reads and reference, agreeing with it on .9997 of V alleles and .9998 of J.
        # It only NOMINATES: every candidate is still aligned against the full scaffold and scored
        # by MMseqs2 below, which is why a different score scale here is not a correctness problem.
        from .. import segmap
        # The same file `_cached_segment_db` compiles for MMseqs2, read directly. Derived from
        # `ref`, not from `seg_db`, so there is one definition of where the segment reference lives.
        rows = segmap.segment_rows(ref.target_fasta.parent / "segments.fasta",
                                   seqs, max_tied=max(_MAX_TIED.values()), threads=threads,
                                   max_indel=segmap.MAX_INDEL_NT if indel_rescue else 0)
    else:
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
        rows = _segment_rows(seg_tsv)

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
    # ...and the same on the J side, which had no handling at all. See `_MAX_TIED_J`: an exact tie
    # between two `J|` targets was being broken lexicographically, which decided the scaffold.
    tied_j: dict[str, list[str]] = {}
    # Best constant-region hit per read, from the `C|` targets. Its own side, NOT the J side: a C
    # hit is evidence about the isotype and none whatsoever about which J the read carries, and
    # folding it into `best_j` is what the old `JC|` targets did -- they won the J side on the
    # strength of their constant half.
    best_c: dict[str, str] = {}
    # Reads whose best segment evidence sits on two diagonals of one target -- the signature of an
    # indel, which a single ungapped extension scores only up to. Such a read's segment score is
    # systematically low, so letting it take the fast path decides its scaffold on truncated
    # evidence. They are demoted to `rescue` below and realigned GAPPED against the full reference.
    split_reads: set[str] = set()
    tied = {"V": tied_v, "J": tied_j}
    for row in rows:
        q, t, bits = row["query"], row["target"], float(row["bits"])
        kind, sep, name = t.partition("|")
        side = _SEGMENT_SIDE.get(kind) if sep else None
        if side is None:
            continue                      # unrecognised target: never silently treat it as a J
        if row.get("split"):
            split_reads.add(q)
        if (q, side) in top:              # rows arrive pre-sorted by (bits desc, target asc)
            # Keep the ties, drop everything below them. Capped: a read that ties against dozens
            # of alleles is degenerate, and the cap bounds the alignment work without affecting
            # the answer in any non-degenerate case.
            bucket = tied.get(side)
            if bucket is not None and bits == top[(q, side)] and len(bucket[q]) < _MAX_TIED[side]:
                bucket[q].append(name)
            continue
        top[(q, side)] = bits
        seg_rows[(q, side)] = row
        allele = ref.segment_j_call(name) if kind == "JC" else name
        {"V": best_v, "J": best_j, "C": best_c}[side][q] = allele
        if side in tied:
            # A `JC|` target is J-side but is named by scaffold id, so it can never be a J ALLELE
            # candidate for `combinations.tsv`; seed the bucket empty rather than with its name.
            tied[side][q] = [] if kind == "JC" else [name]
        if kind == "JC":
            jc_scaffold[q] = name

    # A read with both a J and a C hit names its J+C scaffold by that pair, exactly as a V→J read
    # names its V×J scaffold through `combinations.tsv`. `setdefault` so a pre-split `JC|` hit,
    # which already knows its scaffold, is never overwritten.
    #
    # ⛔ Both calls are expanded to individual ALLELES before the lookup. A `j_call` is a group of
    # alleles arda cannot separate, and the two sides group them by different rules -- see
    # `Reference.jc_combinations`. Matching the comma-joined strings leaves 24 human J+C scaffolds
    # unreachable, including every IGLJ2/IGLJ3 read, which turns the contest off for exactly the
    # reads it exists to protect.
    if best_c:
        jc_combos = ref.jc_combinations()
        for q, c in best_c.items():
            sid = next((s for j in _alleles(best_j.get(q, "")) for ca in _alleles(c)
                        if (s := jc_combos.get((j, ca))) is not None), None)
            if sid:
                jc_scaffold.setdefault(q, sid)

    # Before the shortlist, so the dump sees every read carrying both anchors -- including those the
    # shortlist will send to rescue for a reason unrelated to the junction (an unknown V*J pair, an
    # indel). Dumping after would silently under-report fast-path yield.
    # `indel_rescue` is exactly the flag that makes `segment_rows` compute `split`, so it IS the
    # answer to "did the indel check run" -- passed explicitly rather than re-derived.
    _dump_projection(best_v, best_j, seg_rows, seqs, ref, split_checked=indel_rescue)

    if combos is None:
        combos = load_combinations(ref.target_fasta.parent / "combinations.tsv")
    sl = shortlist(best_v, best_j, combos)

    # Demote indel-bearing reads from the fast path to the rescue set. This is a REROUTE, never a
    # drop: `rescue` is realigned against the full reference by MMseqs2, which gaps, so these reads
    # get a better answer than the fast path could give them -- not a worse one. Doing it here,
    # after `shortlist` has asserted its partition is total, keeps that invariant intact and keeps
    # the whole mechanism to one place.
    if split_reads:
        demoted = [q for q in sl.implied if q in split_reads]
        for q in demoted:
            del sl.implied[q]
            sl.rescue.append(q)
        seg_report_split = len(demoted)
    else:
        seg_report_split = 0

    best, failed = ({}, set())
    if sl.implied:
        # Expand each read's single implied scaffold to every scaffold its exactly-tied V alleles
        # imply against the same best J. `dict.fromkeys` keeps order and de-duplicates; the
        # shortlist's own choice stays first so a read with no ties is unchanged.
        #
        # Resolved through `shortlist._lookup`, not a bare `combos[...]`: a tied V can carry a
        # composite (comma-joined) allele name, which `combinations.tsv` only ever registers by
        # individual member -- 23 of 775 human `V|` targets do, `IGHV3-23*01,IGHV3-23D*01` among
        # them. A bare lookup drops exactly those siblings from the candidate set, so the read is
        # aligned against fewer scaffolds than it ties against and the allele call is decided by
        # which ones happened to survive.
        # Two AXES, not the cross product: the tied V alleles against the best J, and the tied J
        # alleles against the best V. Bounded by `_MAX_TIED_V + _MAX_TIED_J` rather than their
        # product, and a read tied simultaneously on both sides is vanishingly rare -- either axis
        # surfaces the scaffold that then wins on whole-scaffold bit score.
        def _cands(q: str, sid: str) -> list[str]:
            out = [sid]
            bj, bv = best_j[q], best_v[q]
            out += [s for v in tied_v.get(q, ()) if (s := _lookup(combos, v, bj)) is not None]
            out += [s for j in tied_j.get(q, ()) if (s := _lookup(combos, bv, j)) is not None]
            return list(dict.fromkeys(out))

        candidates = {q: _cands(q, sid) for q, sid in sl.implied.items()}
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
    #
    # ⛔ **Nominate from the J, not from a C hit.** Before the constant region became its own
    # target, a J->C read's best J-side hit WAS a `JC|` scaffold, so it named its own contestant.
    # It no longer does, and requiring a `C|` hit instead is not equivalent: such a read spans the
    # J/C boundary, so its constant overlap is often too short to clear the search threshold on its
    # own while the concatenated J+C scaffold it used to hit cleared it easily. Measured on the
    # real-read fixture, gating the contest on C evidence let exactly the bug this contest exists
    # to prevent back in on SRR5233639.12648/1 -- `TRBV12-3*02` invented, `c_call` TRBC2*01
    # destroyed, and `junction_aa` CASSFAGLVNIDEQFF fabricated on a read the one-pass calls V-less.
    #
    # So every implied read offers the J+C scaffolds of its best J, and bit score decides, exactly
    # as the one-pass does. That is 1 extra candidate for TRA/TRD/IGK, 2 for TRB/TRG, 7 for IGL and
    # 11 for IGH. The diagonal transfers cleanly because a J+C scaffold starts AT its J
    # (`j_sequence_start` = 1), which is the same frame the `J|` segment hit is measured in -- the
    # per-scaffold offset that makes this wrong on a V×J scaffold (`len(V) + pad`) does not exist
    # here.
    jc_by_j: dict[str, list[str]] = {}
    for (j_allele, _c), sid in ref.jc_combinations().items():
        jc_by_j.setdefault(j_allele, []).append(sid)
    contest: dict[str, str | list[str]] = {}
    for q in sl.implied:
        named = jc_scaffold.get(q)            # a pre-split `JC|` hit already knows its scaffold
        if named:
            contest[q] = [named]
            continue
        # Per allele, for the same reason as the lookup above; `dict.fromkeys` de-duplicates while
        # keeping order, so the candidate list is deterministic.
        cands = list(dict.fromkeys(
            s for j in _alleles(best_j.get(q, "")) for s in jc_by_j.get(j, ())))
        if cands:
            contest[q] = cands
    if contest:
        jc_hits, _ = _align_implied(query_db, full_db, contest, seg_rows, seqs, tmp,
                                    threads=threads, side="J", ref=ref, tag_prefix="jc")
        for q, row in jc_hits.items():
            cur = best.get(q)
            if cur is None:
                # The V×J alignment produced NOTHING for this read, so there is no score to
                # compete against and this is not a contest -- it is a walkover. Taking the J+C
                # row here and dropping the read from `failed` would hand it a V-less answer (no
                # `v_call`, no junction, no clonotype) that the full reference was never asked
                # about, and the one-pass may well score its V×J scaffold higher. The read stays
                # in `failed`, so `_full_rescue` decides against the WHOLE reference -- which
                # contains these J+C scaffolds too, and will pick one if it really is the best
                # home. The "no read lost" invariant holds either way; what this protects is the
                # clonotype, which no counter in the report would have shown going missing.
                continue
            if float(row["bits"]) > float(cur["bits"]):
                best[q] = row
                failed.discard(q)         # a real contest: the J+C scaffold outscored the V×J one

    # A read whose ONLY segment evidence is a constant-region hit has no (V, J) pair, so the
    # shortlist never sees it. It is real receptor mRNA carrying no V(D)J -- rescue it against the
    # full reference rather than letting it fall out of `seen`, which is how a read gets lost with
    # no error: the partition assertion below only checks what the shortlist was told about.
    c_only = set(best_c) - set(best_v) - set(best_j)
    rescue = sorted(set(sl.rescue) | failed | c_only)

    # A `v_only` read carries NO J -- that is the definition of the class, not a search failure --
    # so searching it against 15,414 V×J scaffolds asks a question whose J half the read cannot
    # answer. It is 77 % of the amplicon rescue set and ~44 % of the amplicon wall, at 338 µs/read
    # against 31 µs for a named-target alignment (results/round18).
    #
    # ⛔ This is NOT round 5's refuted narrowing, which kept the read on SCAFFOLDS and narrowed
    # which ones (so a read whose true home scored higher elsewhere was trapped, and the narrowed
    # best then fell under `--min-score 75`). Here the read is aligned against its own V SEGMENT,
    # by MMseqs2, producing a real bit score over exactly the nucleotides a whole-scaffold
    # alignment of a J-less read would have covered -- so `--min-score` keeps its meaning and the
    # segmap-scale calibration this was blocked on is not needed at all.
    #
    # Anything that fails here stays in `rescue`: the "no read lost" invariant is unchanged.
    #
    # ⛔ The class is gated by GEOMETRY, not by the shortlist reason alone. `v_only` means "no J
    # segment hit", which on a read carrying SHM is not the same statement as "no J in the read":
    # the segment pass misses short hypermutated IGHJ and the full reference then finds it.
    # Routing on the reason code alone cost 77 of 213 bulk junctions.
    #
    # ⛔ And the gate must be measured from the READ's full extent, not from where the ALIGNMENT
    # stopped. An ungapped extension breaks early under mismatch load, so on a hypermutated
    # library the alignment end systematically understates how far the read reaches: gating on it
    # passed every local test (100 k TRA amplicon, 100 k bulk, both zero-loss) and then lost
    # **45 of 88,697 junctions on IGH_repertoire** (median 91.77 % V identity) at cluster scale.
    # Projecting the read's remaining 3' bases along the same diagonal is SHM-proof -- unaligned
    # bases still occupy target positions -- and it is still pure geometry: if the whole read,
    # laid on its diagonal, cannot reach Cys104, no junction can be in it.
    v_seg_hits: dict[str, dict] = {}
    if segment_only_v and seg_db is not None:
        def _cannot_reach_cys104(q: str) -> bool:
            row = seg_rows.get((q, "V"))
            if row is None:
                return False
            entry = ref.get(row["target"])
            cdr3_start = entry.starts[REGIONS.index("cdr3")] if entry else -1
            if cdr3_start <= 0:
                return False               # no CDR3 markup on this segment: do not risk it
            qs, qe = int(row["qstart"]), int(row["qend"])
            t_end = int(row["tstart"]) + abs(qe - qs)
            # Query bases past the alignment end, walking the target FORWARD. segmap emits
            # `qstart > qend` for a minus-strand hit, where target-forward is query-backward.
            tail = (len(seqs[q]) - qe) if qs < qe else (qe - 1)
            return t_end + tail < cdr3_start - 3   # the Cys codon is unreachable, not just unread

        v_only = {q: list(dict.fromkeys([f"V|{best_v[q]}"]
                                        + [f"V|{a}" for a in tied_v.get(q, ())]))
                  for q in rescue
                  if sl.reason_of.get(q) == "v_only" and q in best_v
                  and _cannot_reach_cys104(q)}
        if v_only:
            v_seg_hits, v_seg_failed = _align_implied(
                query_db, seg_db, v_only, seg_rows, seqs, tmp,
                threads=threads, side="V", ref=ref, tag_prefix="vseg")
            rescue = [q for q in rescue if q not in v_seg_hits or q in v_seg_failed]

    if rescue:
        best.update(_full_rescue(query_db, full_db, rescue, tmp, threads=threads,
                                 sensitivity=sensitivity, mm_strand=mm_strand,
                                 max_seqs=max_seqs, kmer=kmer, search_type=search_type))
        _dump_vonly(rescue, best_v, best_j, seg_rows, best, seqs)
    best.update(v_seg_hits)

    seen = set(best_v) | set(best_j) | set(best_c)
    # ⛔ Assert what actually matters: every read the segment pass SAW either came back with a hit
    # or was handed to the full-reference rescue. The obvious form of this check --
    # `implied | rescue | c_only == seen` -- is a TAUTOLOGY: `shortlist` already asserts
    # `implied | rescue == set(best_v) | set(best_j)` internally, and `c_only` is defined three
    # lines up as `set(best_c) - set(best_v) - set(best_j)`, so the two sides are equal by
    # construction and it can never fire. It was standing in for this one, which can.
    unaccounted = seen - set(best) - set(rescue)
    assert not unaccounted, (
        f"{len(unaccounted)} read(s) hit a segment target but were neither aligned nor rescued, "
        f"e.g. {sorted(unaccounted)[:3]}")
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
              "reasons": {**sl.reasons,
                          **({"second_pass_failed": len(failed)} if failed else {}),
                          **({"indel_rescued": seg_report_split} if seg_report_split else {}),
                          **({"c_only": len(c_only)} if c_only else {})},
              **({"v_only_on_segment": len(v_seg_hits)} if v_seg_hits else {})}
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
    fast_segments: bool = False,
    indel_rescue: bool = False,
    segment_only_v: bool = False,
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
                search_type=_SEARCH_TYPE[seqtype], seqs=dict(records), combos=combos,
                fast_segments=fast_segments, indel_rescue=indel_rescue,
                segment_only_v=segment_only_v)
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
