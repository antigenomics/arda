"""Thin wrapper around the ``mmseqs`` binary.

Inspired by pymmseqs (MIT) but deliberately dependency-free: we only need
binary discovery, a subprocess runner, and the ``createdb`` / ``search`` /
``convertalis`` (and ``easy-search``) pipeline used by the annotator.

Discovery order for the binary: ``$ARDA_MMSEQS`` → the optional ``arda-mmseqs``
companion wheel → ``<project>/bin/mmseqs`` → ``mmseqs`` on ``PATH``. Candidates
after the explicit override are **version-matched against the precompiled indexes
in** ``database/``: an index is only reusable by the mmseqs release that built it,
so accepting an arbitrary PATH binary silently discards the shipped index and
rebuilds a private cache. If nothing matches, a known-good static binary is
auto-fetched into ``<project>/bin/mmseqs`` (one-time, transparent) unless
``$ARDA_NO_AUTO_FETCH`` is set — so neither pip nor conda users need to install
mmseqs manually -- the binary is fetched on first use.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from functools import lru_cache

from .paths import bin_dir

__all__ = [
    "MMseqsError",
    "mmseqs_binary",
    "run",
    "version",
    "version_key",
    "versions_compatible",
    "createdb",
    "search",
    "convertalis",
    "top_hit",
    "easy_search",
    "DEFAULT_FORMAT_OUTPUT",
]

# Fields needed to transfer reference markup onto a query. 1-based inclusive
# coords for qstart/qend/tstart/tend; cigar + qaln/taln drive the projection.
DEFAULT_FORMAT_OUTPUT = (
    "query,target,qstart,qend,tstart,tend,qlen,tlen,"
    "alnlen,mismatch,gapopen,cigar,qaln,taln,evalue,bits,pident"
)

# search-type values (see `mmseqs search --help`).
SEARCH_TYPE_AUTO = 0
SEARCH_TYPE_PROTEIN = 1   # aa query vs aa target
SEARCH_TYPE_TRANSLATED = 2  # nt query vs aa target (blastx-like)
SEARCH_TYPE_NUCLEOTIDE = 3  # nt query vs nt target


class MMseqsError(RuntimeError):
    """Raised when an ``mmseqs`` invocation exits non-zero."""


def _bundled_binary() -> str | None:
    """The binary shipped by the optional ``arda-mmseqs`` companion wheel, if installed."""
    try:
        from arda_mmseqs import binary  # type: ignore[import-not-found]
    except Exception:  # noqa: BLE001 — not installed is the normal case
        return None
    try:
        p = Path(binary())
        return str(p) if p.exists() else None
    except Exception:  # noqa: BLE001
        return None


@lru_cache(maxsize=None)
def _version_of(path: str) -> str | None:
    """``mmseqs version`` for one specific binary, or None if it will not run.

    Separate from :func:`version` because that one resolves through
    :func:`mmseqs_binary`, and resolution needs to interrogate candidates first.
    """
    try:
        proc = subprocess.run([path, "version"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return proc.stdout.strip() or None if proc.returncode == 0 else None


def committed_index_version() -> str | None:
    """The mmseqs version the precompiled indexes in ``database/`` were built with.

    ``None`` when no index ships (the packaged reference asset deliberately omits them, so a
    plain ``pip install`` has nothing to match) — in that case any mmseqs will do, because the
    first run builds its own cache regardless.

    Deliberately does *not* go through :func:`~arda.paths.vdj_dir`: that resolves
    ``database_dir()``, which will auto-*download* the reference. Binary discovery must never
    trigger a network fetch.
    """
    try:
        from .paths import _source_root, cache_root

        roots = []
        src = _source_root()
        if src is not None:
            roots.append(src / "database")
        roots.append(cache_root() / "database")
        for root in roots:
            for ver in sorted((root / "vdj").glob("*/mmseqs/*/VERSION")):
                text = ver.read_text().strip()
                if text:
                    return text
    except Exception:  # noqa: BLE001 — discovery must not fail on a layout surprise
        return None
    return None


@lru_cache(maxsize=1)
def mmseqs_binary() -> str:
    """Locate an mmseqs executable that can actually use the shipped indexes.

    Resolution: ``$ARDA_MMSEQS`` → the ``arda-mmseqs`` companion wheel →
    ``<project>/bin/mmseqs`` → ``mmseqs`` on ``PATH`` → auto-fetched static build.

    **Version-matched, not merely present.** Taking whatever `mmseqs` happened to be on PATH
    was a silent correctness and performance bug: an index is only reusable by the release it
    was compiled with, so a mismatched binary makes every run reject ``database/``'s
    precompiled DBs and rebuild a private cache instead — no error, just a slow start and, if
    the two releases align differently, results that are not comparable with anyone else's.
    Found in the wild: a cluster with a bare-git-hash build ahead of conda's on PATH.

    ``$ARDA_MMSEQS`` is never version-checked — an explicit override is the user's call.
    If nothing matches, the known-good static build is fetched (unless ``$ARDA_NO_AUTO_FETCH``);
    if that also fails we fall back to the best candidate and warn, naming the consequence.
    """
    env = os.environ.get("ARDA_MMSEQS")
    if env:
        return env

    candidates = [c for c in (_bundled_binary(),
                              str(bin_dir() / "mmseqs") if (bin_dir() / "mmseqs").exists() else None,
                              shutil.which("mmseqs")) if c]

    want = committed_index_version()
    if want is None:
        # No shipped index to be compatible with; any working mmseqs is fine.
        if candidates:
            return candidates[0]
    else:
        for cand in candidates:
            got = _version_of(cand)
            if got and versions_compatible(got, want):
                return cand

    if "ARDA_NO_AUTO_FETCH" not in os.environ:
        fetched = _auto_fetch()
        if fetched is not None:
            got = _version_of(fetched)
            if want is None or (got and versions_compatible(got, want)):
                return fetched
            candidates.append(fetched)

    if candidates:
        import warnings

        got = _version_of(candidates[0]) or "unknown"
        warnings.warn(
            f"mmseqs {got} does not match the version the shipped indexes were built with "
            f"({want}); the precompiled reference index will be ignored and a private cache "
            f"rebuilt on first use. Let arda auto-fetch a matching build (unset "
            f"$ARDA_NO_AUTO_FETCH) or set $ARDA_MMSEQS to silence this.",
            RuntimeWarning, stacklevel=2,
        )
        return candidates[0]

    raise MMseqsError(
        "mmseqs binary not found. Allow auto-fetch (unset $ARDA_NO_AUTO_FETCH), install it "
        "(conda install -c bioconda mmseqs2), or set $ARDA_MMSEQS."
    )


def _auto_fetch() -> str | None:
    """Download a static mmseqs binary into ``bin/``; return its path or None."""
    try:
        from ._mmseqs_fetch import fetch

        return str(fetch(bin_dir()))
    except Exception as exc:  # noqa: BLE001 - network/layout failures are non-fatal here
        import warnings
        warnings.warn(f"mmseqs auto-fetch failed: {exc}", RuntimeWarning, stacklevel=2)
        return None


def run(args: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    """Run ``mmseqs <args>`` capturing stdout/stderr."""
    cmd = [mmseqs_binary(), *map(str, args)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise MMseqsError(
            f"`{' '.join(cmd)}` failed (exit {proc.returncode}):\n{proc.stderr}"
        )
    return proc


@lru_cache(maxsize=1)
def version() -> str:
    """Return the mmseqs version string (e.g. ``18-8cc5c``)."""
    return run(["version"]).stdout.strip()


def version_key(v: str) -> str:
    """Canonical form of a version string, for comparing an index marker to the running mmseqs.

    ``mmseqs version`` prints the same release+commit with different punctuation across builds:
    the official static binary says ``18-8cc5c``, the bioconda build ``18.8cc5c``. They are the
    same mmseqs and produce byte-compatible indexes; only the separator differs. Comparing the raw
    strings with ``==`` therefore rejected every committed index the moment the toolchain moved from
    conda's mmseqs to the static one -- so the precompiled DBs shipped in ``database/`` were never
    used and every run rebuilt a private cache instead.

    Fold each run of separator characters to a single ``-`` and lowercase. This bridges the cosmetic
    difference while still distinguishing genuinely different versions (``17-b804f`` != ``18-8cc5c``),
    so an incompatible index is never accepted.

    Not sufficient on its own -- see :func:`versions_compatible`, which is what callers should use.
    """
    return re.sub(r"[\s._-]+", "-", v.strip().lower())


def _commit_token(v: str) -> str:
    """The build's git commit from a version string, or ``""`` if it carries none."""
    toks = re.split(r"[\s._-]+", v.strip().lower())
    hexes = [t for t in toks if len(t) >= 5 and all(c in "0123456789abcdef" for c in t)]
    return hexes[-1] if hexes else ""


def versions_compatible(a: str, b: str) -> bool:
    """Do two ``mmseqs version`` strings denote builds with interchangeable index formats?

    Punctuation is not the only way the same build spells itself. The official **static release
    asset prints its full 40-char commit hash** (``8cc5ce367b5638c4306c2d7cfc652dd099a4643f``)
    while the bioconda build and the committed index marker print release+short-commit
    (``18.8cc5c`` / ``18-8cc5c``). Release 18 *is* commit ``8cc5c...``, so those are one build --
    but no amount of separator folding makes the strings equal.

    That mattered concretely: arda's own auto-fetched binary is the static asset, so a pure
    :func:`version_key` comparison rejected the index arda itself ships, on every macOS install.

    So: compare commit hashes when both carry one, accepting a prefix match in either direction
    (short vs full form). Fall back to :func:`version_key` when one has no hash at all. A genuinely
    different build has a different commit (``76da68ad...`` is not ``8cc5c...``) and is still
    rejected.
    """
    ca, cb = _commit_token(a), _commit_token(b)
    if ca and cb:
        return ca.startswith(cb) or cb.startswith(ca)
    return version_key(a) == version_key(b)


def createdb(fasta: str | Path, db: str | Path, *, dbtype: int | None = None) -> Path:
    """Create an mmseqs sequence DB from a FASTA file.

    ``dbtype``: ``None`` auto-detect, ``1`` amino-acid, ``2`` nucleotide.
    """
    args = ["createdb", str(fasta), str(db)]
    if dbtype is not None:
        args += ["--dbtype", str(dbtype)]
    run(args)
    return Path(db)


def search(
    query_db: str | Path,
    target_db: str | Path,
    result_db: str | Path,
    tmp_dir: str | Path,
    *,
    search_type: int = SEARCH_TYPE_AUTO,
    sensitivity: float = 5.7,
    evalue: float = 1e-3,
    max_seqs: int = 300,
    threads: int = 1,
    kmer: int | None = None,
    extra: list[str] | None = None,
) -> Path:
    """Run ``mmseqs search`` with backtrace enabled (``-a``).

    Args:
        kmer: MMseqs2 ``-k``. **This is the memory knob.** The nucleotide prefilter allocates a
            k-mer index table of 4**k entries, so the default k=15 costs 4**15 * 8 B ~ 8.6 GB
            regardless of database size, thread count or chunk size. k=13 costs ~0.7 GB. Leave
            ``None`` for MMseqs2's own default.
    """
    args = [
        "search", str(query_db), str(target_db), str(result_db), str(tmp_dir),
        "--search-type", str(search_type),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "--threads", str(threads),
        "-a",  # keep backtrace so convertalis can emit cigar/qaln/taln
    ]
    if kmer is not None:
        args += ["-k", str(kmer)]
    if extra:
        args += extra
    # Escape hatch for A/B-ing aligner options against a LIVE reference, without a rebuild per leg.
    # Not a supported interface -- MMseqs2 flags are not part of arda's contract, and several of
    # them (`--min-seq-id`, `-e`) silently change which reads are reported. It exists because every
    # recorded flag measurement in this project was taken against a reference that has since been
    # rebuilt twice, so re-measuring has to be cheaper than editing the source.
    if opts := os.environ.get("ARDA_MMSEQS_SEARCH_OPTS", "").split():
        args += opts
    run(args)
    return Path(result_db)


def top_hit(result_db: str | Path, out_db: str | Path) -> Path:
    """Reduce an alignment DB to the single best-scoring hit per query.

    MMseqs2 already stores each query's results sorted by descending score, so taking the first line
    per entry *is* the best hit -- this is the idiom mmseqs' own `filterdb` usage message shows.
    Verified on 100 k reads against the human reference: 4,101 queries, identical target and identical
    bit score to a full sort-and-dedupe in polars, on every one.

    Why it matters: with `--max-seqs 300`, 4,101 hitting queries produced **804,341** alignment rows
    (194 MB of TSV, each row carrying `cigar`/`qaln`/`taln`). Parsing that dominated arda's peak RSS --
    877 MB, against 284 MB for the mmseqs subprocess itself. Reducing before `convertalis` writes
    1.0 MB instead, and costs 0.04 s.
    """
    run(["filterdb", str(result_db), str(out_db), "--extract-lines", "1"])
    return Path(out_db)


def convertalis(
    query_db: str | Path,
    target_db: str | Path,
    result_db: str | Path,
    out_tsv: str | Path,
    *,
    format_output: str = DEFAULT_FORMAT_OUTPUT,
    threads: int = 1,
    search_type: int | None = None,
) -> Path:
    """Convert an alignment result DB to a TSV with the requested columns.

    ``search_type`` must be passed for nucleotide results (3) so convertalis can
    interpret the alignment; otherwise mmseqs cannot tell nt from translated.
    """
    args = [
        "convertalis", str(query_db), str(target_db), str(result_db), str(out_tsv),
        "--format-output", format_output,
        "--threads", str(threads),
    ]
    if search_type is not None:
        args += ["--search-type", str(search_type)]
    run(args)
    return Path(out_tsv)


def easy_search(
    query_fasta: str | Path,
    target_fasta_or_db: str | Path,
    out_tsv: str | Path,
    tmp_dir: str | Path,
    *,
    search_type: int = SEARCH_TYPE_AUTO,
    sensitivity: float = 5.7,
    evalue: float = 1e-3,
    max_seqs: int = 300,
    threads: int = 1,
    format_output: str = DEFAULT_FORMAT_OUTPUT,
    strand: int | None = None,
    extra: list[str] | None = None,
) -> Path:
    """One-shot createdb+search+convertalis producing a TSV.

    ``strand`` (nucleotide search only): 1 forward, 2 both strands; ``None`` lets
    mmseqs default (forward).
    """
    args = [
        "easy-search", str(query_fasta), str(target_fasta_or_db),
        str(out_tsv), str(tmp_dir),
        "--search-type", str(search_type),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "--threads", str(threads),
        "--format-output", format_output,
    ]
    if strand is not None:
        args += ["--strand", str(strand)]
    if extra:
        args += extra
    run(args)
    return Path(out_tsv)
