Installation
============

``arda`` uses a dedicated conda environment for the MMseqs2 binary and the C++
toolchain; the package itself installs with pip and builds a small C++ extension.

Bootstrap
---------

.. code-block:: bash

   bash setup.sh
   conda activate arda

``setup.sh`` flags:

* ``--no-conda`` — use the already-active environment instead of creating ``arda``.
* ``--build-db`` — rebuild the reference database after install (needs IgBLAST).
* ``--tests`` — run the fast unit + synthetic suites.

What gets installed
-------------------

* The ``arda`` conda env (Python, ``mmseqs2``, a C++ compiler, perl).
* The latest IgBLAST release into ``bin/`` (gitignored) — only needed to rebuild
  references, not at annotation time.
* The ``arda`` package + the ``arda._markup`` C++ extension (editable install).

MMseqs2 without conda
---------------------

You do not have to install MMseqs2. Either of these works with nothing else set up::

    pip install arda-mapper              # auto-fetches a static binary on first use
    pip install 'arda-mapper[mmseqs]'    # ships the binary in the wheel; no download at all

Resolution order: ``$ARDA_MMSEQS`` → the ``arda-mapper[mmseqs]`` wheel → ``<project>/bin/mmseqs``
→ ``mmseqs`` on ``PATH`` → auto-fetch.

**Candidates are version-matched, not merely found.** An mmseqs index is only reusable by the
release that compiled it, so an unrelated ``mmseqs`` on ``PATH`` makes arda discard the
precompiled indexes shipped in ``database/`` and rebuild a private cache — silently, costing a
slow first run and, if the two releases align differently, results that are not comparable with
anyone else's. arda therefore checks each candidate against the shipped index marker and falls
back to a known-good build rather than accepting a mismatch. If it cannot find or fetch a
matching one it warns and says what the consequence will be.

Two deliberate exceptions: ``$ARDA_MMSEQS`` is never version-checked (an explicit override is
your call), and when no precompiled index ships — which is the case for a plain ``pip install``,
since the packaged reference omits the indexes — there is nothing to match, so any working
mmseqs is accepted.

Controls:

* ``$ARDA_MMSEQS`` — use a specific mmseqs binary (highest priority, unchecked).
* ``$ARDA_MMSEQS_ASSET`` — override the release asset (e.g.
  ``mmseqs-linux-sse41.tar.gz`` on pre-AVX2 CPUs).
* ``$ARDA_NO_AUTO_FETCH`` — disable auto-fetch (then install mmseqs yourself).

Fetch eagerly with ``python scripts/fetch_mmseqs.py`` (``setup.sh`` does this for you).

The committed ``database/vdj/<organism>/`` references — including **precompiled
MMseqs2 indexes** under ``mmseqs/`` — mean a source checkout needs no build. The shipped
indexes are used when the local MMseqs2 version matches; otherwise arda rebuilds a private
cache on first run. ``arda build-index`` (re)builds the shipped indexes for your MMseqs2
version.

**PyPI install (no source tree).** ``pip install arda-mapper`` ships code only. On first use
it **auto-fetches** the curated ``vdj/`` references (the ``arda-reference-vdj.tar.gz`` release
asset, ~3 MB) into ``$XDG_CACHE_HOME/arda`` (default ``~/.cache/arda``) and builds the MMseqs2
index there — **no ``$ARDA_HOME`` and no reference build required**. Set
``ARDA_NO_AUTO_FETCH`` to disable the download (air-gapped runs with a pre-populated cache).

Everything ``arda rnaseq`` needs -- including ``seqtree``, the clonotype neighbour search used by
``arda rnaseq correct`` -- comes with a plain ``pip install arda-mapper``. No extra. (``seqtree``
was an optional extra before 2.5.5, which meant a plain install could map and assemble a whole
sample and only then fail, before writing any clonotype table.)
