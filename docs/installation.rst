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

If you install with plain ``pip`` (no conda env) and ``mmseqs`` is not on
``PATH``, arda **auto-fetches** a static MMseqs2 binary into ``bin/mmseqs`` on
first use — no manual install needed. Controls:

* ``$ARDA_MMSEQS`` — use a specific mmseqs binary (highest priority).
* ``$ARDA_MMSEQS_ASSET`` — override the release asset (e.g.
  ``mmseqs-linux-sse41.tar.gz`` on pre-AVX2 CPUs).
* ``$ARDA_NO_AUTO_FETCH`` — disable auto-fetch (then install mmseqs yourself).

Fetch eagerly with ``python scripts/fetch_mmseqs.py`` (``setup.sh --no-conda``
does this for you).

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
