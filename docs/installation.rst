Installation
============

``arda`` installs with pip and builds three small C++ extensions (``_markup``, ``_segmap``,
``_denoise``). The MMseqs2 binary is fetched and version-matched at runtime — there is nothing to
install by hand.

.. code-block:: bash

   pip install arda-mapper        # imports as `arda`; binary wheels ship the extensions

Development bootstrap
---------------------

``setup.sh`` is **uv-based, not conda**. Conda is used only by the Nextflow integration, which
ships its own ``environment.yml``.

.. code-block:: bash

   bash setup.sh
   source .venv/bin/activate

``setup.sh`` flags:

* ``--build-db`` — rebuild the reference database after install (needs IgBLAST).
* ``--tests`` — run the unit + synthetic suites, and **fail** if they fail.

What it does, and what it verifies
-----------------------------------

* Wipes any stale ``build/``. ⛔ Not cosmetic: scikit-build-core caches CMake's configuration
  including the **absolute path** of the interpreter it configured against, so a ``build/`` left
  by a venv that no longer exists makes every later on-import rebuild fail with
  *"Could NOT find Python"* — and arda then falls back to its pure-Python markup path.
* Creates ``.venv`` with ``uv`` and installs ``-e .[test,dev]`` with ``--no-build-isolation``, so
  the editable on-import rebuild can find ``pybind11`` (pinned ``>=3.0.2,<4``, matching
  ``pyproject.toml``'s build-system: ``PYBIND11_MODULE`` changed to multi-phase init in 3.0.0).
* Downloads the latest IgBLAST release into ``bin/`` (gitignored) — needed only to rebuild
  references and to run ``arda igblast``, never for annotation.
* Fetches a static MMseqs2 binary into ``bin/`` unless one is already on ``PATH``.

Then it checks four things, because a bare ``import arda`` checks none of them:

#. ``arda._markup``, ``arda._segmap`` and ``arda._denoise`` actually imported. arda **falls back
   to a pure-Python markup path** when they are missing, so a failed C++ build looks like a
   successful install and surfaces much later as a silent slowdown.
#. ``arda.__version__`` — which a unit test now pins to ``pyproject.toml``'s ``version``, because
   they are two independent literals and a release once had them disagree.
#. ``mmseqs`` resolves and reports its version.
#. Every mode and stage command resolves on the **CLI** — ``rnaseq``, ``amplicon``,
   ``singlecell``, ``map``, ``correct``, ``assemble``, ``shm``, ``cluster``, ``annotate``. A
   deploy into the wrong environment prints a correct version and still lacks the commands.

MMseqs2 without conda
---------------------

You do not have to install MMseqs2. This is enough::

    pip install arda-mapper              # auto-fetches a static binary on first use

Resolution order: ``$ARDA_MMSEQS`` → an ``arda_mmseqs`` package if one is installed →
``<project>/bin/mmseqs`` → ``mmseqs`` on ``PATH`` → auto-fetch.

The ``[mmseqs]`` extra still resolves, but it installs nothing: the companion ``arda-mmseqs``
distribution that used to bundle the binary is not published, and an extra naming a package PyPI
does not have is a hard failure rather than a soft one. Auto-fetch pulls the same static
binary.

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

.. _segments.fasta:

The segment reference is generated, not shipped
-----------------------------------------------

``--two-pass`` (and therefore ``--fast-segments`` / ``--v-only-on-segment``) nominates
candidates from a second, much smaller reference: ``segments.fasta``, the 924 collapsed per-allele
V / J / C targets. That file is **generated, not shipped** — it is not in the auto-fetched
reference tarball.

**It is now built automatically on first use**, under the same build lock the MMseqs2 index
uses, taking about 0.3 s once per organism per cache. Running ``arda build-index`` beforehand is
no longer required for the two-pass configurations.

.. note::

   Before this, a plain ``pip install`` had no ``segments.fasta``, so every ``--two-pass`` run
   silently degraded to the one-pass search behind a single log line: the flagship amplicon
   configuration was unreachable out of the box, with correct output and none of the speedup.
   If you are upgrading, nothing needs doing — a ``segments.fasta`` written before arda 2.8.0
   is detected **by format** (it carries the old ``JC|`` targets) and regenerated, rather than
   used.

Both the generation and the regeneration take ``.segments.build.lock`` in the reference
directory and build into a staging path, so concurrent runs — a Nextflow process per sample, a
SLURM task per shard — cannot race each other into a partial file. If generation fails, arda
warns and falls back to the one-pass search rather than raising.

Export the generated reference, or any part of it, with :doc:`arda export-ref
<reference_export>`.

.. _IgBLAST without a checkout:

IgBLAST without a checkout
--------------------------

``arda igblast`` runs IgBLAST and emits AIRR — it is how the gold-standard comparisons in the
benchmark are produced, and it is the only part of arda that needs IgBLAST at all. Annotation
does not.

It needs no setup. arda downloads the current NCBI release on first use, the same way it
handles MMseqs2::

    pip install arda-mapper
    arda igblast -i reads.fq -o truth.tsv    # fetches IgBLAST once, then runs

Resolution order: ``$ARDA_IGBLAST`` → ``<project>/bin`` if a checkout already has one (what
``setup.sh`` produces) → ``<cache>/igblast``, auto-fetched.

Controls:

* ``$ARDA_IGBLAST`` — a directory holding ``igblastn`` and ``internal_data/`` (highest
  priority, never fetched). Point this at a conda or system IgBLAST to reuse it.
* ``$ARDA_IGBLAST_ASSET`` — override the NCBI asset suffix (``x64-linux``, ``x64-macosx``,
  ``x64-win64``).
* ``$ARDA_NO_AUTO_FETCH`` — refuse to download; arda then errors and names the fix rather
  than proceeding without IgBLAST.

``arda.igblast.igblast_version()`` reports which NCBI release is installed, so a benchmark can
record it. Fetch eagerly with ``python scripts/fetch_igblast.py --dest <dir>``.

**PyPI install (no source tree).** ``pip install arda-mapper`` ships code only. On first use
it **auto-fetches** the curated ``vdj/`` references (the ``arda-reference-vdj.tar.gz`` release
asset, ~3 MB) into ``$XDG_CACHE_HOME/arda`` (default ``~/.cache/arda``) and builds the MMseqs2
index there — **no ``$ARDA_HOME`` and no reference build required**. Set
``ARDA_NO_AUTO_FETCH`` to disable the download (air-gapped runs with a pre-populated cache).

Everything ``arda rnaseq`` needs -- including ``seqtree``, the clonotype neighbour search used by
``arda correct`` -- comes with a plain ``pip install arda-mapper``. No extra. (``seqtree``
was an optional extra before 2.5.5, which meant a plain install could map and assemble a whole
sample and only then fail, before writing any clonotype table.)
