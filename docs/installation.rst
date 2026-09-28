Installation
============

.. code-block:: bash

   pip install arda-mapper

That is the whole installation for normal use. The distribution is ``arda-mapper``; it imports as
``arda``. Binary wheels are published for CPython 3.10–3.13 on Linux, macOS (arm64) and Windows,
and carry arda's four C++ extensions prebuilt.

Requires Python 3.10 or later.

.. contents::
   :local:
   :depth: 1

Verify
------

.. code-block:: bash

   arda --version
   arda info

``arda info`` prints the resolved reference path, the cache location and which external tools arda
can see. Run it first whenever anything looks wrong.

What arda fetches at runtime
----------------------------

Two things arda needs are not in the wheel, and both are fetched automatically on first use into
``$XDG_CACHE_HOME/arda`` (by default ``~/.cache/arda``):

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - What
     - When
   * - The germline reference
     - On the first command that annotates anything. About 3 MB, taken from the matching GitHub
       release asset. The MMseqs2 index is then built in the cache.
   * - The MMseqs2 binary
     - On the first search, unless a suitable ``mmseqs`` is already resolvable. A static build, so
       there is no conda requirement.

IgBLAST is fetched only if you run ``arda igblast`` or rebuild the reference. Annotation never
needs it.

.. note::

   A source checkout uses the committed ``database/vdj/<organism>/`` references instead, including
   precompiled MMseqs2 indexes, so it needs no download and no build step.

Supported organisms
-------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Organism
     - Loci
     - Note
   * - ``human``
     - IG + TR
     -
   * - ``mouse``
     - IG + TR
     -
   * - ``rat``
     - IG only
     - IMGT carries no TR reference for these three, so their TR loci build empty and are
       recorded as such in ``loci_manifest.tsv``.
   * - ``rabbit``
     - IG only
     -
   * - ``rhesus_monkey``
     - IG only
     -

Development checkout
--------------------

``setup.sh`` bootstraps a development environment with `uv <https://docs.astral.sh/uv/>`_. Conda is
used only by the Nextflow integration, which ships its own ``environment.yml``.

.. code-block:: bash

   git clone https://github.com/antigenomics/arda
   cd arda
   bash setup.sh
   source .venv/bin/activate

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Flag
     - Effect
   * - ``--build-db``
     - Rebuild the reference database after install. Needs IgBLAST.
   * - ``--tests``
     - Run the unit and synthetic suites, and fail if they fail.

What ``setup.sh`` does:

#. Removes any stale ``build/`` directory.
#. Creates ``.venv`` with ``uv`` and installs ``-e .[test,dev]`` with ``--no-build-isolation``, so
   the editable on-import rebuild can find ``nanobind``.
#. Downloads the current IgBLAST release into ``bin/`` (gitignored).
#. Fetches a static MMseqs2 binary into ``bin/`` unless one is already on ``PATH``.

It then checks four things that a bare ``import arda`` does not:

#. ``arda._markup``, ``arda._segmap`` and ``arda._denoise`` actually imported. arda falls back to a
   pure-Python markup path when they are missing, so a failed C++ build otherwise looks like a
   successful install and surfaces later as an unexplained slowdown.
#. ``arda.__version__`` agrees with ``pyproject.toml``'s ``version``.
#. ``mmseqs`` resolves, and reports its version.
#. Every mode and stage command resolves on the CLI. A deploy into the wrong environment prints a
   correct version and still lacks the commands.

The editable install rebuilds the C++ extensions on import, so editing a ``.cpp`` file takes effect
on the next ``import arda`` with no manual build step.

.. note::

   ``build/`` is not a cache you can safely delete on its own. scikit-build-core stores CMake's
   configuration there, including the absolute path of the interpreter it configured against. If
   you remove it, reinstall (``pip install -e . --no-build-isolation``) rather than relying on the
   on-import rebuild, which builds but does not configure.

Running the test suites:

.. code-block:: bash

   python -m pytest tests/unit tests/synthetic tests/realworld -q   # the CI gate
   ruff check src/
   make -C docs html                                                # zero warnings required

``tests/realworld`` compares against IgBLAST on committed fixtures and runs offline.
``tests/benchmark`` is skipped unless ``RUN_BENCHMARK=1``. Optional extras gate optional suites:
``.[groundtruth]`` (``olga``) for the generative ground-truth tests and ``.[test]`` for ``airr``
schema validation — without them those tests skip, so install ``-e '.[test]'`` before reading a
green suite as full coverage.

Choosing the MMseqs2 binary
---------------------------

You do not have to install MMseqs2. Resolution order:

#. ``$ARDA_MMSEQS``
#. an ``arda_mmseqs`` package, if one is installed
#. ``<project>/bin/mmseqs``
#. ``mmseqs`` on ``PATH``
#. auto-fetch a static build

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Effect
   * - ``ARDA_MMSEQS``
     - Use this binary. Highest priority, and deliberately not version-checked — an explicit
       override is your call.
   * - ``ARDA_MMSEQS_ASSET``
     - Override the release asset, for example ``mmseqs-linux-sse41.tar.gz`` on a pre-AVX2 CPU.
   * - ``ARDA_NO_AUTO_FETCH``
     - Refuse to download. arda then errors and names the fix.

Fetch eagerly with ``python scripts/fetch_mmseqs.py``; ``setup.sh`` does this for you.

.. important::

   **Candidates are version-matched, not merely found.** An MMseqs2 index is only reusable by the
   release that compiled it. An unrelated ``mmseqs`` on ``PATH`` would make arda discard the
   precompiled indexes in ``database/`` and rebuild a private cache — costing a slow first run and,
   if the two releases align differently, results not comparable with anyone else's. arda checks
   each candidate against the shipped index marker, falls back to a known-good build rather than
   accepting a mismatch, and warns with the consequence named if it can find neither.

   Two exceptions: ``$ARDA_MMSEQS`` is never version-checked, and where no precompiled index ships
   — which is the case for a plain ``pip install``, since the packaged reference omits the indexes
   — there is nothing to match and any working mmseqs is accepted.

The ``[mmseqs]`` and ``[rnaseq]`` extras still resolve but install nothing. They are kept as empty
aliases so existing pins continue to resolve; everything the bulk RNA-seq pipeline needs, including
``seqtree`` and ``dnaio``, is a core dependency.

IgBLAST without a checkout
--------------------------

``arda igblast`` runs IgBLAST and emits AIRR. It is how the gold-standard comparisons are produced,
and it is the only part of arda that needs IgBLAST at all.

.. code-block:: bash

   arda igblast -i reads.fq -o truth.tsv      # fetches IgBLAST once, then runs

Resolution order: ``$ARDA_IGBLAST`` → ``<project>/bin`` if a checkout already has one →
``<cache>/igblast``, auto-fetched.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Effect
   * - ``ARDA_IGBLAST``
     - A directory holding ``igblastn`` and ``internal_data/``. Highest priority, never fetched.
       Point it at a conda or system IgBLAST to reuse one.
   * - ``ARDA_IGBLAST_ASSET``
     - Override the NCBI asset suffix: ``x64-linux``, ``x64-macosx``, ``x64-win64``.
   * - ``ARDA_NO_AUTO_FETCH``
     - Refuse to download; arda errors and names the fix rather than proceeding without IgBLAST.

``arda.igblast.igblast_version()`` reports which NCBI release is installed, so a benchmark can
record it. Fetch eagerly with ``python scripts/fetch_igblast.py --dest <dir>``.

.. note::

   IgBLAST needs its J-frame table, ``optional_file/<org>_gl.aux``, to report a junction. Without
   it, IgBLAST returns V and J normally and an empty ``junction_aa`` on every read, at exit 0 —
   indistinguishable from a truth that genuinely has none. arda therefore raises when the file is
   missing rather than degrading. All five ``*_gl.aux`` tables ship with the auto-fetched release.

Air-gapped and offline installs
-------------------------------

.. code-block:: bash

   export ARDA_NO_AUTO_FETCH=1

With auto-fetch disabled, arda uses a pre-populated cache and errors with the missing path named
rather than reaching for the network. Populate the cache by copying
``$XDG_CACHE_HOME/arda/database`` from a machine that has run once, or by unpacking the
``arda-reference-vdj.tar.gz`` release asset there.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Variable
     - Effect
   * - ``ARDA_HOME``
     - Point at a **source checkout**. This also disables auto-fetch. To simulate a pip user from a
       checkout, set ``XDG_CACHE_HOME`` instead — the cache is
       ``$XDG_CACHE_HOME/``\ **arda**\ ``/database``.
   * - ``ARDA_NO_AUTO_FETCH``
     - Disable every download.

``arda qc report`` is designed for this setting: the HTML it writes inlines its data and has no
external reference of any kind, so it opens on a login node with no network.

.. _segments.fasta:

The segment reference is generated, not shipped
-----------------------------------------------

``--two-pass`` — and therefore ``--fast-segments`` and ``--v-only-on-segment`` — nominates
candidates from a second, much smaller reference: ``segments.fasta``, the 924 collapsed per-allele
V/J/C targets. That file is generated rather than shipped, and it is not in the auto-fetched
reference tarball.

It is built automatically on first use, under the same build lock the MMseqs2 index takes, in about
0.3 s per organism per cache. Running ``arda build-index`` beforehand is not required.

Both generation and regeneration take ``.segments.build.lock`` in the reference directory and build
into a staging path, so concurrent runs — a Nextflow process per sample, a SLURM task per shard —
cannot race each other into a partial file. If generation fails, arda warns and falls back to the
one-pass search rather than raising.

.. note::

   Before this was automatic, a plain ``pip install`` had no ``segments.fasta``, so every
   ``--two-pass`` run silently degraded to the one-pass search behind a single log line: correct
   output, exit 0, and none of the speedup. If you are upgrading, nothing needs doing — a
   ``segments.fasta`` written before arda 2.8.0 is detected by format and regenerated rather than
   used.

Export the generated reference, or any part of it, with :doc:`arda export-ref <reference_export>`.

Troubleshooting
---------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Symptom
     - Cause and fix
   * - Annotation is far slower than expected
     - The C++ extensions did not build and arda is on its pure-Python markup path. Check with
       ``python -c "import arda._markup"``; reinstall from a wheel, or run ``setup.sh``.
   * - ``--two-pass`` produces correct output with no speedup
     - Check the log for the one-pass fallback line, then confirm ``segments.fasta`` exists in the
       reference directory.
   * - ``0/N reads mapped``, exit 0
     - Usually a reference or index problem rather than the data. Run ``arda info``, and check
       whether a concurrent run was building the index.
   * - A version of arda without the commands you expect
     - A deploy into the wrong environment. ``arda --version`` succeeds without the C++ extensions
       and without the reference, so check ``arda info`` too.
   * - MMseqs2 index rebuilt on every run
     - The resolved ``mmseqs`` version does not match the shipped index. Either pin the matching
       binary with ``ARDA_MMSEQS``, or run ``arda build-index`` once for your version.
