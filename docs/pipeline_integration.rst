Pipeline integration
====================

arda is a plain CLI over named files, so it embeds in any workflow engine without glue code.
This page covers the one-shot RNA-seq command and the ready-made Nextflow module.

One-shot command
----------------

``arda rnaseq run`` runs ``map``, ``assemble`` and ``correct`` in a single call — the three steps
a bulk RNA-seq pipeline almost always wants together:

.. code-block:: bash

   arda rnaseq run --r1 R1.fq.gz --r2 R2.fq.gz --out-prefix SAMPLE --out-dir results/

It writes four files under ``--out-dir``:

===============================  ============================================================
``SAMPLE.clones.tsv``            corrected clonotype table (``junction, junction_aa, v_call,
                                 j_call, c_call, locus, duplicate_count, consensus_count``,
                                 plus ``d_call, d2_call, d_support, d2_support``)
``SAMPLE.airr.tsv``              mapped reads, AIRR Rearrangement schema
``SAMPLE.assembled.airr.tsv``    Stage-3 long-CDR3 reads rescued by contig assembly
                                 (omitted with ``--no-assemble``)
``SAMPLE.arda.json``             merged run report (map + assemble + correct: reads mapped,
                                 per-locus counts, isotype/constant fragments, timing, peak RSS)
===============================  ============================================================

Single-end input drops ``--r2``. The defaults match the individual commands (``--min-score 75``,
``--kmer 12`` for ~300 MB peak RSS, complete-junction clonotypes, D mapping on); use ``map``,
``assemble`` and ``correct`` separately when you need to tune their individual knobs, and
``--no-map-d`` to skip D in all three stages.

Nextflow module
---------------

A drop-in, nf-core-style local module ships in ``integrations/nextflow/arda/`` (``main.nf``,
``environment.yml``, ``Dockerfile``, ``nextflow.config``, ``README.md``). It wraps one
``arda rnaseq run`` call per sample, emits a ``versions.yml``, and publishes to
``${params.outdir}/arda/``.

Requirements
~~~~~~~~~~~~~

arda (PyPI ``arda-mapper``) plus the ``mmseqs2`` binary, both declared in ``environment.yml``.
``-profile conda`` builds the environment automatically; for ``-profile docker``/``singularity``,
build the image from the provided ``Dockerfile`` and push it to your registry.

Drop into an nf-core/rnaseq pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module consumes the same per-sample FASTQ channel the aligners do, so the sample sheet is
unchanged. Five edits, each mirroring how an existing tool is wired:

#. Copy ``integrations/nextflow/arda/`` to ``modules/local/arda/`` in your pipeline checkout.
#. Include and call it in ``workflows/rnaseq/main.nf`` on the trimmed, aligner-independent channel
   ``ch_strand_inferred_filtered_fastq``, and mix ``ARDA.out.versions`` into ``ch_versions``.
#. ``includeConfig "../../modules/local/arda/nextflow.config"`` from
   ``workflows/rnaseq/nextflow.config`` (sets ``ext.args`` and the publishDir).
#. Register a ``run_arda = false`` param in ``nextflow.config`` and ``nextflow_schema.json`` (strict
   schema validation).
#. For container profiles, pin ``withName: 'ARDA' { container = '<registry>/arda-mapper:2.5.6' }``
   in your deployment config.

Then run with ``--run_arda``. The module's ``README.md`` has copy-paste snippets and a standalone
one-process test harness.

Runtime & resources
~~~~~~~~~~~~~~~~~~~~~

arda is **CPU-bound** — the MMseqs2 search dominates — so give it cores. Measured on bulk tumor
RNA-seq at 32 cores, ``arda rnaseq run --assemble`` (map + assemble + correct):

.. list-table::
   :header-rows: 1

   * - reads
     - cores
     - wall time
     - throughput
     - peak RSS
     - clonotypes
   * - 104.9 M (2×150)
     - 32
     - 46 min
     - ~38,300 reads/s
     - 2,707 MB
     - 28,444
   * - 153.9 M (2×100)
     - 32
     - 51 min
     - ~50,300 reads/s
     - 511 MB
     - 1,519
   * - 138.5 M (2×100)
     - 32
     - 45 min
     - ~51,900 reads/s
     - 314 MB
     - 30

Throughput scales roughly linearly with cores; a typical full-depth bulk RNA-seq sample takes ~45 min
on 32 cores. The Nextflow module is labelled ``process_high``; raise it with
``withName: 'ARDA' { cpus = 32 }`` for full-depth data. ``--threads`` follows ``task.cpus``.

**Peak memory tracks repertoire richness, not read depth.** The third sample above has *more* reads
than the first yet uses 9× less RAM — it carries 30 clonotypes against 28,444. Mapping alone is flat
and small (~300–400 MB at any depth); it is Stage 3 (assembly + error correction) that holds the clone
set in memory. Budget **~4 GB** for a lymphocyte- or B-cell-rich sample and ~1 GB for a cold one. Under
a hard memory cap, ``--no-assemble`` restores the flat mapping-only profile — at the cost of the long
CDR3s that no single read spans.

Tuning
~~~~~~~

All arda flags pass through the module's ``ext.args`` (``--organism mouse``, ``--reconstruct``,
``--min-score 0``, ``--kmer 11`` …). ``--threads`` is wired to ``task.cpus`` automatically.
