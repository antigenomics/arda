Samples split across files
==========================

Illumina writes one FASTQ per lane (``PT01_S1_L001_R1_001.fastq.gz``, ``…_L002_…``), and in-house
pipelines chunk a run their own way. Those files are **one repertoire** and must give **one**
clonotype table. arda calls each ``(R1, R2)`` pair a **read group** and takes as many as you like
per sample. Do not ``cat`` them first — arda concatenates after Stage 1, which is byte-identical
and skips a full copy of the data.

On the command line
-------------------

``--r1``, ``--r2`` and ``--id`` are repeatable and matched **by position**. Repeat an id to merge
those read groups into one sample; the id becomes the output basename, so ``--out-prefix`` is not
used alongside it.

.. code-block:: bash

   # two lanes, one sample -> results/PT01.clones.tsv
   arda rnaseq -d results/ \
       --r1 PT01_S1_L001_R1_001.fastq.gz --r2 PT01_S1_L001_R2_001.fastq.gz --id PT01 \
       --r1 PT01_S1_L002_R1_001.fastq.gz --r2 PT01_S1_L002_R2_001.fastq.gz --id PT01

Drop ``--r2`` for single-end. A sample may not mix the two, and more than one ``--r1`` without
``--id`` is **refused rather than guessed**: given ``RNA-SAMPLE_ID:12:00XX919:3_1.fastq.gz`` no
rule can say which field is the sample, and guessing wrong splits one repertoire into four with no
error anywhere.

.. _samples-sheet:

From a sample sheet
-------------------

``--samples sheet.tsv`` (or ``.csv``) replaces ``--r1``/``--r2``/``--id``. The columns are
nf-core's, so an existing nf-core samplesheet works unmodified:

.. code-block:: text

   sample	fastq_1	fastq_2
   PT01	PT01_S1_L001_R1_001.fastq.gz	PT01_S1_L001_R2_001.fastq.gz
   PT01	PT01_S1_L002_R1_001.fastq.gz	PT01_S1_L002_R2_001.fastq.gz
   PT02	PT02_S2_L001_R1_001.fastq.gz	PT02_S2_L001_R2_001.fastq.gz

.. code-block:: bash

   arda rnaseq --samples sheet.tsv -d results/

Two samples, three read groups. Repeated ``sample`` values merge **in row order** — nf-core's
re-sequencing rule. ``fastq_2`` may be blank for single-end. Extra columns (``strandedness``,
``seq_platform``, …) are ignored but named once in a warning, because a header spelling
``fastq2`` is not a sheet with no R2 and treating it as single-end would halve the data. Relative
paths resolve against the **sheet's own directory**, so a sheet travels with its data.

.. warning::

   **Never sort read groups by filename.** ``A_L010`` sorts before ``A_L002``, and the clonotype
   fold is not permutation-invariant: ``correct`` collapses an error child onto the parent it meets
   first. Declared order — command-line order, or sheet row order — is the order arda uses.

Why the result is byte-identical
--------------------------------

Stage 1 (``map``) is per-read and shards perfectly; Stages 2–3 (``assemble``, ``correct``) are
global and do not shard at all. So arda maps each read group, concatenates the Stage-1 AIRR **in
declared order**, and runs Stages 2–3 once over the whole sample. Each file is a contiguous run of
the library, and ``map`` output does not depend on where chunk boundaries fall, because a fragment
is never split across one.

Verified on ``tests/data/rnaseq_real`` cut into four read groups: ``<id>.airr.tsv``,
``<id>.assembled.airr.tsv`` and ``<id>.clones.tsv`` are byte-for-byte identical to the one-file
run.

.. note::

   Splitting one sample into two *samples* instead costs the contigs that tile across the split:
   on the test fixture, 57 reads as one sample against 41 + 17 = 58 as two, because reads that
   together span a long CDR3 never meet. With Stage 3 off the two partition exactly
   (35 + 16 = 51 reads, 34 + 15 = 49 clonotypes), which is what says ``correct`` neither
   double-counts nor drops anything at a sample boundary.

.. _samples-readgroups:

One worker per read group
-------------------------

Several samples from one CLI call run **one at a time**, each with every core — MMseqs2 threads
internally, so *N* samples at ``cores/N`` each is slower than *N* in a row. Parallelism belongs to
the scheduler, and there the unit is the read group: a sheet of 5 samples × 4 lanes is **20
independent jobs, not 5**.

.. code-block:: bash

   # SLURM: map array over read groups, then reduce array over samples
   arda cluster submit-samples --samples sheet.tsv --work-dir work/ -d results/ \
       --regime rnaseq --threads 16 --partition medium --submit

   # any other scheduler: emit the work units and the two commands, filled in
   arda cluster plan --samples sheet.tsv --work-dir work/ -d results/ --regime rnaseq --threads 16

``plan`` writes ``work/readgroups.tsv`` (``sample``, ``read_group``, ``r1``, ``r2``,
``part_airr``, ``part_report``) and ``work/samples.tsv``, creates the per-sample shard
directories, and prints one ``arda map`` per read-group row and one ``arda cluster reduce`` per
sample row.

.. important::

   The flags are printed filled in rather than as ``<flags>`` for a reason. ``arda rnaseq`` runs
   both stages in one call, so it can wire them together itself: ``--ec-mode rnaseq`` reads a
   column that Stage 1 writes **only** when asked with ``--junction-quality``. A scheduler running
   ``map`` and ``reduce`` as separate jobs has no such link — the column is simply absent, the
   quality-directed rescue silently never runs, and the exit code is 0. ``--regime`` is what picks
   both halves; do not hand-assemble them.

Both shipped workflow integrations schedule at the same granularity off the same sheet —
Snakemake (``integrations/snakemake/arda/``) and Nextflow
(``integrations/nextflow/arda/``, feed it a ``groupTuple`` so ``reads`` arrives pair-adjacent).
Both ask :func:`arda.cluster.regime_flags` for the two halves rather than restating the presets.
See :doc:`pipeline_integration`.

For one sample that arrived as one very large pair, shard it instead: :doc:`cluster`.
