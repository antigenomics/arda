Samples split across files
==========================

A sample arrives split across files more often than not. Illumina writes one FASTQ per lane —

.. code-block:: text

   PT01_S1_L001_R1_001.fastq.gz    PT01_S1_L001_R2_001.fastq.gz
   PT01_S1_L002_R1_001.fastq.gz    PT01_S1_L002_R2_001.fastq.gz

where ``PT01`` is the sample name from the sample sheet, ``S1`` its position in that sheet,
``L001`` the flow-cell lane, ``R1``/``R2`` the read designation and ``001`` a segment number that
is always ``001``. (``bcl2fastq --no-lane-splitting`` drops the ``L00#`` field and merges the lanes
for you; most sites do not use it.) In-house pipelines chunk a run their own way:

.. code-block:: text

   RNA-SAMPLE_ID:12:00XX919:1_1.fastq.gz    RNA-SAMPLE_ID:12:00XX919:1_2.fastq.gz
   RNA-SAMPLE_ID:12:00XX919:2_1.fastq.gz    RNA-SAMPLE_ID:12:00XX919:2_2.fastq.gz

Either way these are **one repertoire** and must give **one** clonotype table. arda calls each
``(R1, R2)`` pair a **read group**, takes as many as you like per sample, and writes one set of
outputs per sample.

Do not ``cat`` them first. arda maps each read group and concatenates *after* Stage 1, which is
byte-identical to the same reads in one file and skips a full copy of the data.

The grouping is declared, never guessed
---------------------------------------

arda does **not** infer which files belong together from their names. Given
``RNA-SAMPLE_ID:12:00XX919:3_1.fastq.gz``, no rule can say which colon-separated field is the
sample without being told — and a rule that guesses wrong splits one repertoire into four with no
error anywhere. So more than one ``--r1`` without ``--id`` is refused rather than guessed at.

.. _samples-options:

On the command line
-------------------

``--r1``, ``--r2`` and ``--id`` are repeatable and matched **by position**: the *i*-th ``--r2`` is
the mate of the *i*-th ``--r1``, and the *i*-th ``--id`` names the sample it belongs to. Repeat an
id to merge those read groups into one sample.

One pair is the ordinary case and is unchanged:

.. code-block:: bash

   arda rnaseq --r1 S_1.fq.gz --r2 S_2.fq.gz --out-prefix S -d results/

Five pairs, five samples:

.. code-block:: bash

   arda rnaseq -d results/ \
       --r1 a_1.fq.gz --r2 a_2.fq.gz --id A \
       --r1 b_1.fq.gz --r2 b_2.fq.gz --id B \
       --r1 c_1.fq.gz --r2 c_2.fq.gz --id C \
       --r1 d_1.fq.gz --r2 d_2.fq.gz --id D \
       --r1 e_1.fq.gz --r2 e_2.fq.gz --id E

Five pairs, three samples — the first three read groups are one repertoire:

.. code-block:: bash

   arda rnaseq -d results/ \
       --r1 'RNA-SAMPLE_ID:12:00XX919:1_1.fastq.gz' --r2 'RNA-SAMPLE_ID:12:00XX919:1_2.fastq.gz' --id A \
       --r1 'RNA-SAMPLE_ID:12:00XX919:2_1.fastq.gz' --r2 'RNA-SAMPLE_ID:12:00XX919:2_2.fastq.gz' --id A \
       --r1 'RNA-SAMPLE_ID:12:00XX919:3_1.fastq.gz' --r2 'RNA-SAMPLE_ID:12:00XX919:3_2.fastq.gz' --id A \
       --r1 d_1.fq.gz --r2 d_2.fq.gz --id B \
       --r1 e_1.fq.gz --r2 e_2.fq.gz --id C

That writes ``results/A.clones.tsv``, ``results/B.clones.tsv`` and ``results/C.clones.tsv``. The
sample id is the output basename, so ``--out-prefix`` is not used with ``--id``.

Drop ``--r2`` entirely for single-end input. A sample may not mix the two: half its reads silently
losing their mate is not a configuration, and it is refused at parse time rather than surfacing as
a length mismatch several minutes into the run.

.. _samples-sheet:

From a sample sheet
-------------------

``--samples sheet.tsv`` (or ``.csv``) replaces ``--r1``/``--r2``/``--id``. The columns are
nf-core's, so an existing nf-core samplesheet works here unmodified:

.. code-block:: text

   sample	fastq_1	fastq_2
   PT01	PT01_S1_L001_R1_001.fastq.gz	PT01_S1_L001_R2_001.fastq.gz
   PT01	PT01_S1_L002_R1_001.fastq.gz	PT01_S1_L002_R2_001.fastq.gz
   PT02	PT02_S2_L001_R1_001.fastq.gz	PT02_S2_L001_R2_001.fastq.gz

.. code-block:: bash

   arda rnaseq --samples sheet.tsv -d results/

Two samples, three read groups. Repeated ``sample`` values merge **in row order** — nf-core's
re-sequencing rule, and row order is what makes the result reproducible. ``fastq_2`` may be blank
for single-end. Extra columns (``strandedness``, ``seq_platform``, …) are ignored, but named once
in a warning, because a header spelling ``fastq2`` is not a sheet with no R2 and treating it as
single-end would halve the data. Relative paths resolve against the **sheet's own directory**, so
a sheet travels with its data.

Why the result is byte-identical
--------------------------------

Stage 1 (``map``) is per-read and shards perfectly. Stages 2 and 3 (``assemble``, ``correct``) are
global and do not shard at all — so arda maps each read group, concatenates the Stage-1 AIRR **in
declared order**, and runs Stages 2–3 once over the whole sample. Two properties make that exact:

* each file is a contiguous run of the library, so concatenating in order reproduces the original
  row order; and
* ``map`` output does not depend on where chunk boundaries fall, because a fragment is never split
  across one.

Verified on ``tests/data/rnaseq_real`` cut into four read groups with ``arda cluster split``:
``<id>.airr.tsv``, ``<id>.assembled.airr.tsv`` and ``<id>.clones.tsv`` are byte-for-byte identical
to the one-file run.

.. warning::

   **Never sort read groups by filename.** ``A_L010`` sorts before ``A_L002``, and the clonotype
   fold is not permutation-invariant: ``correct`` collapses an error child onto the parent it meets
   first, and coverage assignment is first-with-longest-overlap-wins. Declared order — command-line
   order, or sheet row order — is the order arda uses.

.. note::

   Merging read groups is not optional bookkeeping. Splitting one sample's reads into two
   *samples* costs the contigs that tile across the split: on the test fixture, 57 reads as one
   sample against 41 + 17 = 58 as two, because reads that together span a long CDR3 never meet.
   With Stage 3 off the two do partition exactly (35 + 16 = 51 reads, 34 + 15 = 49 clonotypes),
   which is what says ``correct`` neither double-counts nor drops anything at a sample boundary.

.. _samples-readgroups:

One worker per read group
-------------------------

Several samples from one CLI call run **one at a time**, each with every core. That is
deliberate: MMseqs2 threads internally and is given all of them, so *N* samples at ``cores/N``
each is slower than *N* in a row, plus dispatch. Parallelism belongs to the scheduler — and there
the unit is the read group, not the sample. A sheet of 5 samples × 4 lanes is **20 independent
jobs, not 5**.

``arda cluster plan`` emits those work units for any scheduler:

.. code-block:: bash

   arda cluster plan --samples sheet.tsv --work-dir work/ -d results/ --regime rnaseq --threads 16

It writes ``work/readgroups.tsv`` (one row per read group: ``sample``, ``read_group``, ``r1``,
``r2``, ``part_airr``, ``part_report``), ``work/samples.tsv`` (one row per sample), creates the
per-sample shard directories, and prints the two commands **filled in**:

.. code-block:: bash

   # one worker per row of readgroups.tsv:
   arda map --r1 "$r1" ${r2:+--r2 "$r2"} -o "$part_airr" --report "$part_report" \
       --organism human --threads 16 --kmer 12 --min-score 75.0 --prefilter --junction-quality
   # then one per row of samples.tsv, after ALL of its read groups finish:
   arda cluster reduce --shard-dir "$shard_dir" --out-dir results/ --out-prefix "$out_prefix" \
       --organism human --threads 16 --ec-mode rnaseq --assemble --complete-only

.. important::

   The flags are printed filled in rather than as ``<flags>`` for a reason. ``arda rnaseq`` runs
   both stages in one call, so it can wire them together itself: ``--ec-mode rnaseq`` reads a
   column that Stage 1 writes **only** when asked with ``--junction-quality``. A scheduler running
   ``map`` and ``reduce`` as separate jobs has no such link — the column is simply absent, the
   quality-directed rescue silently never runs, and the exit code is 0. ``--regime`` is what picks
   both halves; do not hand-assemble them.

For SLURM, ``arda cluster submit-samples`` renders exactly that DAG as two arrays — ``map`` over
read groups, then ``reduce`` over samples gated on the whole map array — with no split step,
because the files already are the shards:

.. code-block:: bash

   arda cluster submit-samples --samples sheet.tsv --work-dir work/ -d results/ \
       --regime rnaseq --threads 16 --partition medium --submit

See :doc:`cluster` for the single-sample sharded path (``arda cluster submit --shards N``), which
is what you want when one sample arrives as one very large pair.

Workflow engines
----------------

Both shipped integrations schedule at read-group granularity and read the same sheet:

* **Snakemake** — ``integrations/snakemake/arda/``. ``snakemake -s Snakefile --config
  samples=sheet.tsv outdir=results -c 32``, or ``--profile slurm``. One ``map_read_group`` job per
  read group, one ``reduce_sample`` per sample.
* **Nextflow** — ``integrations/nextflow/arda/``. Feed it a ``groupTuple`` so ``reads`` arrives
  pair-adjacent; the module's ``README.md`` has the channel snippet. A single-lane sample is
  unaffected.

Both ask :func:`arda.cluster.regime_flags` for the two halves rather than restating the presets,
so they cannot drift from what the mode commands do.
