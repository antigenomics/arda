Command reference
=================

Every arda command, grouped by what you would reach for it to do. Each command's own
``--help`` is authoritative for its flags and their defaults; this page is the map.

.. tip::

   Run ``--help`` through a wide terminal. Long flag names are truncated at 80 columns, so
   ``--v-only-on-segment`` renders as ``--v-only-on-...`` and greps for it fail::

      COLUMNS=200 arda map --help

Global options
--------------

These come **before** the subcommand.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Effect
   * - ``--version``
     - Print the arda version and exit.
   * - ``-v``, ``--verbose``
     - Raise console verbosity. The default already prints stage lines and a throttled progress
       line; ``-v`` adds DEBUG with the level and module on each line. Repeatable.
   * - ``-q``, ``--quiet``
     - Console warnings and errors only. Does not quieten ``--log-file``.
   * - ``--log-file PATH``
     - Also write a DEBUG log here, timestamped, with the process peak RSS on every line. Always
       DEBUG whatever the console level is, so a quiet cluster job still leaves a full record.

.. code-block:: bash

   arda -v --log-file run.log amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/

Reads to clonotypes
-------------------

The three pipeline modes. Each runs mapping, contig assembly and error correction in one call and
carries the speed configuration its library type needs. See :doc:`usage` for how to choose.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Command
     - What it does
   * - .. _cli-rnaseq:

       ``arda rnaseq``
     - Bulk RNA-seq to clonotypes, for whole-transcriptome libraries where 0.02–3 % of reads are
       receptor-derived.
   * - .. _cli-amplicon:

       ``arda amplicon``
     - Targeted RepSeq / 5'RACE amplicon to clonotypes, for libraries whose reads span V into J.
   * - .. _cli-cells:

       ``arda cells``
     - Single cell. Takes one UMI consensus per molecule with the barcode in the record name,
       assembles each cell's contigs with no germline reference, annotates, pairs the chains and
       flags multiplets. See :doc:`singlecell`.

``arda singlecell`` exists but is **reserved**: it is a placeholder kept so the name cannot be
taken for a third speed preset, and the single-cell work is in ``arda cells``.

``rnaseq`` and ``amplicon`` share their interface:

.. code-block:: bash

   arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 8
   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 8
   arda rnaseq   --samples sheet.tsv -d out/           # a whole sheet, or a lane-split sample
   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --exact   # no speed levers at all

Flags shared by both modes, beyond input and output:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Flag
     - Effect
   * - ``--organism``
     - Reference organism. Default ``human``. Read per sample from a sheet's ``species`` column.
   * - ``--threads``
     - Search threads. ``0`` uses all cores.
   * - ``--exact``
     - Turn every speed lever off. Use when you want the unaccelerated result.
   * - ``--assemble`` / ``--no-assemble``
     - Stage 3 contig assembly for CDR3s no single read spans. On by default.
   * - ``--ec-mode``
     - Denoising preset: ``fast`` (default), ``accurate``, ``amplicon``, ``rnaseq``. See
       :doc:`error_correction`.
   * - ``--clonotype-key``
     - ``full`` (default: locus, V, J, junction) or ``junction``.
   * - ``--call-level``
     - ``allele`` (default) or ``gene`` — drop the allele suffix before the clonotype key.
   * - ``--shm``
     - ``framework`` (default) or ``both``. See :doc:`shm`.
   * - ``--isotype`` / ``--no-isotype``
     - IGH constant-region class per clonotype. On by default.
   * - ``--map-d`` / ``--no-map-d``
     - D and tandem D-D alignment into the junction. On by default. See :doc:`d_segments`.
   * - ``--min-junction-q``
     - Phred gate on the base that discriminates a clonotype from its parent. Requires
       ``map --junction-quality``.
   * - ``--cell-from``, ``--cell-regex``
     - Lift a cell barcode out of ``sequence_id`` into ``cell_id``.
   * - ``--limit``
     - Stop after N reads. Useful for a fast look at ``fast_fraction`` before committing.

``--indel-rescue`` is accepted by ``amplicon`` only, and refused rather than ignored elsewhere.

The three stages
----------------

The modes run these for you. Run them separately to inspect one, replace one, or shard one.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Command
     - What it does
   * - ``arda map``
     - **Stage 1, per read.** Stream FASTQ, keep the reads that map to a receptor scaffold, write
       them as AIRR. This is where all per-read work happens and the only stage that shards.
   * - ``arda assemble``
     - **Stage 3, global.** Reconstruct clonotypes whose CDR3 no single 100–150 bp read spans, by
       greedy overlap extension anchored on Stage 1's per-read ``cdr3_start``.
   * - ``arda correct``
     - **Stage 2, global.** Aggregate reads into clonotypes and collapse sequencing-error
       variants of a junction onto their parent.

.. code-block:: bash

   arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --report map.json
   arda assemble -i mapped.airr.tsv -o assembled.airr.tsv
   arda correct  -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv

.. important::

   ``--extra-airr`` is what folds Stage 3 back into the clonotype table. Without it the assembled
   reads are written and then discarded. The modes wire it for you.

Notable ``arda map`` flags:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Flag
     - Effect
   * - ``--junction-quality``
     - Emit a Phred+33 string over exactly the bases of ``junction``. Required by
       ``correct --min-junction-q`` and by the junction-quality QC metrics.
   * - ``--mutation-quality``
     - Emit the Phred behind each ``v_mutations`` / ``j_mutations`` entry, comma-joined and
       one-for-one.
   * - ``--reconstruct``
     - Merge each overlapping mate pair into one fragment, resolving overlap mismatches by the
       higher-Phred base.
   * - ``--two-pass``, ``--fast-segments``, ``--v-only-on-segment``,
       ``--prefilter``, ``--indel-rescue``
     - The five speed levers, individually, for A/B work. All off by default. See
       :ref:`speed levers <usage-speed-levers>`.
   * - ``--min-score``, ``--max-seqs``, ``--kmer``, ``--sensitivity``
     - Search parameters passed through to MMseqs2.
   * - ``--report PATH``
     - Write the stage's run report as JSON.

Annotating sequences
--------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Command
     - What it does
   * - ``arda annotate``
     - Annotate FR/CDR regions for sequences you already have and write an AIRR TSV. Streamed and
       memory-bounded. ``--seqtype nt`` or ``aa``; ``--strand both`` by default.
   * - ``arda markup``
     - Mark up and repair bare ``(junction_aa, V, J)`` records that have no read behind them — a
       VDJdb row. Emits a vdjdb-style ``cdr3fix`` report.
   * - ``arda shm``
     - Recount somatic hypermutation outside the junction on a table already written. Needs no
       reference and no re-map.
   * - ``arda resolve-ties``
     - Widen ``v_call`` / ``j_call`` to every germline the read's alignment cannot rule out, or
       narrow them to a donor's genotype with ``--genotype``.
   * - ``arda genotype``
     - Infer which V alleles a donor carries, from reads arda has already mapped. See
       :doc:`genotype`.
   * - ``arda igblast``
     - Run IgBLAST across all annotatable loci and emit AIRR, as a gold standard to compare
       against. IgBLAST is fetched on first use.

.. code-block:: bash

   arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
   arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa
   arda markup -i junctions.tsv -o marked.tsv --report -
   arda shm -i mapped.airr.tsv -o rescoped.airr.tsv
   arda resolve-ties -i mapped.airr.tsv -o widened.airr.tsv --loci IGK,IGL
   arda genotype -i mapped.airr.tsv -o donor.genotype.tsv --loci TRB
   arda igblast -i reads.fastq -o truth.airr.tsv

Quality control
---------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Command
     - What it does
   * - ``arda stats``
     - Per-sample run QC as one long-format TSV — reads, clonotypes, chains, genes, candidate
       alleles and the distributions. Every mode writes this automatically; the command reproduces
       it from any subset of a run's artifacts.
   * - ``arda qc batch``
     - Combine every sample's ``stats.tsv`` in a directory into one cohort table, with each
       metric's group median, MAD and robust z. Reads only ``stats.tsv`` files.
   * - ``arda qc report``
     - Render a QC JSON as one self-contained interactive HTML page, with no external reference of
       any kind.

.. code-block:: bash

   arda stats -i SAMPLE.airr.tsv -c SAMPLE.clones.tsv -r SAMPLE.arda.json -o SAMPLE.stats.tsv
   arda qc batch  -d results/ -o results/batch --samples sheet.tsv
   arda qc report -i results/batch.qc.json -o results/batch.qc.html

Details, including what each scope means and why no threshold ships, in :doc:`qc`.

Running at scale
----------------

``arda cluster`` holds every sharding and SLURM helper. A sharded run is byte-identical to a
single-node one: Stage 1 shards, Stages 2–3 run once over the merged output.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - What it does
   * - ``arda cluster submit``
     - Write, and optionally submit, a SLURM script chaining split → array ``map`` → reduce with
       ``afterok`` dependencies. The usual entry point for one large paired FASTQ.
   * - ``arda cluster submit-samples``
     - The same for a whole sample sheet, one array task per read group.
   * - ``arda cluster plan``
     - Print the read-group work units a scheduler would allocate workers to, without submitting.
   * - ``arda cluster split``
     - Split paired FASTQ into contiguous blocks of read pairs.
   * - ``arda cluster reduce``
     - Merge a sharded Stage 1, then run Stages 2–3 once over the whole thing.
   * - ``arda cluster merge``
     - Concatenate per-shard AIRR TSVs into one, with a single header.
   * - ``arda cluster split-fasta``, ``arda cluster submit-fasta``
     - The single-end / FASTA chain. These drop quality and separate mates, so paired RNA-seq must
       not use them.

.. code-block:: bash

   arda cluster submit --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE --shards 32 \
       --work-dir work/ -d results/ --threads 8 --time 02:00:00 --mem 16G --partition medium

See :doc:`cluster` for resource sizing, and :doc:`pipeline_integration` for the Nextflow module
and Snakemake workflow.

The reference database
----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Command
     - What it does
   * - ``arda info``
     - Print the resolved reference path, cache location and external tool availability. Run this
       first when anything looks wrong.
   * - ``arda export-ref``
     - Export reference sequences with their FR/CDR markup, in three kinds (``scaffolds``,
       ``segments``, ``anchors``) × four formats (``tsv``, ``fasta``, ``gff3``, ``airr``).
   * - ``arda build-db``
     - Rebuild the curated reference from IMGT. Needs IgBLAST. Not required for annotation.
   * - ``arda build-index``
     - Rebuild the precompiled MMseqs2 indexes shipped under ``database/``.

.. code-block:: bash

   arda info
   arda export-ref --kind scaffolds --locus TRB --format gff3 -o trb.gff3
   arda export-ref --kind segments  --locus TRB --format fasta
   arda build-db --organism all

See :doc:`reference_export` and :doc:`reference_build`.

Model estimation
----------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Command
     - What it does
   * - ``arda scenarios``
     - Estimate a generative model of V(D)J recombination from nucleotide junctions — trimming and
       insertion distributions, ``P(D | J)`` — by EM over the recombination scenario set. See
       :doc:`scenarios`.
   * - ``arda shm-model``
     - Fit ``P(substitution | 5-mer germline context, region)`` from arda's own mutation calls.
       See :doc:`shm`.

.. code-block:: bash

   arda scenarios -i clones.tsv -o d_prior.tsv
   arda shm-model -i mapped.airr.tsv -o shm.tsv --locus IGH

Streams and exit behaviour
--------------------------

* Progress, stage lines and warnings go to **stderr**.
* The paths of the files written go to **stdout**, one per line, so
  ``$(arda rnaseq ... | head -1)`` is a usable shell idiom.
* A missing input, an unknown organism, an unparseable sample sheet or an allele name absent from
  the reference **raises** rather than degrading to a partial answer.
