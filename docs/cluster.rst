Running on a cluster (SLURM)
============================

arda ships three commands for splitting a large input across a scheduler and putting the pieces
back together: ``arda split``, ``arda slurm`` and ``arda merge``. They are a thin layer over the
same per-shard CLI you would run by hand, so nothing about the result depends on the scheduler.

.. important::

   **Sharding must not change the answer.** arda's Stage 2 sorts its groups on a total key and
   sorts the reads inside each group, so a sharded run and a single-node run over the same input
   produce a **byte-identical** clonotype table. That property is not free — it exists because
   ``polars``' ``group_by`` is a multithreaded hash aggregation, and without the sorts three
   different runs of the same input produced three different tables. Anything you add around these
   commands must preserve it: shard **contiguously**, merge **once**, and run Stage 2 **globally**
   rather than per shard.

One command for the whole chain
-------------------------------

.. code-block:: bash

   arda slurm reads.fq annotated.tsv work/ \
       --shards 32 --threads 8 --time 02:00:00 --mem 16G --partition medium

That writes ``work/submit.sh``, which chains three steps with an ``afterok`` dependency so the
merge runs only once every shard has succeeded:

1. ``arda split`` — one cheap pass over the input;
2. ``sbatch --array=0-31`` — one ``arda annotate`` per shard;
3. ``arda merge`` — concatenate the per-shard AIRR TSVs under a single header.

Add ``--submit`` to submit it instead of only writing it. ``--arda-mmseqs`` exports
``ARDA_MMSEQS`` into the array tasks.

.. warning::

   **Pin the aligner.** An MMseqs2 index is only reusable by the release that built it. If array
   tasks resolve different ``mmseqs`` binaries, some will silently reject the precompiled reference
   index and rebuild a private cache — no error, and results that are not comparable across shards.
   Pass ``--arda-mmseqs /path/to/mmseqs`` or export ``ARDA_MMSEQS`` in your job script.

Splitting by hand
-----------------

.. code-block:: bash

   arda split reads.fq shards/ --shards 32               # single-end / amplicon, FASTA out
   # paired input, contiguous blocks of PAIRS, quality preserved:
   python -c "from arda.cluster import split_pairs; split_pairs('r1.fq','shards',shards=32,r2='r2.fq')"

.. warning::

   ``split`` and ``split_pairs`` are **not** interchangeable. ``split`` writes FASTA — which drops
   the quality string that ``merge_pair``'s per-base tie-break needs under ``--reconstruct`` — and
   round-robins *records*, which puts the two mates of one fragment in different shards. Use
   ``split_pairs`` for paired FASTQ: it writes contiguous blocks of read **pairs**, byte for byte.

Very large inputs: shard Stage 1, run Stage 2 once
--------------------------------------------------

For a full-depth library the pattern is:

.. code-block:: bash

   # Stage 1, per shard, in an array task
   arda map --r1 shards/shard_${SLURM_ARRAY_TASK_ID}_1.fq \
                   --r2 shards/shard_${SLURM_ARRAY_TASK_ID}_2.fq \
                   -o out/shard_${SLURM_ARRAY_TASK_ID}.airr.tsv \
                   --junction-quality --prefilter --threads ${SLURM_CPUS_PER_TASK}

   # then ONCE, over the merged per-read table
   arda merge out/*.airr.tsv all.airr.tsv
   arda correct -i all.airr.tsv -o clones.tsv --ec-mode rnaseq

.. important::

   **Stage 2 runs once, globally — never per shard.** Error correction compares a clonotype against
   its neighbours by abundance, so running it per shard asks the question against a fraction of the
   evidence and gets a different answer in each one. The same applies to ``assemble``.

Choosing the regime
-------------------

The two speed levers do **not** compose, and picking the wrong one is a slowdown reported as a
result:

.. list-table::
   :header-rows: 1
   :widths: 22 40 38

   * - library
     - flags
     - why
   * - amplicon / RepSeq
     - ``--two-pass --fast-segments --v-only-on-segment``
     - Primer-anchored reads hit both a V and a J (~85 %), which is what the two-pass needs.
   * - bulk RNA-seq
     - ``--prefilter``
     - Bulk reads land anywhere in a transcript and hit both only ~5 % of the time; the cost there
       is the k-mer scan, which is what the prefilter removes.

``--two-pass`` on its own is **slower** than the default on both bulk and IGH amplicon.

Resource sizing
---------------

Right-size the request to the real read count — an oversized request can queue behind quota far
longer than it saves.

.. list-table::
   :header-rows: 1
   :widths: 30 22 22 26

   * - workload
     - threads
     - memory
     - note
   * - 100 k amplicon reads
     - 8
     - 4 GB
     - ~5 s wall in the amplicon regime
   * - 5 M bulk reads
     - 16
     - 16 GB
     - ~35 s wall with ``--prefilter``
   * - full-depth bulk (50–100 M)
     - 16–32
     - 32–64 GB
     - shard Stage 1 if wall time is the constraint

.. note::

   ``sacct``'s ``MaxRSS`` counts **page cache** on some clusters, so it scales with input size
   rather than with footprint — a run yielding zero clonotypes has been reported at 2.70 GB. Use
   ``/usr/bin/time %M`` when comparing tools, and ``arda.json``'s ``peak_rss_mb`` when comparing
   arda's own modes.

Prefilter threads saturate at 16 and **regress at 32** (4.99× → 4.35× on a Xeon Silver 4210R): it is
bound by memory bandwidth over a small index, not by CPU. A 32-thread run of that stage pays ~4.6 %
more than a 16-thread one.
