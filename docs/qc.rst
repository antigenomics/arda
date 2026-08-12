Run QC, verbosity and logging
=============================

Every ``arda rnaseq`` / ``arda amplicon`` run writes ``<prefix>.stats.tsv`` alongside its
outputs: the numbers that decide whether a sample is usable, derived from the artifacts the run
already produced and **without re-reading the FASTQ**. ``arda stats`` builds the same table from
any subset of those artifacts, so it also runs on a bare ``arda annotate`` output.

.. code-block:: bash

   arda stats -i SAMPLE.airr.tsv -c SAMPLE.clones.tsv -r SAMPLE.arda.json \
              --r1 R1.fq.gz --r2 R2.fq.gz -o SAMPLE.stats.tsv

Every input is optional and each contributes its own scopes:

``--airr`` (``-i``)
   Stage-1 or ``annotate`` AIRR — the per-read ``chain`` rows, per-gene read counts, and the
   candidate-allele shortlist.
``--clones`` (``-c``)
   the clonotype table — the per-chain clonotype rows, chimera counts, per-gene clonotypes.
``--report`` (``-r``)
   ``<prefix>.arda.json`` or a single-stage ``--report`` JSON. **The only source** of total and
   mapped reads, FASTQ size, read length, pairedness, threads, wall time and peak RSS: the AIRR
   holds the mapped subset, so its row count and ``sequence`` lengths describe receptor reads
   rather than the library.
``--r1`` / ``--r2``
   read only for their size on disk and to record that the library is paired. Use these when you
   have no run report.

The table
---------

Four columns — ``scope``, ``key``, ``metric``, ``value`` — one value per cell.

.. list-table::
   :header-rows: 1
   :widths: 18 24 58

   * - scope
     - key
     - what
   * - ``run``
     - ``map`` / ``correct`` / ``assemble``
     - the run report, flattened verbatim: reads, ``input_bytes``, ``read_length_*``, ``paired``,
       ``threads``, ``wall_seconds``, ``peak_rss_mb``, ``per_locus.*``, the prefilter and
       segment-search accounting
   * - ``sample``
     - *(blank)*
     - library-wide totals, junction lengths and quality, SHM rate, V/J gene coverage,
       ``allele_candidates`` / ``shm_variants`` and their mean Phred
   * - ``chain``
     - ``TRB``, ``IGH``, …
     - per locus, **reads and clonotypes**: productive / non-functional, stop codons,
       out-of-frame, truncated junctions, min/max/mean junction length in nt and aa,
       junction quality, SHM rate, mutations per read, chimeras
   * - ``v_gene`` / ``j_gene``
     - ``TRBV19``
     - ``reads``, ``clonotypes``, ``reads_in_clonotypes`` per germline gene
   * - ``allele_candidate``
     - ``TRBV19*01:G45A``
     - ``reads``, ``allele_reads``, ``frequency``, ``mean_quality`` for a recurrent V mutation

.. code-block:: text

   $ awk -F'\t' '$1=="chain" && $2=="IGH"' SAMPLE.stats.tsv
   chain  IGH  reads                       104
   chain  IGH  reads_with_junction         5
   chain  IGH  reads_truncated_junction    1
   chain  IGH  junction_nt_min             42
   chain  IGH  junction_nt_max             63
   chain  IGH  junction_nt_mean            48.75
   chain  IGH  junction_quality_mean       35.6718
   chain  IGH  shm_rate                    0.0413
   chain  IGH  clonotypes_chimeric         2

Long, not wide, and deliberately: the metric set differs per scope (a gene has no junction
length, a chain has no allele frequency), so a wide table would be mostly empty cells. Long
format is what ``grep``, ``join`` and a per-metric plot across samples want.

A metric with **no input is omitted, never emitted as 0**. A run without ``--junction-quality``
has no ``junction_quality_mean`` row rather than a zero that reads like a terrible library.

Gene coverage
-------------

``v_gene_coverage_reads`` is the fraction of the organism's V genes seen on at least one read;
``..._multi`` is the fraction seen on more than one; ``..._clonotypes`` is the same over the
clonotype table. The reference universe comes from the shipped ``cdr3_anchors.tsv``, so coverage
is measured against the germline set arda actually maps to rather than a hand-kept list.

Quality columns
---------------

Two of the QC metrics need a column Stage 1 only writes when asked. Both are opt-in on ``map``,
because both append non-schema columns:

``--junction-quality``
   the read's Phred+33 string over exactly the bases of ``junction``, same orientation. Also what
   ``correct --min-junction-q`` gates on. See :doc:`error_correction`.
``--mutation-quality``
   ``v_mutation_quality`` / ``j_mutation_quality``: the Phred of the read base behind each entry
   of ``v_mutations`` / ``j_mutations``, comma-joined, **one-for-one and in the same order**.

.. warning::

   The two encodings differ. ``junction_quality`` is raw Phred+33 *characters* (so it lines up
   byte-for-byte with ``junction``); ``v_mutation_quality`` is comma-joined *integers* (there is
   no string to line up with). Reading one as the other gives plausible numbers off by 33.

Alleles versus SHM
------------------

A novel allele, somatic hypermutation and a base miscall are the **same string** in the mutation
list. What separates them is how often the mutation recurs across the reads calling that allele —
a germline the reference does not carry is in essentially every one of them, while hypermutation
is per-clone — and how good the base is. ``arda stats`` reports both, per variant:

.. code-block:: text

   allele_candidate  IGHV3-53*01:G257C  reads          3
   allele_candidate  IGHV3-53*01:G257C  allele_reads   4
   allele_candidate  IGHV3-53*01:G257C  frequency      0.75
   allele_candidate  IGHV3-53*01:G257C  mean_quality   33

.. warning::

   This is a **shortlist to look at, never a genotype call** — arda does not genotype. The
   thresholds are exposed (``--allele-min-frac``, default 0.5; ``--allele-min-reads``, default 10)
   precisely so the number can be re-derived rather than trusted.

Likewise the chimera, non-functional and stop-codon counts are **flags, not filters**: nothing in
``stats`` removes a row from any output, and the chimera signature cannot separate a true PCR
template-switch artefact from two real clones sharing a prefix and a suffix (see
``correct --flag-chimeras``).

Verbosity and logging
---------------------

Three **global** options, placed before the subcommand:

.. code-block:: bash

   arda -v --log-file run.log amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/
   arda -q rnaseq --r1 R1.fq.gz -p SAMPLE -d out/          # warnings and errors only

``(default)``
   the stage lines plus a throttled progress line — reads seen, reads mapped, reads/s and peak
   RSS, at most one every 30 s. Time is the right axis, not chunk count: a bulk sample flushes
   hundreds of chunks and chunk wall time varies ~50× with the receptor fraction.
``-v`` / ``--verbose``
   DEBUG, with the level and module name on each line. Repeatable.
``-q`` / ``--quiet``
   warnings and errors only.
``--log-file PATH``
   **always DEBUG whatever the console level is**, with a timestamp and the process peak RSS on
   every line. ``-q`` does not silence it — a quiet cluster job should still leave a full record.

.. code-block:: text

   2026-08-12 16:26:42 DEBUG   arda                        18.9 MB  arda 2.20.0 | python 3.12.13 | darwin | 16 cores | pid 52823
   2026-08-12 16:26:44 INFO    arda.rnaseq.map            187.2 MB  map: 453/1,320 reads mapped (34.32 %) in 1.4 s (927 reads/s), peak 187 MB

Streams
-------

**Progress goes to stderr; results go to stdout.** A mode run prints its output paths one per
line and nothing else on stdout, so ``arda amplicon ... | tail -1`` and ``$(arda map ...)`` are
usable and ``arda export-ref ... > out.tsv`` cannot interleave a progress line into the data.

Peak RSS
--------

``peak_rss_mb`` is the **whole-process** (plus reaped children) high-water mark as of the end of
that stage, and it is monotone: ``getrusage`` reports high-water marks only and offers no
per-stage reset, so a stage cannot be charged its own peak when all three run in one process.
``rss_gain_mb`` is how much that stage raised the mark. For per-stage attribution, run
``arda map`` / ``assemble`` / ``correct`` separately.

Children are included on purpose: 92 % of a ``map`` run's wall time is inside the ``mmseqs``
subprocess, whose nucleotide prefilter allocates the ``4**k`` index table that dominates the
footprint (``--kmer`` is the memory knob). Reporting the Python process alone understated peak
RSS by roughly an order of magnitude.
