Quality control
===============

Every run writes ``<prefix>.stats.tsv`` and ``<prefix>.stats.json`` alongside its outputs: the
numbers that decide whether a sample is usable, derived from the artifacts the run already
produced and **without re-reading the FASTQ**. ``arda rnaseq``, ``arda amplicon`` and
``arda cells`` all write them, in the same scopes, so one cohort table holds every kind of
sample. ``arda stats`` builds the same table from any subset of those artifacts, so it also runs
on a bare ``arda annotate`` output.

Three commands, and the split between them is *reduce this run* against *compare these runs*:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - command
     - what
   * - ``arda stats``
     - one sample's QC table. Written for you by every run; run it by hand for an ``annotate``
       output or to change the allele thresholds.
   * - ``arda qc batch``
     - many samples' tables as one cohort view, with each metric's group median, MAD and
       robust z.
   * - ``arda qc report``
     - a cohort (or one sample) as a single self-contained, interactive HTML page.

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
``--cells``
   an ``arda cells`` output **prefix**, not a file. Its report and ``.chains.tsv`` become the
   same ``sample`` / ``chain`` / ``*_gene`` / ``junction_aa_len`` scopes a bulk run writes.
``--json``
   also write the rows as nested, typed JSON. A run writes ``<prefix>.stats.json`` without
   being asked.

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
   * - ``junction_aa_len``
     - ``IGH:17``
     - ``reads``, ``clonotypes`` (or ``chains``) at each junction length in residues
   * - ``read_len``
     - ``IGH:90``
     - ``reads`` at each aligned read length, 10-nt buckets
   * - ``clone_size``
     - ``IGH:8``
     - ``clonotypes`` and ``reads`` in each power-of-two clonotype size
   * - ``isotype``
     - ``IGH:IGHG``
     - ``clonotypes`` and ``reads`` per constant-region class
   * - ``chain_support``
     - ``TRB:4``
     - ``chains`` at each power-of-two molecule count (``arda cells`` only)

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
has no ``junction_quality_mean`` row rather than a zero that reads like a terrible library. The
same holds for an aggregation with nothing to aggregate: if no IGH read spans both anchors there
is no ``chain IGH junction_nt_min`` row, rather than a blank that casts to zero.

Distributions
-------------

min, max and mean cannot show a **shape**, and shape is most of the diagnosis. A bimodal junction
length is two primer sets in one tube; a read-length cliff is an adapter left on; a clone-size
distribution with no singletons is a library amplified before it was sequenced. Each of those
reads as an unremarkable mean.

Five scopes carry the distributions, all keyed ``locus:bucket`` so one parser reads them all, and
all **sparse** — only occupied buckets get a row, the same rule ``v_gene`` uses. An amplicon
occupies perhaps 30 of 200 possible junction lengths.

.. code-block:: text

   $ awk -F'\t' '$1=="clone_size" && $3=="clonotypes"' S0.stats.tsv
   clone_size  IGH:1  clonotypes  2
   clone_size  IGH:2  clonotypes  1
   clone_size  IGK:1  clonotypes  7
   clone_size  IGK:2  clonotypes  3

   $ awk -F'\t' '$1=="isotype"' S0.stats.tsv
   isotype     IGH:IGHG  clonotypes  1
   isotype     IGH:IGHG  reads       1
   isotype     IGK:IGKC  clonotypes  4
   isotype     IGK:IGKC  reads       8

``clone_size`` buckets are powers of two keyed by their lower bound, and come weighted both by
clonotypes and by reads — the read-weighted form is the one that shows over-amplification.
``read_len`` is the length arda actually **aligned**, which is not the run report's
``read_length_mean``: that one is measured over every read before mapping.

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

A cohort: ``arda qc batch``
---------------------------

The long format exists so a metric can be joined across samples. ``arda qc batch`` performs that
join, reading **only** the per-sample ``*.stats.tsv`` — never an AIRR, never a clonotype table —
so a thousand-sample cohort is one concatenation and runs on results copied off a cluster with
the bulk data left behind.

.. code-block:: bash

   arda qc batch -d results/ -o results/batch --samples sheet.tsv --report results/batch.html

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - file
     - what
   * - ``<prefix>.qc.tsv``
     - long: ``sample``, ``project``, ``batch`` in front of the four ``stats`` columns, plus
       ``median``, ``mad``, ``z`` and ``outlier``
   * - ``<prefix>.qc.wide.tsv``
     - one row per sample, ``sample`` scope only — the table you read
   * - ``<prefix>.qc.json``
     - typed and nested, what ``arda qc report`` inlines

``--samples`` is read **only** for two optional sheet columns, ``project`` and ``batch`` (see
:ref:`samples-sheet`). They are labels: nothing in the pipeline looks at them and no output
changes because of them. They name the group a sample is compared within.

.. code-block:: text

   $ cut -f1-3,5,7 results/batch.qc.wide.tsv | column -t
   sample  project  batch  reads  clonotypes
   S0      TRIAL9   RUN1   210    19
   S1      TRIAL9   RUN1   103    16
   S2      TRIAL9   RUN1   24     3

.. important::

   **Flags, never filters, and no shipped threshold.** There is no pass/fail column and no
   constant anywhere saying what a good ``mapped_fraction`` is — that depends on the library, the
   organism and the depth, and a number that looked calibrated would be worse than none. What a
   batch *can* say is whether a sample looks like the batch it came in with. So within each
   ``(project, batch)`` group, every numeric metric carries its group's ``median``, its ``mad``,
   and a robust z, ``0.6745 · (x − median) / mad``; ``outlier`` is 1 at ``|z| ≥ 3.5``
   (Iglewicz–Hoaglin). You decide what to do about it.

   A group of fewer than **5** samples gets a median and **no z** — MAD over four points is not a
   scale estimate, and a z derived from one would flag whichever sample happened to sit furthest
   from the middle. A metric whose MAD is 0 gets none either. Non-numeric values (a version, a
   file list) never acquire a median.

.. code-block:: text

   $ awk -F'\t' '$11==1 && $4=="sample"' results/batch.qc.tsv | cut -f1,6,7,8,10
   S0  clonotype_junction_aa_max  21  16  6.745
   S4  clonotype_junction_aa_max  21  16  6.745

.. note::

   Repertoire biology is deliberately **not** here: no diversity, clonality, rarefaction, overlap
   or cross-sample clonotype matching. Those belong to ``vdjtools``, which takes arda as a
   dependency. This answers "did this run work, and is this sample like its batch".

The dashboard: ``arda qc report``
---------------------------------

.. code-block:: bash

   arda qc report -i results/batch.qc.json -o results/batch.qc.html
   arda qc report -i results/S0.stats.json -o S0.qc.html      # one sample is a cohort of one

One HTML file with the data inlined and the charts drawn by plain JavaScript. It has **no
external reference of any kind** — no CDN, no bundle, no fetch — so it opens on an air-gapped
login node, off a USB stick, or as an email attachment, and it keeps working after the results
directory is gone. arda gains no plotting dependency for it, which is the same stance as
``arda cells --plot``: the gnuplot script is written whether or not gnuplot exists.

Four panels, all filtered together by project, batch and sample name:

**Samples**
   one row per sample, sortable by any column, cells outlined where ``|z| ≥ 3.5``. Nineteen
   metrics by default; *all metrics* shows every one the cohort carries.
**One metric across the cohort**
   any numeric metric as a bar chart, with each group's median drawn across it.
**Distributions**
   junction length, read length, clone size, chain support, isotype and gene usage, one line per
   sample, per locus, as a fraction of each sample's own total — so one deep sample does not
   flatten the rest. Click a legend entry to drop a curve.
**Provenance**
   arda and mmseqs versions, the reference and its path, threads, ``min_score`` and the input
   files, per sample.

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
