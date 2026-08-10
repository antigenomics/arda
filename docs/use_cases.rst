Use cases and analysis guidelines
=================================

What to run, and what the numbers mean afterwards. Each section is a real question with a
configuration that answers it and a read-out that does not overstate what was measured.

.. contents::
   :local:
   :depth: 1

Bulk RNA-seq: extract a repertoire from a transcriptome
-------------------------------------------------------

The primary case. A bulk library is 0.02–3 % receptor, so almost all the work is proving that
reads are *not* receptor reads — which is what the mode's ``--prefilter`` removes (1.99×, at a
measured ~0.15 % of mapped reads; ``--exact`` turns it off). ``--ec-mode rnaseq`` is this mode's
denoising default.

.. code-block:: bash

   arda rnaseq --r1 R1.fq.gz --r2 R2.fq.gz -d out -p sample --threads 16

Read-outs, and what each is worth:

``clones.tsv``
   One row per clonotype. ``duplicate_count`` is **reads encompassing the junction** (coverage, the
   true expression estimate); ``consensus_count`` is distinct fragments. Use ``consensus_count``
   when you need molecules, ``duplicate_count`` when you need expression.
``arda.json``
   ``mapped_reads / total_reads`` is the receptor fraction. On bulk that is normally 0.02–3 %; a
   value far above it on a non-lymphoid tissue is a reason to look at the input, not a result.

.. warning::

   **Depth is not a nicety here.** A 200,000-read subsample of a bulk library yields a handful of
   clonotypes — measured: 1 for GM12878, 3 for Jurkat — and no diversity, overlap or expansion
   statistic computed on that means anything. Full depth on the same GM12878 library gives a real
   repertoire. If you must subsample, say so beside every number.

Targeted amplicon / RepSeq
--------------------------

.. code-block:: bash

   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -d out -p sample --threads 16

The mode carries ``--two-pass --fast-segments --v-only-on-segment`` and ``--ec-mode amplicon``.
Those are **not** interchangeable with the bulk ones and do not compose, which is why the regime is
the command name; see :doc:`cluster`. On an IGH RepSeq amplicon the amplicon mode is **4.15×**
faster than ``--exact`` at 2.7× less memory.

Monoclonal QC: is my cell line what it says it is?
---------------------------------------------------

A monoclonal line has one productive rearrangement per expressed locus, so the read-out is
**purity** — the dominant clonotype's share of its locus's reads — and, if the clone is published,
**how many reads land on it**.

.. code-block:: bash

   arda rnaseq --r1 R1.fq.gz -d out -p Jurkat \
       --ec-mode amplicon --clonotype-key junction

Measured on Jurkat (14,531 junction-bearing reads):

.. list-table::
   :header-rows: 1
   :widths: 34 16 16 18 16

   * - configuration
     - clonotypes
     - reads
     - TRB purity
     - reads on the 2 published clones
   * - ``fast`` (default)
     - 90
     - 14,531
     - .98963
     - 14,177 (.9756)
   * - ``amplicon`` + ``junction``
     - **10**
     - **14,531**
     - **.99990**
     - **14,313 (.9850)**

.. important::

   **A clonotype count falling is the job; a read count falling is a bug.** Every arda denoising
   mode *moves* reads onto a parent and never discards them, so ``duplicate_count`` summed over the
   table is invariant. If it drops when you change ``--ec-mode``, that is a defect — report it.

.. warning::

   Do **not** reach for "drop clonotypes with 1 read" instead. Measured on the same library it
   loses 82 reads outright and plateaus at .99398 purity, because the largest error class *by
   reads* is call splits carrying 33 and 30 reads — no abundance rule separates those. Raising the
   threshold from 1 to 2 buys **0.00000** purity and loses 2 more reads.

Negative controls
-----------------

K562 and HepG2 carry no V(D)J recombinase, so **any** clonotype is a false positive by
construction. Run them exactly as your real samples and expect zero.

.. note::

   Passing them validates *cross-lineage specificity* and nothing else. In this project four
   different reference configurations — including one with a ``--min-score`` low enough to take
   precision from 94.3 % to 65.5 % — all yielded zero clonotypes on both. A negative control that
   passes is not evidence that a threshold is right.

Low-frequency variants: spike-ins, MRD, minor clones
-----------------------------------------------------

.. code-block:: bash

   arda correct -i s1.airr.tsv -o clones.tsv --ec-mode accurate --error-rate 1e-5

Two things decide whether a rare real variant survives:

``--error-rate``
   At the shipped ``1e-3`` the abundance model erases both published MIGEC spike-in variants. At
   ``1e-5`` both are recovered exactly. ⛔ This is a **per-library calibration**, not a default
   change: the right value depends on the library's error-cloud abundance ratio.
``--ec-mode accurate``
   Adds a Phred gate on the base that discriminates a clonotype from its parent — evidence the
   abundance model does not have. A real low-frequency variant is a *good read of a rare molecule*:
   the MIGEC variants sit at median mean junction Phred 37.2–37.5 against the parent's 37.8.

.. warning::

   At signal-to-noise ≈ 1 **no** abundance threshold can separate a real variant from PCR error.
   On the MIGEC data the published V2 variant is *less* abundant than the worst 2-substitution PCR
   error (V2/Err2 = 0.28 on raw reads). That is not a tool limitation; it is why UMI consensus
   exists. If your variant is at that level, the answer is chemistry, not parameters.

Somatic hypermutation
---------------------

``v_mutations`` / ``j_mutations`` record SHM in germline coordinates; see :doc:`shm`.

.. warning::

   **Local zero-loss is not zero-loss.** Three times in this project a clean result turned out to
   be a property of the test library's SHM level. Validate anything SHM-sensitive on a
   hypermutated library (median V identity ~92 %), not on a naive one — and report both, because a
   change that helps a hypermutated repertoire can cost a naive one.

Comparing arda against another tool
------------------------------------

.. important::

   **Name the stage.** TRUST4's read and false-positive counts are scored at candidate
   *extraction*; arda's are its final post-``--min-score`` clonotype table. Those are different
   quantities and differencing them produces a number that means nothing.

   **Benchmark every tool at its best configuration.** MiXCR's shipped ``rna-seq`` preset discards
   reads it aligned (``minSumScore=200``), understating its recall roughly ninefold. Measured on a
   200 k Jurkat subsample: 154 aligned reads at the default preset, **5,136** with
   ``-OallowNoCDR3PartAlignments=true -OminSumScore=40``, at the same wall time.

   **Compare at gene level across tools.** MiXCR suffixes every allele ``*00``, so an allele-level
   comparison scores it zero by construction. And treat an ambiguous call as an abstention:
   comparing only the first element of a comma-joined tie list scores two tools that agree the
   answer is ambiguous *between the same genes* as a disagreement.

   **Give each call metric its own denominator.** ``j_gene`` agreement divided by all truth reads
   measures how many reads reach the J — a property of the library. On one IGH library that read
   .4678; over the reads where both the truth and the tool named a J it is **.9418**.

The one test that needs no external truth
------------------------------------------

If you have the same material sequenced two ways — a deep targeted amplicon and bulk RNA-seq — the
strongest validation available is that they agree:

* what fraction of the amplicon repertoire is recovered from RNA-seq, **weighted by amplicon
  abundance** (did the tool find the clones that are actually there?);
* rank correlation of abundances on the shared clones — Spearman, not Pearson, because the two
  protocols have completely different amplification biases;
* the RNA-seq-only rate, as a **candidate** false-positive rate and never as an FP count: the
  amplicon covers only its target loci and misses low expressers.
