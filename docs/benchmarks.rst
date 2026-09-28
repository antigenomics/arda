Benchmarks
==========

Every performance and accuracy figure quoted anywhere in arda's documentation is collected here
with the method behind it. Each table names what was run, on what, how many observations it rests
on, and what the comparison is and is not entitled to claim.

.. contents::
   :local:
   :depth: 2

How these are measured
----------------------

**Comparisons are end to end.** Each leg of a head-to-head runs a whole pipeline and emits
clonotypes. The only latitude a tool gets is picking its best-fitting preset for the library type.
A stage-against-stage ratio is not reported, because arda's AIRR-emitting stage and another tool's
non-AIRR-emitting stage do different amounts of work.

**Legs alternate and are replicated.** Three repetitions per leg, medians reported, all legs of a
tier on the same staged input.

**Coverage is printed before any rate.** A per-tool inner join hands each tool its own
denominator, so a truth read a tool emitted nothing for vanishes instead of counting as a miss.
Where that matters, both denominators are shown.

**Alleles are scored as tie lists.** IgBLAST and arda both return an ambiguous allele as a
comma-joined set; scoring that as a miss is a scoring artifact. Across 25 datasets the median is
``v_allele_exact`` .8328 against ``v_allele_resolved`` .9763.

The measurements live in a separate repository so they can be re-run independently of arda's
release history.

Speed and memory
----------------

TRA amplicon, 100,000 reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``arda amplicon`` · ``mixcr analyze generic-amplicon --rna`` · TRUST4 defaults. M3 Mac, 8 threads,
arda 2.27.0, MiXCR 4.7.0. TRUST4's rows count only its **complete** CDR3s (``C…[FW]``, no ``_`` or
``?``), so the clonotype and read columns mean the same thing in all three rows.

.. list-table::
   :header-rows: 1
   :widths: 30 14 14 16 14 12

   * - Pipeline
     - Wall (s)
     - CPU (s)
     - Peak RSS (MB)
     - Clonotypes
     - Reads in clonotypes
   * - MiXCR ``generic-amplicon``
     - **7.82**
     - 42.91
     - 3,052
     - 19,697
     - 42,712
   * - arda ``amplicon``
     - 14.40
     - **30.49**
     - 965
     - **19,841**
     - **43,503**
   * - TRUST4
     - 73.77
     - 117.96
     - **490**
     - 18,559
     - 37,688

MiXCR is 1.84× faster on wall clock on a primer-anchored amplicon; that is its regime. arda spends
1.41× less CPU and 3.16× less RSS, and returns the most clonotypes (+0.7 % over MiXCR, +6.9 % over
TRUST4) over the most reads (+1.9 %, +15.4 %). TRUST4 is 5.1× slower than arda on this arm.

Bulk RNA-seq, 660,000 pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

SRR5233639. ``arda rnaseq`` · ``mixcr analyze rna-seq`` · TRUST4 defaults. Same machine and thread
count.

.. list-table::
   :header-rows: 1
   :widths: 30 14 14 16 14 12

   * - Pipeline
     - Wall (s)
     - CPU (s)
     - Peak RSS (MB)
     - Clonotypes
     - Reads in clonotypes
   * - TRUST4
     - **13.68**
     - **54.53**
     - **460**
     - 1,941
     - 6,247
   * - arda ``rnaseq``
     - 21.95
     - 219.68
     - 1,033
     - **2,213**
     - **8,484**
   * - MiXCR ``rna-seq``
     - 30.23
     - 233.11
     - 2,849
     - 1,732
     - 4,288

On the regime arda exists for it finds the most: +27.8 % clonotypes and +97.9 % reads assigned
against MiXCR, +14.0 % and +35.8 % against TRUST4, at 1.38× MiXCR's wall clock and 2.76× less RSS
for comparable CPU. TRUST4 is genuinely 1.61× faster on wall at 4.0× less CPU here, and the
cheapest of the three on memory in both arms; it reaches 87.7 % of arda's clonotypes and 73.6 % of
its assigned reads. Walls were stable across three repetitions (arda 21.88–24.04, MiXCR
30.09–30.83, TRUST4 13.65–14.22).

IGH RepSeq amplicon: what the amplicon configuration is worth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

100,000 pairs, 32 threads, against the shipped one-pass default on hypermutated IGH.

.. list-table::
   :header-rows: 1
   :widths: 24 40 18 18

   * - Dataset
     - Configuration
     - Wall (s)
     - Peak RSS (MB)
   * - IGH_repertoire
     - one-pass default
     - 316.44
     - 4,018
   * - IGH_repertoire
     - ``--two-pass --fast-segments --v-only-on-segment``
     - **76.25**
     - **1,479**
   * - IGH_naive
     - one-pass default
     - 305.32
     - 3,736
   * - IGH_naive
     - ``--two-pass --fast-segments --v-only-on-segment``
     - **64.86**
     - **1,363**

4.15× and 4.71× on wall clock, at roughly 2.7× less memory.

Against TRUST4 on IGH amplicon, 32 threads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wall clock only. Every leg of a tier ran on the same staged input.

.. list-table::
   :header-rows: 1
   :widths: 28 16 28 28

   * - Dataset
     - Reads
     - arda ``amplicon`` wall (s)
     - TRUST4 wall (s)
   * - IGH_repertoire
     - 100,000
     - **201.05**
     - 615.04
   * - IGH_naive
     - 100,000
     - **136.51**
     - 359.66
   * - migec_exp1_TCR
     - 500,000
     - **316.25**
     - 423.96
   * - migec_exp1_IGH
     - 500,000
     - 223.77
     - **225.10**

A full-depth, hours-scale head-to-head is still cluster work and is not claimed here.

No regression across eight releases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

arda 2.18.0 against 2.27.0, same job, legs alternating, three repetitions, the same committed
reference and the same MMseqs2 binary — ``database/`` has not changed since ``v2.18.0``, which is
what makes the comparison valid.

.. list-table::
   :header-rows: 1
   :widths: 22 16 16 14 14 20

   * - Arm
     - 2.18.0 wall (s)
     - 2.27.0 wall (s)
     - 2.18.0 RSS
     - 2.27.0 RSS
     - Clonotypes / reads
   * - amplicon 100 k
     - 13.94
     - **13.66**
     - 939 MB
     - 939 MB
     - 19,841 / 43,503
   * - bulk 660 k pairs
     - **21.45**
     - 21.74
     - 1,033 MB
     - 1,033 MB
     - 2,213 / 8,484

Bulk is 1.4 % slower on wall and 0.8 % higher on CPU, inside the repetition spread (2.18.0
21.06–22.26, 2.27.0 21.52–22.07) and accounted for by the per-run QC stage 2.20.0 added; amplicon
is 2.0 % faster. Both arms are **byte-identical** between the two versions, checked as a call
digest over ``(locus, v_call, j_call, junction, duplicate_count)`` rather than a row count.

Memory
~~~~~~

arda is CPU-bound, and large FASTQ is streamed in bounded chunks with a background reader
prefetching the next chunk. **Mapping is flat at 300–650 MB at any read depth.**

Peak RSS tracks repertoire **richness**, not read count, because Stage 3 holds the clone set:

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Sample
     - Reads
     - Clonotypes
     - Stage-3 ``correct`` peak RSS
   * - B-cell-rich tumour
     - 105 M
     - 28,444
     - **2,071.7 MB**
   * - colder sample, more reads
     - 139 M
     - —
     - 549 MB

Budget about 4 GB, and size a SLURM ``--mem`` from Stage 3 rather than Stage 1.

Scaling against IgBLAST (synthetic)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   These two tables are **synthetic** — generated human IGH sequences, not a real library — from
   ``scripts/bench_vs_igblast.py`` and ``scripts/bench_prefilter.py``. They measure scaling
   *shape*, not head-to-head standing. Use the tables above for standing.

16 threads.

.. list-table::
   :header-rows: 1
   :widths: 20 18 18 22 22

   * - Sequences
     - arda wall
     - arda rate
     - Speedup vs IgBLAST
     - Region concordance
   * - 10,000
     - 5.5 s
     - ~1.8 k/s
     - 4.4×
     - 98.9 %
   * - 50,000
     - 16 s
     - ~3.0 k/s
     - 7.3×
     -
   * - 100,000
     - 30 s
     - ~3.3 k/s
     - 7.9×
     -

Bulk RNA-seq is faster *per read* than amplicon, because MMseqs2 prefilters by k-mer matching and
reads carrying no receptor k-mer are rejected before alignment. 150 nt reads, 16 threads:

.. list-table::
   :header-rows: 1
   :widths: 40 30

   * - Receptor content
     - Throughput
   * - 100 % (amplicon)
     - ~5.7 k reads/s
   * - 10 %
     - ~19 k reads/s
   * - 1 % (blood RNA-seq)
     - ~25 k reads/s

Accuracy
--------

Recall and precision on bulk RNA-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three-way comparison on the **same 16 datasets where every tool ran**: 5,273 real fragments,
Wilson 95 % confidence intervals. Best value bold; ties bold together.

.. list-table::
   :header-rows: 1
   :widths: 14 30 16 20 12 14 12

   * - Tool
     - Configuration
     - Recall
     - Precision (lower bound)
     - False positives
     - FP / 1 M reads
     - Peak RSS
   * - **arda**
     - shipped defaults
     - **0.986** [.982–.989]
     - **0.889** [.881–.897]
     - **70**
     - **21.9**
     - **254 MB**
   * - **TRUST4**
     - defaults
     - **0.987** [.984–.990]
     - 0.160 [.156–.164]
     - 22,185
     - 6,932.8
     - 306 MB
   * - MiXCR
     - ``-OallowNoCDR3PartAlignments=true -OminSumScore=40``
     - 0.964 [.958–.969]
     - 0.051 [.050–.052]
     - 91,229
     - 28,509.1
     - 1,213 MB
   * - MiXCR
     - ``align --preset rna-seq`` as shipped
     - 0.193 [.183–.204]
     - 0.565 [.542–.588]
     - 91
     - 28.4
     - 1,213 MB
   * - arda
     - ``--reconstruct`` (6 paired datasets)
     - 0.995
     - 0.905
     - 14
     - 11.7
     - 249 MB

**Recall is a statistical tie between arda and TRUST4** — a 6-fragment gap with overlapping
intervals. Neither should be quoted as having higher recall than the other. **Precision is not a
tie**: arda's lower bound, 0.881, sits above every competitor's *upper* bound.

Three caveats belong with this table.

MiXCR must be benchmarked at its best configuration, not its default: ``align --preset rna-seq``
as shipped gives 0.193 recall, and the two free options above take it to 0.964.

Precision here is a **lower bound** — the grey band ``30 ≤ v_score < 70`` is scored under neither
metric.

TRUST4's recall and false positives are scored on its **candidate extraction**, a different and
much less filtered stage than arda's post-``--min-score`` output. The two are not like for like on
false positives, which is why that column carries a per-million normalisation rather than a bare
ratio.

The J→C class
~~~~~~~~~~~~~

What explains the table above is the J→C class: 22 % of real fragments on bulk RNA-seq.

.. list-table::
   :header-rows: 1
   :widths: 34 32 34

   * - Tool
     - V-covered reads (n = 4,117)
     - J→C agreement (n = 1,156)
   * - arda
     - 0.9864
     - **0.9844**
   * - TRUST4
     - 0.9944
     - 0.9611
   * - MiXCR (recall configuration)
     - **0.9990**
     - 0.8382
   * - MiXCR (default)
     - 0.1482
     - 0.3538

On V-covered reads every tool is at or above 0.986. arda's overall recall rests on the 345
:term:`J+C scaffold` entries in its reference, which took this class from 0.0606 to 0.9844 and
overall recall from 0.78 to 0.986.

Reported as *agreement* rather than recall: arda now ships C scaffolds, so the adjudicator is no
longer independent of it on this class.

Gene calls on a targeted amplicon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Against an IgBLAST truth on the TRA amplicon, 100,000 reads, arda 2.27.0 and MiXCR 4.7.0 at its
best amplicon preset, both scored per read from the same truth file in the same job.

Coverage first. Over the 48,033 truth reads at ``v_score ≥ 70``, **arda emits a row for 48,030
(99.99 %) and MiXCR for 46,503 (96.81 %)**, so both denominators are shown.

.. list-table::
   :header-rows: 1
   :widths: 44 28 28

   * - Metric
     - arda 2.27.0
     - MiXCR 4.7.0
   * - *all 48,033 truth reads*
     -
     -
   * - ``v_gene`` recall
     - **.9867**
     - .9660
   * - ``v_gene`` precision
     - **.9996**
     - .9977
   * - ``j_gene`` recall
     - **.9892**
     - **.9892**
   * - ``j_gene`` precision
     - .9953
     - **.9996**
   * - ``junction`` recall (nt, exact)
     - .9473
     - **.9708**
   * - *common subset, 46,502 reads*
     -
     -
   * - ``v_gene`` recall
     - .9869
     - **.9977**
   * - ``v_gene`` precision
     - **.9997**
     - .9977
   * - ``j_gene`` recall
     - .9959
     - **.9996**
   * - ``junction`` recall (nt, exact)
     - .9533
     - **.9778**

The two views say different and equally true things. **On the reads it emits, MiXCR is the more
accurate caller.** **Over the whole library, arda recalls more V genes** — .9867 against .9660 —
because MiXCR emits nothing at all for 1,530 truth reads and arda for 3. arda's V calls are the
more precise of the two under either denominator: it declines rather than guessing.

Junction, stated the same way: of the 46,787 truth junctions, arda emits one for 94.81 % at .99919
precision among emitted, and MiXCR for 97.09 % at .99989. The 5.19 % arda declines are reads that
*have* an anchor pair and lost the projection — a specific known gap rather than a calling error.

MiXCR suffixes every allele ``*00``, that is, makes no allele call, so there is no ``v_allele``
comparator. arda's is .9868 resolved (.9461 by exact string, the four-point difference being
ambiguous-allele tie lists).

These five arda figures are unchanged from 2.11.1, fifteen releases earlier: ``v_gene`` recall
.9867, precision .9996, ``j_gene`` recall .9892, precision .9953, and junction precision among
emitted .99919 all reproduce to every digit published.

.. note::

   A V/J boundary disagreement *inside* the junction is not scored, because V(D)J recombination is
   probabilistic: exonuclease chew-back and N/P addition mean the V-end / N-D-N / J-start
   partition is often not identifiable from sequence alone, so the ground truth for it is unknown.
   What the table scores is the junction's **outer bounds**, the **gene and allele calls**, and
   whether a tool invents a junction it has no anchor for.

Region concordance with IgBLAST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On roughly 7,300 real GenBank mRNA records spanning all five organisms and their loci — committed,
gzipped test fixtures — region concordance with IgBLAST on productive records is **98–99.7 % per
organism**, and ``junction_aa`` / ``cdr3_aa`` match IgBLAST about 99 % while satisfying the AIRR
invariants exactly.

GenBank also carries genomic, partial and non-productive entries that confuse both tools; those are
excluded.

What an IG V call is worth, by read coverage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Measured against an ``arda igblast`` truth (``v_score ≥ 70``) on a human IGH 5'RACE library of
98,639 truth reads, stratified by how much V germline the **truth** alignment covers.

.. list-table:: ``v_gene`` recall by V germline span — human IGH 5'RACE, 98,639 truth reads
   :header-rows: 1
   :widths: 28 22 16 22

   * - V germline span
     - Truth reads
     - Share
     - ``v_gene`` recall
   * - < 60 nt
     - 5,802
     - 5.9 %
     - **.1170**
   * - 60–90 nt
     - 978
     - 1.0 %
     - .6748
   * - 90–120 nt
     - 4,224
     - 4.3 %
     - .5727
   * - 120–160 nt
     - 1,845
     - 1.9 %
     - .9507
   * - 160–200 nt
     - 33,903
     - 34.4 %
     - .9647
   * - ≥ 200 nt
     - 51,887
     - 52.6 %
     - **.9896**

At 200 nt or more of V germline arda scores .9896, which is the TRA amplicon's .9867 — there is no
IG-specific accuracy deficit at that coverage. The 5.9 % of reads carrying under 60 nt produce
about 56 % of every V miss on the library.

The same measurement on bulk RNA-seq, where reads land across the V rather than at its 3' end, over
two libraries and three loci (``SRR5233639`` / ``SRR5233640``, 660,000 read pairs each, 100 nt):

.. list-table:: ``v_gene`` recall by V germline span — human IG bulk RNA-seq
   :header-rows: 1
   :widths: 10 12 16 12 12 14 14

   * - Locus
     - Sample
     - Truth reads
     - All
     - < 60 nt
     - 60–90 nt
     - 90–120 nt
   * - IGH
     - 639
     - 5,458
     - .9331
     - .8472
     - .8592
     - **.9661**
   * - IGH
     - 640
     - 5,104
     - .9373
     - .8426
     - .8789
     - **.9689**
   * - IGK
     - 639
     - 4,534
     - .9828
     - .9587
     - .9805
     - **.9879**
   * - IGK
     - 640
     - 4,190
     - .9802
     - .9492
     - .9690
     - **.9897**
   * - IGL
     - 639
     - 3,142
     - .9494
     - .6966
     - .8932
     - **.9938**
   * - IGL
     - 640
     - 2,942
     - .9470
     - .6630
     - .9149
     - **.9903**

Monotonic in span for every locus and both samples, and the two samples agree within 0.5 points
everywhere. 100 nt reads cannot reach the 120 nt-and-above bins, which is why bulk tops out around
.97–.99 rather than at the .9896 the long bin reaches.

Four independent arms — two IGH libraries, a multiplex V-primer amplicon and a 5'RACE, each read
from both ends — put ``v_gene`` recall on the ``≥ 200 nt`` bin at .9891, .9935, .9923 and .9930.
The protocol name is not the risk factor: the same 5'RACE protocol scores .1170 on a short read and
.9645 at 251 nt. What matters is how much V the read carries, and from which end.

The largest single confusion on the 5'RACE library, ``IGHV1-69`` called as ``IGHV1-18``, is 4,745 of
roughly 9,080 misses, on reads where IgBLAST aligns germline positions 242–296 — 55 nt of the V's 3'
end. The two genes are only 0.9054 identical, so they are genuinely separable and this is missing
evidence rather than homology. IgBLAST lists 1.64 V genes per read on that library itself: its
answer there is a different tie-break on the same missing evidence, not a better one.

**Position beats length.** IGHV genes diverge in FR1/CDR1/CDR2 and are conserved through FR3 near
Cys104, so a short 5' alignment identifies the gene and a short 3' one does not — which is why the
same ``< 60 nt`` bin scores .1170 on 5'RACE, whose reads run C → J → V, against .8472 on bulk
RNA-seq of the same locus.

Somatic hypermutation does not order this. Stratified by IgBLAST's own ``v_identity`` the deficit
is non-monotonic: unmutated .9425, 97–99 % **.7518** (the worst of six bins), 92–95 % .9697 (the
best).

arda ships **no threshold** on this. A span gate was measured at 7.6 : 1 in favour on IGH 5'RACE
(.92585 → .97693) and a clear loss on bulk IGH, IGK and IGL (−3.24, −6.22 and −2.39 points of
recall for no precision gain), so ``v_germline_start`` and ``v_germline_end`` are emitted for you
to filter on instead.

Widening the V call, per locus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``arda resolve-ties`` measured against an ``arda igblast`` truth (``v_score ≥ 70``) on eleven arms
across five geometries, scored on **exact ``v_gene``-set agreement** rather than intersection.

.. list-table::
   :header-rows: 1
   :widths: 14 26 22 22 20

   * - Locus
     - Arms
     - Exact before
     - Exact after
     - Change
   * - **IGK**
     - 2 bulk
     - .5707 / .5558
     - **.9373 / .9308**
     - **+36.7 / +37.5 pt**
   * - **IGL**
     - 2 bulk
     - .8586 / .8577
     - **.9137 / .9115**
     - **+5.5 / +5.4 pt**
   * - **TRB**
     - 1 amplicon
     - .9299
     - **.9694**
     - **+4.0 pt**
   * - IGH
     - 6 (bulk, 5'RACE, multiplex)
     - .8255–.9755
     - .6860–.9608
     - **−0.6 to −14.9 pt**

The deciding quantity is the **ambiguity deficit**: IgBLAST's genes per read minus arda's, before
widening. arda's IGK call is 0.42 genes/read too narrow (1.18 against IgBLAST's 1.60), and closing
that gap *is* the 37-point win. arda's IGH call is already as wide as IgBLAST's (deficit 0.00–0.07
on every arm, 1.589 against 1.64 on a 5'RACE library), so there is nothing to recover and each arm
only pays the overshoot.

The cost is that fewer reads name a single gene: on IGK that share falls .8195 → .3883. Clonotype
cost is +2 of 362 on IGK, +3 of 23,559 on TRB, and zero on IGL. Intersection recall rises on all
eleven arms including IGH — only the exact-set number shows the IGH regression.

Single cell against Cell Ranger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On ``sc5p_v2_hs_PBMC_1k`` VDJ-T, **98.2 % of Cell Ranger's 943 CDR3s appear verbatim inside one of
arda's contigs**, with no germline reference used in the assembly at all. Method, the doublet
criteria and where the agreement curve breaks are in :doc:`singlecell`.

Denoising, measured
-------------------

Full derivations are in :doc:`error_correction`; these are the headline numbers.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Configuration
     - Measured effect
   * - ``--error-rate 1e-3`` (default)
     - Erases both published MIGEC spike-in variants. Not a defect: on the paper's own metric over
       raw reads, V1/Err1 = 1.35 and V2/Err2 = 0.28, so no abundance-based method at any threshold
       can separate them.
   * - ``--error-rate 1e-5``
     - Recovers both spike-in variants exactly. Alone, it cost 3.5 points of monoclonal purity
       (.99540 → .96034).
   * - ``--error-rate 1e-5 --ec-mode accurate``
     - Both variants kept **and** purity back to .99530. Spurious junctions 297 → 62, distinct
       error clonotypes 1,630 → 124. The published spike-ins read median Q 34–35 at the
       discriminating base; the error cloud around them reads median Q 24.
   * - ``--error-rate 1e-4``
     - Kept both variants while removing 72 % of real PCR errors on an independent error cloud.
   * - ``--ec-mode amplicon --clonotype-key junction``
     - Jurkat 90 → 10 clonotypes, TRB purity .98963 → .99990, reads unchanged at 14,531, with
       98.50 % of them on the two published clones.
   * - ``--ec-mode amplicon`` on hypermutated IGH
     - Removes 178 clonotypes carrying 179 reads, 177 of them singletons. On the matched naive
       library it removes zero — which is why the modes are off by default.
   * - ``--d-max-evalue 0.01``
     - D gene agreement with IgBLAST .9765 → .9985 on a TRB amplicon, at about one third the call
       rate.

.. important::

   Every denoising mode **moves** reads onto a parent and never discards them. The sum of
   ``duplicate_count`` is invariant across modes, and a clonotype with no qualifying parent keeps
   its reads and is reported as an orphan. On a polyclonal hypermutated repertoire a plain quality
   *filter* at the same threshold would strand 3.70 % of all junction-bearing reads with no parent
   to inherit them. If the read total moves when you change ``--ec-mode``, that is a defect worth
   reporting.

Costs of the speed levers
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Lever
     - Measured cost or gain
   * - ``--prefilter``
     - 1.99× on bulk. Costs ~0.15 % of mapped reads over 122 bulk datasets, up to 2.46 % on one
       library, concentrated in J→C and hypermutated IGH. Roughly 1.4 GB of RAM, and about 3 % on
       amplicon — where it is the wrong lever.
   * - ``--two-pass`` alone
     - **A loss**: 0.762× on bulk, 0.87× on an IGH amplicon. It pays only with
       ``--fast-segments``.
   * - ``--two-pass --fast-segments --v-only-on-segment``
     - 4.15× and 4.71× on two IGH RepSeq amplicons (table above).
   * - ``--adaptive``
     - 1.84× wall, 3.04× CPU, 1.34× RSS on 660 k real bulk pairs, read set preserved exactly — but
       ``junction`` moves on 93 of 35,795 rows, and the movement is one-directional: 89 to empty
       against 4 the other way, a net −85 against 3,856 reads carrying one. 91 of the 93 sit at
       90–150 bits, above the score-only trigger, so the trigger is not calibratable at this
       scale. Stays off.
   * - ``--indel-rescue``
     - Value tracks SHM load, so it remains a per-library call and never rides a preset. Accepted
       on ``amplicon`` only, and refused rather than ignored elsewhere.

See also
--------

* :doc:`usage` — which configuration your library needs, and the ``fast_fraction`` rule.
* :doc:`error_correction` — the denoising model and its calibration.
* :doc:`genotype` — why allele-level calls need read length.
* :doc:`how_it_works` — the mechanism these numbers are measuring.
