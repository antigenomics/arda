Usage
=====

arda turns receptor reads into a clonotype table in one command. Pick the mode that matches the
library; the mode name carries the speed configuration, so there are no flags to get right.

For copy-paste recipes — the one-liners, annotating sequences you already have, and the polars
snippets that turn a clonotype table into clonal fractions, gene usage, repertoire overlap and
cohort QC — see :doc:`examples`.

.. _choosing a mode:

The three modes
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 44 36

   * - command
     - the library it is for
     - what makes it fast
   * - ``arda rnaseq``
     - Bulk RNA-seq, where only 1–5 % of reads are receptor-derived.
     - ``--prefilter``: a C++ 16-mer screen answers "is this read receptor-derived" before
       MMseqs2 is invoked at all.
   * - ``arda amplicon``
     - Targeted RepSeq / amplicon, where reads span V into J.
     - ``--two-pass --fast-segments --v-only-on-segment``: a structural pass over a 924-target
       segment reference names the scaffold, so the 15,414-scaffold search mostly never runs.
   * - ``arda cells``
     - Single cell: **one UMI consensus per molecule**, barcode in the record name.
     - Reference-free per-cell assembly, then chain pairing and doublet calls — :doc:`singlecell`.

``rnaseq`` and ``amplicon`` take raw paired FASTQ, run ``map`` → ``assemble`` → ``correct`` in one
call, and write four files:

.. code-block:: bash

   arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 8
   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 8

   # out/SAMPLE.airr.tsv        one AIRR Rearrangement row per mapped read
   # out/SAMPLE.clones.tsv      the clonotype table
   # out/SAMPLE.arda.json       run report: counts, wall, peak RSS, versions
   # out/SAMPLE.stats.tsv       run QC, long format

``arda cells`` starts one step later — see :ref:`single cell <usage-singlecell>` below.

.. _usage-speed-levers:

.. warning::

   **The two speed configurations do not compose, and neither half is optional.** This is why
   the mode name picks the preset and ``--exact`` turns every lever off.

   * ``--two-pass`` **alone is a loss**: 0.762× on bulk and 0.87× on an IGH amplicon. It pays
     only together with ``--fast-segments``.
   * ``--fast-segments`` and ``--v-only-on-segment`` are ignored without ``--two-pass``.
   * ``--prefilter`` does **not** compose with ``--fast-segments``. Measured, the combination is
     slower than either lever on its own.

   The predictor is not the library type but ``fast_fraction`` in the ``--report`` JSON: the
   fraction of reads that hit **both** a V and a J segment. A primer-anchored TCR amplicon sits
   near .85; an IGH RepSeq library sits near .50; a bulk library sits near .05. Run 100 k reads,
   read ``fast_fraction``, then choose.

``arda map`` still exposes all five tuning flags (``--two-pass``, ``--fast-segments``,
``--v-only-on-segment``, ``--prefilter``, ``--indel-rescue``) individually for A/B work; every
one is off by default.

Bulk RNA-seq
------------

.. code-block:: bash

   # one shot -- the mode carries the configuration
   arda rnaseq --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 8

   # or stage by stage
   arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --report map.json \
       --threads 8 --prefilter
   arda assemble -i mapped.airr.tsv -o assembled.airr.tsv --report assemble.json
   arda correct  -i mapped.airr.tsv --extra-airr assembled.airr.tsv \
       -o clones.tsv --report correct.json

``--prefilter`` drops reads that share no exact 16-mer with the reference before ``createdb``
runs, so the FASTA write and the DB build are skipped along with the search. It costs ~0.5 % of
real reads — concentrated in J→C and hypermutated IGH — which is why it is off by default.

What the mode costs against MiXCR and TRUST4, end to end on a 660,000-pair bulk library, is under
`How it compares`_ below and in full in :doc:`benchmarks`.

Amplicon / RepSeq
-----------------

.. code-block:: bash

   # one shot -- the mode carries the configuration
   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 8

   # or stage by stage
   arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --report map.json \
       --threads 8 --two-pass --fast-segments --v-only-on-segment
   arda assemble -i mapped.airr.tsv -o assembled.airr.tsv --report assemble.json
   arda correct  -i mapped.airr.tsv --extra-airr assembled.airr.tsv \
       -o clones.tsv --report correct.json

.. important::

   ``correct`` takes the Stage-3 output through ``--extra-airr``. Omit it and the contigs
   ``assemble`` just built are silently discarded — the clonotypes whose CDR3 no single read
   spans never reach the table. ``arda rnaseq`` / ``arda amplicon`` wire this up for you.

How it compares
~~~~~~~~~~~~~~~

Against MiXCR and TRUST4, end to end to a clonotype table, each tool on its best-fitting preset:

* **TRA amplicon, 100,000 reads.** MiXCR is 1.84× faster on wall clock — that is its regime. arda
  reaches the same point on 1.41× less CPU and 3.16× less RSS, and returns the most clonotypes
  (19,841) over the most reads (43,503). TRUST4 is 5.1× slower than arda here.
* **Bulk RNA-seq, 660,000 pairs.** arda returns +27.8 % clonotypes and +97.9 % reads assigned
  against MiXCR, at 1.38× MiXCR's wall clock and 2.76× less RSS. TRUST4 is 1.61× faster on wall at
  4.0× less CPU, reaching 87.7 % of arda's clonotypes.
* **IGH RepSeq, 100,000 pairs, 32 threads.** The amplicon configuration is 4.15× and 4.71× faster
  than arda's own one-pass default on two libraries, at roughly 2.7× less memory.

Full tables, with the method, repetition spread and what each comparison is and is not entitled to
claim: :doc:`benchmarks`.

Accuracy in the amplicon configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Against an IgBLAST truth on the same 100,000-read TRA amplicon, scored per read from one truth
file in one job, arda and MiXCR at its best amplicon preset:

* **arda recalls more V genes over the whole library** — ``v_gene`` recall .9867 against .9660 —
  because MiXCR emits nothing at all for 1,530 of the 48,033 truth reads and arda for 3.
* **On the reads it emits, MiXCR is the more accurate caller**: ``v_gene`` recall .9977 and exact
  junction recall .9778 on the common subset, against arda's .9869 and .9533.
* **arda's V calls are the more precise of the two under either denominator** — .9996 against
  .9977 — because it declines rather than guessing.

Two habits make this kind of comparison meaningful, and both are worth adopting in your own.

**Print the coverage before the rates.** A per-tool inner join gives each tool its own
denominator, so a truth read the tool emitted no row for vanishes instead of counting as a miss.

**Score alleles as tie lists, not as exact strings.** IgBLAST and arda both report an ambiguous
allele call as a comma-joined set (``TRAV8-4*01,TRAV8-4*04,TRAV8-4*05``); an exact-string
comparison marks a correct-but-ambiguous call wrong. Across 25 datasets the median ``v_allele``
score is .8328 scored exactly against .9763 scored as tie-list membership — the same output, a
14-point difference in the scorer.

Both denominators, every metric and the allele-comparator caveat: :doc:`benchmarks`.

Junction correctness is discussed in :ref:`what a junction disagreement means`, which also says
which kinds of disagreement are *not* errors.

.. _ig-coverage-floor:

What an IG V call means when the read is short
----------------------------------------------

A ``v_call`` is only as good as the V germline the read actually covers, and on IG that is the
single largest thing separating a usable call from an unusable one. This is a property of the
library, not of the caller: no tool can name a gene from bases that are not in the read.

Measured against an ``arda igblast`` truth on an IGH 5'RACE library of 98,639 truth reads,
``v_gene`` recall rises monotonically with the V germline span the read covers, from **.1170 below
60 nt to .9896 at 200 nt or more** — and .9896 is the TRA amplicon's .9867, so there is no
IG-specific accuracy deficit at that coverage. The 5.9 % of reads carrying under 60 nt produce
about 56 % of every V miss. The full stratification, on 5'RACE and on two bulk libraries across
three loci, is in :doc:`benchmarks`.

Three things follow for how you read a call.

.. important::

   **Position beats length.** A short **5'** alignment identifies the gene; a short **3'** one does
   not. IGHV genes diverge in FR1/CDR1/CDR2 and are conserved through FR3 near Cys104, so 55 nt
   taken from the 3' end carries almost none of what separates one gene from another. The same
   ``< 60 nt`` bin scores **.1170** on a 5'RACE library, whose reads run C → J → V, against
   **.8472** on a bulk RNA-seq library of the same locus.

   The protocol name is not the risk factor either: the same 5'RACE protocol scores .1170 on a short
   read and .9645 at 251 nt.

.. note::

   **Somatic hypermutation does not order this.** Stratified by IgBLAST's own ``v_identity`` on the
   same 98,639 IGH reads the deficit is *non-monotonic* — essentially unmutated reads score .9425,
   the 97–99 % bin **.7518** (the worst of six), and the 92–95 % bin **.9697** (the best). Mutation
   load is not what makes an IG V call hard here; read geometry is.

**arda ships no threshold on this, and the span is yours to filter on.** Dropping a V call whose
alignment starts at germline position ≥ 150 and spans < 60 nt was measured: on IGH 5'RACE it is a
7.6 : 1 win (``v_gene`` precision .92585 → .97693 for −0.69 points of recall), and on bulk it is a
clear loss — −3.24 points of IGH recall, −6.22 IGK, −2.39 IGL, for essentially no precision. One
constant would be wrong more often than right. ``v_sequence_start`` / ``v_sequence_end`` and
``v_germline_start`` / ``v_germline_end`` are written on every row for you to gate on yourself.

.. warning::

   **Map the pair, not one mate — especially on a C-anchored amplicon.** On a multiplex V-primer
   IGH amplicon (251 nt paired), mapping **R1 alone** returns ``j_call`` and ``c_call`` on 98 % of
   reads and ``v_call`` on **1.9 %**, with ``junction`` on **40 of 49,036 rows** — while exiting 0
   and reporting 98.07 % of reads mapped. Nothing in the report moves.

   The cause is reference geometry rather than the library. A V·J scaffold carries no constant
   region, so a read that runs C → J → V and reaches ~80 nt into the constant region scores 244 bits
   on a :term:`J+C scaffold` against 241 on its own V·J scaffold — and the J+C target has no V to
   report. Trimming the constant region off the same reads moves ``v_gene`` recall from .0265 to
   .8775 and junctions from 40 to 40,005.

   ``arda amplicon --r1 --r2`` is the answer and needs no flag: it annotates each mate and Stage 3
   bridges them. On those same 50,000 pairs that is **24,655 contigs, 100 % of which carry**
   ``v_call``, ``j_call``, ``c_call``, ``junction`` **and** ``junction_aa``. This is why
   ``--no-assemble`` is not a preset on any mode.

.. _resolve-ties-loci:

Widening the V call, and the loci where that helps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``arda resolve-ties`` widens ``v_call`` to every germline the read's alignment cannot rule out. It
is off by default and it is **not** uniformly good: measured against an ``arda igblast`` truth on
eleven arms across five geometries, scored on **exact** ``v_gene``-set agreement rather than
intersection —

.. list-table:: exact ``v_gene``-set agreement with IgBLAST, before and after ``resolve-ties``
   :header-rows: 1
   :widths: 12 30 18 18 22

   * - locus
     - arms
     - before
     - after
     - change
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

.. code-block:: sh

   arda resolve-ties -i mapped.airr.tsv -o widened.airr.tsv --loci IGK,IGL

``--loci`` widens only the named loci and copies every other row through untouched, so one command
on a mixed library takes the win on IGK and IGL and leaves IGH byte-identical to the input.

The reason is mechanical. The deciding quantity is the **ambiguity deficit** — IgBLAST's genes per
read minus arda's, *before* widening. arda's IGK call is **0.42 genes/read too narrow** (1.18
against 1.60) and closing that gap is the whole 37-point win; arda's IGH call is already as wide as
IgBLAST's (deficit 0.00–0.07 on every arm, 1.589 against 1.64 on 5'RACE), so there is nothing to
recover and each IGH arm only pays the overshoot. **You can check this on your own library with no
truth file**: compare the mean number of genes per ``v_call`` before and after, and if it moves far
more than a few hundredths the rule is overshooting on that locus.

.. warning::

   Two things this trades away. **Fewer reads name a single gene** — on IGK the share falls
   ``.8195 → .3883``, which is exactly what an IGK read supports, and why this stays an opt-in
   command. And **intersection recall rises on all eleven arms, IGH included**, so a
   scorer that only asks "does arda's list contain the truth gene" reports the 10-point IGH
   regression as a 5-point IGH win.

   Clonotype cost, since ``v_call`` is part of the clonotype key: **+2 of 362 on IGK, +3 of 23,559
   on TRB, and zero on IGL**.

Allele-level separation needs considerably more than gene-level identification does — IGHV needs a
median **230 nt** from the 3' end to separate a gene's alleles, against 175 for TRAV and 150 for
TRBV. See :doc:`genotype` for that table and what it means for ``arda genotype``.

.. _usage-singlecell:

Single cell
-----------

``arda cells`` takes **one UMI consensus per molecule** with the cell barcode in the record name —
what ``migec assemble`` writes — and assembles each cell's contigs with **no germline reference**,
then annotates them, pairs the chains and flags multiplets:

.. code-block:: bash

   arda cells asm/PBMC.consensus.fq.gz -p out/PBMC --cells ref/PBMC.cells.tsv --plot svg

arda does no demultiplexing, no barcode correction, no UMI collapse and no cell calling; the
upstream tool does all four and leaves the barcode in the record name. The dialect is sniffed
(``migec``, ``cellranger``, ``prefix``) or named with ``--cell-from`` / ``--cell-regex``. Measured
against Cell Ranger on ``sc5p_v2_hs_PBMC_1k`` VDJ-T, **98.2 %** of its 943 CDR3s appear verbatim
inside one of arda's contigs. Full write-up, including the QC notebook: :doc:`singlecell`.

The three stages
----------------

The modes run these for you; run them separately when you want to shard Stage 1 across a cluster.

**map** streams paired FASTQ and writes only the reads that map to a receptor scaffold. The
reference includes ``J + C`` constant-region scaffolds, so a read spanning the J→C splice still
maps and carries a ``c_call`` (CH1 exon) and a ``c_class`` isotype (``IGHG``/``IGHM``/``IGHA`` …
— the class, never the subclass). With ``--reconstruct``, overlapping mates are merged into one
fragment; a mismatch in the overlap is resolved base-by-base in favour of the higher-Phred call.

**assemble** reconstructs clonotypes with a CDR3 too long for any single 100–150 bp read to span
(V(DD)J ultralong, ~20–40 aa): it anchors on Stage-1's per-read ``cdr3_start`` and grows contigs
by greedy overlap-extension, then folds the recovered reads back into ``correct``.

**correct** aggregates reads into clonotypes keyed by ``(locus, v_call, j_call, junction)`` and
collapses sequencing-error CDR3 variants. Abundance is the AIRR ``duplicate_count`` — **every
read that encompasses the junction** (spanning or partial, assigned by alignment), the true
expression estimate — with ``consensus_count`` giving the distinct-fragment count. Error
correction uses a per-base sequencing-error model whose threshold scales with junction length
(``error_rate * junction_len`` per substitution; ~1/20 at a 45 nt junction), tolerant of
somatic-hypermutation indels, and tested only over reads that observe the discriminating
position. Each clonotype's D is mapped once into its error-corrected junction
(``d_call``/``d2_call``/``d_support``), not voted over reads — D is a function of the junction.
See :doc:`error_correction` for the abundance model, what ``--error-rate``/``--max-subs``
actually do, and where the method reaches its limit.

.. tip::

   A sample delivered as several FASTQ pairs (lanes, chunks) is handled in one call: ``--r1`` /
   ``--r2`` / ``--id`` are repeatable, or pass an nf-core-shaped ``--samples`` sheet. The result
   is byte-identical to the same reads concatenated. See :doc:`samples`.

What a run writes
-----------------

The output is a **spec-valid AIRR Rearrangement** TSV (it passes ``airr.schema`` validation) with
1-based, closed region coordinates (``fwr1_start``/``fwr1_end`` … ``cdr3_start``/``cdr3_end``),
region nucleotide and amino-acid sequences, ``v_call``/``d_call``/``d2_call``/``j_call``, the
constant-region ``c_call``/``c_class`` (isotype), per-segment CIGARs
(``v_cigar``/``d_cigar``/``j_cigar``/``c_cigar``) and the matching V/J/D germline coordinates,
``sequence_alignment`` / ``germline_alignment``, ``v_identity``, the per-segment SHM lists
``v_mutations`` / ``j_mutations``, ``stop_codon``, ``vj_in_frame``, ``junction``, and
``productive``.

On a score tie ``d_call``/``d2_call`` are comma-separated allele ambiguity lists, and ``d2_call``
is the second (3′) segment of a D-D fusion — called in every D locus (IGH, TRB, TRD). The
``sequence`` field holds the read **as submitted**; ``rev_comp`` = ``T`` signals that the other
output fields describe its reverse complement (per the AIRR spec).

The D call is accepted on a Karlin–Altschul E-value (``d_support``, shipped so a consumer can
re-threshold), and is constrained by germline geometry: TRBD2 lies 3′ of the whole TRBJ1 cluster,
so no TRBJ1 rearrangement is ever assigned TRBD2. D mapping also runs on ``--seqtype aa`` input,
against each D germline's three translated frames — informative for IGH (a D call on ~36 % of real
records, agreeing with the nucleotide call on 98 % of them) and mostly silent for the TR loci,
whose D is too short to survive trimming into protein. For aa input the ``d_germline_*`` columns
and ``d_cigar`` are left empty on purpose: the alignment offsets index a reading frame, not the D
germline.

Optional columns
----------------

All four are off by default, so the shipped output never moves under you.

Somatic hypermutation
~~~~~~~~~~~~~~~~~~~~~

``v_mutations`` and ``j_mutations`` are the read's substitutions against its called germline —
``G45A,C112T``: germline base, 1-based position **in that segment's own allele**, read base. That is
the coordinate frame a lineage or selection-pressure tool needs, so two reads of one clone are
directly comparable and the germline is the root. A read with none is empty; the counterpart
``v_identity`` is the same information as a fraction.

.. note::

   **The list covers the V and J germline-aligned regions only.** A mismatch inside the junction is
   not attributable to a germline: V(D)J recombination trims the segment ends and inserts
   non-templated N/P bases, so the V-end / NDN / J-start partition is frequently not identifiable
   from the sequence. arda aligns to a ``V + N-pad + J [+ C]`` scaffold, and the pad is not a
   segment — an NDN position has no germline coordinate to be recorded under.

   Diffing ``sequence_alignment`` against ``germline_alignment`` by hand does **not** reproduce this
   list: on a real bulk IG library, 20.1 % of the mismatches that diff finds lie in the pad or the
   constant region.

Substitutions only; an indel is in the CIGAR as ``I``/``D``. Germline coordinates after an indel are
still correct. Positions are on the coding strand, so for ``rev_comp = T`` a read-side lookup (a
Phred quality, say) must be made against the reverse complement of ``sequence``.

Accuracy against IgBLAST, and what a lineage-tree builder needs on top of this: :doc:`shm`.

D segments
~~~~~~~~~~

``--d-max-evalue`` moves the gate that accepts a D call (and a tandem second D) — the shipped
operating point is 0.2 for nt and 0.05 for aa, and ``0.01`` is the band where D agrees .9985 with
IgBLAST at gene level on a TRB amplicon, at roughly a third of the call rate. Germline geometry is
applied before the statistics: ``/OR`` orphons cannot rearrange and are excluded, TRBD2 can never
join a TRBJ1, and a tandem D-D must run in genomic order. The bands, the constraints and how to
consume the D-D markup: :doc:`d_segments`.

Quality over the junction, and at each mutation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``arda map --junction-quality`` adds a ``junction_quality`` column — the read's Phred+33 string
over exactly the bases of ``junction``, same orientation. Stage 1 is the only place the FASTQ
quality is still in hand, and it is what ``correct --min-junction-q`` gates on: see
:ref:`quality gate`.

``arda map --mutation-quality`` adds ``v_mutation_quality`` / ``j_mutation_quality`` — the Phred of
the read base behind each entry of ``v_mutations`` / ``j_mutations``, comma-joined, one-for-one and
in the same order. A novel allele, somatic hypermutation and a base miscall are the same string in
the mutation list; the recurrence separates the first from the second and the Phred separates both
from the third. ``arda stats`` reads it to score its ``allele_candidate`` shortlist.

.. warning::

   **The two quality columns use different encodings.** ``junction_quality`` is raw Phred+33
   characters, so it lines up byte for byte with ``junction``; ``v_mutation_quality`` is comma-joined
   integers. Reading one as the other gives plausible numbers off by 33. Both are refused with
   ``--reconstruct``. See :doc:`qc`.

Finishing a truncated junction from the germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A read can reach Cys104, run into the J, and stop before [FW]118 — a junction arda declines
because its 3' boundary was never observed. ``--complete-junctions N`` finishes it from the called
J's germline, taking at most ``N`` nt.

This is sound for one reason, and only in one direction: the J's 5' chew-back and the N/P additions
all lie **upstream** of the read's last aligned J base, so everything from there to [FW]118 is
germline-**templated**. The V side has no counterpart — a read short at the 5' end is missing bases
the V germline does not template either, which is why ``v_anchor_prefix`` refuses rather than
extrapolates.

**The added bases are imputed, not observed.** Every completed row carries the count in
``junction_completed_nt``, so a consumer can filter or weight on that column rather than trusting
the junction blindly. An empty value means the junction is entirely observed — which is what every
junction is unless the flag is passed. Off by default (``0``).

Measured (arda 2.17.0, human, ``--complete-junctions 40``):

.. list-table::
   :header-rows: 1

   * - library
     - junctions
     - completed
     - median nt imputed
   * - bulk RNA-seq, 1 M pairs (SRR5233639)
     - 3,856 → 4,103 (**+6.4 %**)
     - 247
     - 13
   * - TRA amplicon, 100 k reads
     - 44,497 → 44,527 (+0.07 %)
     - 30
     - 14

No junction that was already observed moves in either arm. 246 of the 247 bulk completions close on
[FW]118; the one that does not is ``TRBJ2-7*02``, whose anchor codon really is not [FW] — the same
allele biology as ``TRAJ35*01``'s Cys anchor, and read from ``anchor_nt`` rather than from a motif.

.. warning::

   Two caveats, both measured. On IG the imputed span can hide the somatic hypermutation the read
   would have shown, biasing a completed junction's 3' end toward germline — and IG is where the
   yield is, with 175 of the 247 bulk completions on IGH.

   And a read whose alignment stops more than a partial codon short of its own 3' end is **refused**:
   on the TRA amplicon, 236 of 266 candidates run from the V straight into ``TRAC`` with no J at all
   while the aligner still names a J off a few coincidental bases. Completing those would have
   manufactured one junction per chimera.

Sequences that were never reads
-------------------------------

``arda annotate`` takes FASTA or FASTQ of assembled sequences, nucleotide or protein:

.. code-block:: bash

   arda info                                             # resolved paths and tool availability
   arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
   arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa
   arda stats -i out.airr.tsv -o out.stats.tsv           # run QC, long-format TSV

``-v`` / ``-q`` / ``--log-file`` are global and go **before** the subcommand
(``arda -v --log-file run.log map ...``); progress goes to stderr and results to stdout.

From Python, the same annotation path:

.. code-block:: python

   import arda

   records = arda.annotate_sequences(
       ["GACGTGCAG...", ("clone7", "CAGGTG...")],
       seqtype="nt",
       organism="human",
   )

Each record is a dict keyed by the AIRR fields above.

``arda igblast -i reads.fastq -o truth.airr.tsv`` runs IgBLAST across all loci as a gold-standard
reference for benchmarking. ``--receptor ig`` or ``tr`` skips a whole pass over every read on a
library whose receptor type is known; the default ``both`` runs both passes and keeps whichever
scores higher, which is what a truth file for an unknown library wants.

Bare germline segments
~~~~~~~~~~~~~~~~~~~~~~

There is no coverage filter, so a **V-only** or **J-only** query maps to its scaffold and only the
regions inside the query's coverage come back: a bare V yields ``fwr1`` through ``fwr3``, a bare J
yields ``fwr4``. That makes the annotator usable as a way to pull per-allele FR/CDR subsequences out
of the germline:

.. code-block:: python

   from arda.annotate.mapper import annotate_records

   recs = annotate_records(
       [("TRBV9*01", v_germline_nt), ("TRBJ2-7*01", j_germline_nt)],
       organism="human", seqtype="nt", strand="forward", map_d=False,
   )

mirpy uses exactly this to bake per-allele FR/CDR subsequences into its gene library; see
``tests/synthetic/test_germline_segments.py``. :doc:`arda export-ref <reference_export>` is the CLI
equivalent, and is usually the better route when you want the whole reference rather than a handful
of alleles.

Junction markup and repair
~~~~~~~~~~~~~~~~~~~~~~~~~~

``arda markup`` works on records that have **no read behind them** — a CDR3 amino acid, a V
call and a J call, as in a VDJdb row. It reports which residues each germline templates, where
the submitted junction disagrees with them, and how far the disagreement extends:

.. code-block:: bash

   arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report -
   arda markup -i records.tsv -o marked.tsv --organism human --d-posterior

Coordinates are **junction space** throughout (Cys104 … Phe/Trp118, both anchors included) —
the convention VDJdb's ``cdr3`` column uses, which is *not* arda's ``cdr3`` field. Output adds
``v_end``/``j_start`` in residues and ``v_end_nt``/``j_start_nt`` in nucleotides, a per-error
list (substitution / insertion / deletion, with position and extent), a VDJdb-compatible
``cdr3fix`` JSON blob, and a repaired ``cdr3_repaired``. Repair is
deliberately conservative: only anchor-adjacent edits are *applied* (``--max-replace``), while
errors deeper in the junction are reported and left alone — on 102,990 VDJdb records this
reproduces VDJdb's own repair on 96.4 % of the records it marks as needing one, and rewrites
nothing it should not.

**Read the boundary in nucleotides.** A germline run ends wherever the exonuclease stopped,
which is not a codon boundary, so the residue counts round it — ``v_end`` is exact on 71.8 % of
junctions with external nucleotide truth against ``v_end_nt``'s 92.9 % under the same VDJdb
residue convention (:func:`arda.cdr3fix.boundary_nt` carries the whole measurement). The residue
counts keep their meaning and their callers; they are simply the coarser answer.

``--d-posterior`` adds a D-gene call inferred from the junction *length* — the nucleotide length
pins ``insVD + |D surviving| + insDJ``, so the D can be placed to a median 1–3 nt even when the
protein shows nothing of it. Available for human IGH/TRB/TRD and mouse TRB, the pairs with a
published generative model; everything else returns nothing rather than guessing.

Running at scale
----------------

MMseqs2 runs multi-threaded (``--threads``); inputs may be FASTA or FASTQ, plain or gzipped. To
spread one sample across a cluster, shard **Stage 1 only** — ``arda cluster submit`` does this and
the result is byte-identical to the single-node run. The SLURM scripts, the two cluster adapters
and why they are not interchangeable: :doc:`cluster`. Nextflow and Snakemake:
:doc:`pipeline_integration`.

Run reports
-----------

``--report`` (per stage) and the modes (merged, ``<prefix>.arda.json``) record what
happened and what it cost. Three resource fields appear on every stage:

``wall_seconds``
   Elapsed time in that stage.

``peak_rss_mb``
   The **whole-process** (including the child ``mmseqs``) high-water mark **as of the end of
   that stage** — so it is monotone across stages. ``resource.getrusage`` reports high-water
   marks only and offers no per-stage reset, so when all three stages share one process a stage
   cannot be charged its own peak in isolation. This is deliberately the number to size a SLURM
   ``--mem`` or Nextflow ``memory`` directive from: it is what the process actually required by
   that point.

``rss_gain_mb``
   How much *that* stage raised the mark. For an unambiguous per-stage figure, run the stage in
   its own process (``arda map`` / ``assemble`` / ``correct`` separately); then
   ``peak_rss_mb`` is that stage alone.

.. note::

   Budget for **Stage 3**, not Stage 1. Mapping is flat at ~300-650 MB regardless of read
   depth, but ``assemble``/``correct`` hold the clone set: on a B-cell-rich tumour
   (28,444 clonotypes from 105 M reads) Stage-3 ``correct`` peaked at **2,071.7 MB**, while a
   colder sample with *more* reads (139 M) peaked at **549 MB**. Peak RSS tracks repertoire
   richness, not read depth. Budget ~4 GB.

The report also carries ``arda_version``, ``mmseqs_version`` and a ``reference`` fingerprint
(path, size, mtime). If two runs disagree, compare those first — a different aligner build or a
differently-fetched reference is the usual cause, and neither is visible in the output.

``arda stats`` turns a finished run into one long-format QC TSV — reads and clonotypes per chain,
functional / non-functional / stop-codon / truncated-junction counts, junction length and quality,
SHM rate, chimeras, V/J gene coverage and a candidate-allele shortlist. The mode commands write it
automatically as ``<prefix>.stats.tsv``. Full description, plus the verbosity and ``--log-file``
options: :doc:`qc`.

Supported organisms
-------------------

* **human, mouse** — full IG and TR loci.
* **rat, rabbit, rhesus_monkey** — IG only (IgBLAST ships no TR internal
  annotation for these organisms).

.. note::

   Every command on this page has a runnable counterpart in ``examples/``, built from real data
   committed to the repo and regenerated by ``python examples/regenerate.py``: one mRNA per locus,
   the two human reads that carry a tandem D-D, six VDJdb records covering every junction-repair
   outcome, and a 1,035-read FASTQ that runs the whole bulk RNA-seq pipeline in about six seconds.
