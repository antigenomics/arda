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

Measured on 100 k bulk RNA-seq reads:

.. list-table::
   :header-rows: 1
   :widths: 28 18 18 18

   * - tool
     - wall (s)
     - CPU (s)
     - peak RSS (MB)
   * - arda
     - 2.51
     - 5.40
     - 234
   * - MiXCR 4.7.0
     - 4.54
     - 31.80
     - 3,022
   * - TRUST4
     - **1.91**
     - **4.36**
     - **192**

.. note::

   TRUST4's row is **not the same stage of work**. It measures candidate *extraction* — which
   reads look receptor-derived — while arda's and MiXCR's rows measure a per-read AIRR
   Rearrangement record with gene calls and a junction. TRUST4 does that annotation later, on
   ~1,000 assembled contigs rather than per read. Quote the row only with that caveat attached.

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

⛔ **A stage comparison is not a benchmark.** The legs below run **end to end and emit
clonotypes**; each tool gets its best-fitting preset and nothing else. One job, six legs
alternating, three reps, medians. 8 threads; arda 2.27.0, MiXCR 4.7.0, TRUST4. TRUST4's rows
count only its **complete** CDR3s, so the last two columns mean the same thing in every row.

TRA amplicon, 100,000 reads:

.. list-table::
   :header-rows: 1
   :widths: 26 13 13 15 16 17

   * - pipeline
     - wall (s)
     - CPU (s)
     - peak RSS (MB)
     - clonotypes
     - reads in clonotypes
   * - MiXCR ``generic-amplicon``
     - **7.82**
     - 42.91
     - 3,052
     - 19,697
     - 42,712
   * - **arda** ``amplicon``
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

MiXCR is **1.84× faster on wall** in its own regime; arda gets there on **1.41× less CPU** and
**3.16× less RSS**, and returns the most clonotypes over the most reads. TRUST4 is 5.1× slower
than arda here.

Bulk RNA-seq, 660,000 pairs (``SRR5233639``):

.. list-table::
   :header-rows: 1
   :widths: 26 13 13 15 16 17

   * - pipeline
     - wall (s)
     - CPU (s)
     - peak RSS (MB)
     - clonotypes
     - reads in clonotypes
   * - TRUST4
     - **13.68**
     - **54.53**
     - **460**
     - 1,941
     - 6,247
   * - **arda** ``rnaseq``
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

On the regime arda exists for it returns **+27.8 % clonotypes and +97.9 % reads assigned** against
MiXCR (+14.0 % / +35.8 % against TRUST4), at 1.38× MiXCR's wall and 2.76× less RSS. ⚠ TRUST4 is
genuinely 1.61× faster on wall at 4.0× less CPU on this arm, reaching 87.7 % of arda's clonotypes
and 73.6 % of its assigned reads.

On real IGH RepSeq at 32 threads (aldan3, 100 k pairs), against arda's own shipped one-pass
default:

.. list-table::
   :header-rows: 1
   :widths: 26 18 18 18 20

   * - dataset
     - default (s)
     - amplicon cfg (s)
     - default RSS (MB)
     - amplicon cfg RSS (MB)
   * - IGH_repertoire
     - 316.44
     - **76.25**
     - 4,018
     - **1,479**
   * - IGH_naive
     - 305.32
     - **64.86**
     - 3,736
     - **1,363**

4.15× and 4.71× respectively, at roughly a third of the memory.

Accuracy in the amplicon configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Against an IgBLAST truth on the same 100 k-read TRA amplicon, arda 2.27.0 and MiXCR 4.7.0 at its
best amplicon preset, scored per read from one truth file in one job.

.. important::

   **Print the coverage before the rates.** A per-tool inner join gives each tool its own
   denominator — a truth read the tool emitted no row for vanishes instead of counting as a miss.
   Of the 48,033 truth reads at ``v_score >= 70``, **arda emits a row for 48,030 (99.99 %) and
   MiXCR for 46,503 (96.81 %)**, so both denominators are reported below.

.. list-table::
   :header-rows: 1
   :widths: 40 16 16 16

   * - metric
     - arda, all truth
     - MiXCR, all truth
     - arda / MiXCR, common
   * - ``v_gene`` recall
     - **.9867**
     - .9660
     - .9869 / **.9977**
   * - ``v_gene`` precision
     - **.9996**
     - .9977
     - **.9997** / .9977
   * - ``j_gene`` recall
     - **.9892**
     - **.9892**
     - .9959 / **.9996**
   * - ``j_gene`` precision
     - .9953
     - **.9996**
     - .9979 / **.9996**
   * - ``junction`` recall (nt, exact)
     - .9473
     - **.9708**
     - .9533 / **.9778**

The two views say different and equally true things. **On the reads it emits MiXCR is the more
accurate caller**; **over the whole library arda recalls more V genes** (.9867 vs .9660), because
MiXCR emits nothing at all for 1,530 truth reads against arda's 3. arda's V calls are the more
precise of the two under either denominator — it declines rather than guessing. Of the 46,787
truth junctions arda emits one for **94.81 % at .99919 precision among emitted**, MiXCR for
97.09 % at .99989; the 5.19 % arda declines are reads that *have* an anchor pair and lost the
projection.

MiXCR emits ``*00`` for every allele, i.e. it makes no allele call, so there is no ``v_allele``
comparator. arda's is **.9868 resolved** and .9461 by exact string.

.. note::

   These five arda figures are **unchanged from 2.11.1**, fifteen releases back — ``v_gene``
   recall .9867, precision .9996, ``j_gene`` recall .9892, precision .9953 and junction precision
   among emitted .99919 all reproduce to every published digit.

.. note::

   **Score alleles as tie lists, not as exact strings.** IgBLAST and arda both report ambiguous
   allele calls as comma-joined sets (``TRAV8-4*01,TRAV8-4*04,TRAV8-4*05``); an exact-string
   comparison marks a correct-but-ambiguous call wrong. Across 25 cluster datasets the median
   ``v_allele`` score is **.8328** scored exactly against **.9763** scored as tie-list
   membership — the same output, a 14-point difference in the scorer.

Junction correctness is discussed in :ref:`what a junction disagreement means`, which also says
which kinds of disagreement are *not* errors.

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

Never: The V and J germline-aligned regions only, **by construction**. A mismatch inside the junction is
not attributable to a germline: V(D)J recombination trims the segment ends and inserts non-templated
N/P bases, so the V-end / NDN / J-start partition frequently is not identifiable from the sequence.
arda aligns to a ``V + N-pad + J [+ C]`` scaffold, and the pad is not a segment — an NDN position has
no germline coordinate to be recorded under. Diffing ``sequence_alignment`` against
``germline_alignment`` by hand does **not** give you this: on a real bulk IG library 20.1 % of the
mismatches that diff finds lie in the pad or the constant region.

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

Never: The two quality columns use **different encodings**: ``junction_quality`` is raw Phred+33
characters (it lines up byte-for-byte with ``junction``), and ``v_mutation_quality`` is comma-joined
integers. Both are refused with ``--reconstruct``. See :doc:`qc`.

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

**Never: The added bases are imputed, not observed.** Every completed row carries the count in
``junction_completed_nt``, so a consumer filters or weights on that column instead of trusting the
junction; an empty value means the junction is entirely observed, which is what every junction is
unless the flag is passed. Off by default (``0``).

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

⚠ **Two caveats, both measured.** On IG the imputed span can hide the SHM the read would have
shown, biasing a completed junction's 3' end toward germline — and IG is where the yield is (175 of
the 247 bulk completions are IGH). And a read whose alignment stops more than a partial codon short
of its own 3' end is **refused**: on the TRA amplicon 236 of 266 candidates run from the V straight
into ``TRAC`` with no J at all, while the aligner still names a J off a few coincidental bases.
Completing those would have manufactured one junction per chimera.

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
reference for benchmarking.

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
``v_end``/``j_start``, a per-error list (substitution / insertion / deletion, with position and
extent), a VDJdb-compatible ``cdr3fix`` JSON blob, and a repaired ``cdr3_repaired``. Repair is
deliberately conservative: only anchor-adjacent edits are *applied* (``--max-replace``), while
errors deeper in the junction are reported and left alone — on 102,990 VDJdb records this
reproduces VDJdb's own repair on 96.4 % of the records it marks as needing one, and rewrites
nothing it should not.

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
