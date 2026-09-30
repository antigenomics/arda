How arda works
==============

Aligning a receptor read is expensive mostly because deciding *where the regions are* is
expensive. arda pays that cost once, offline, and then reduces every runtime query to a search
plus a coordinate projection.

.. contents::
   :local:
   :depth: 2

The idea
--------

IgBLAST is the reference implementation for receptor annotation, and its region calls are what
arda reproduces. What makes it awkward at scale is that it re-derives the V–D–J structure for
every query, per invocation.

arda inverts that. Before any data arrives, it enumerates every germline **combination** a
rearrangement could use and annotates each one with IgBLAST, once. The result is a reference in
which the regions are already known:

.. code-block:: text

   offline, once                          at runtime, per query
   ─────────────                          ─────────────────────
   IMGT germlines                         query sequence
        │                                      │
        ├─ enumerate V·J scaffolds             ├─ MMseqs2 search against the scaffolds
        ├─ annotate with igblastn              ├─ best hit → its known region coordinates
        └─ store FR1–4 / CDR1–3 coords         ├─ project those coordinates onto the query (C++)
                                               ├─ local-align the D germlines in the junction
                                               └─ AIRR Rearrangement row

Annotating a query therefore costs one homology search and one linear walk of a CIGAR string.
IgBLAST is never invoked at runtime — it is needed only to rebuild the reference, or to produce a
gold-standard comparison with ``arda igblast``.

The reference
-------------

Scaffolds
~~~~~~~~~

A **scaffold** is one V allele joined in frame to one J allele, with a short frame-neutral N
spacer where D would sit. D is deliberately not enumerated: it only affects the interior of the
junction, which is query-specific, so enumerating it would multiply the reference without adding
information.

For human, arda ships:

.. list-table::
   :header-rows: 1
   :widths: 14 20 20 20 26

   * - Locus
     - V·J scaffolds
     - J+C scaffolds
     - D germlines
     - Status
   * - TRA
     - 6,767
     - 69
     - 0
     - ok
   * - TRB
     - 2,080
     - 32
     - 3
     - ok
   * - TRG
     - 105
     - 12
     - 0
     - ok
   * - TRD
     - 73
     - 5
     - 3
     - ok
   * - IGH
     - 4,225
     - 154
     - 48
     - ok
   * - IGK
     - 819
     - 10
     - 0
     - ok
   * - IGL
     - 1,000
     - 63
     - 0
     - ok
   * - **Total**
     - **15,069**
     - **345**
     - **54**
     -

Only in-frame combinations are kept, and IMGT ORFs and pseudogenes are excluded, so arda's call
vocabulary is narrower than IgBLAST's germline database. What that costs, and the three checks
applied before a gene is added or removed, is recorded with the reference build.

J+C scaffolds
~~~~~~~~~~~~~

A read spanning the J→C splice ends in the constant region and has no V to anchor on. Against a
V·J reference it has nowhere to land. arda therefore also ships **J+C scaffolds** — the CH1 exon
spliced onto each J, with no V — which gives such a read a target and lets arda report its
``c_call`` and, for IGH, its isotype ``c_class``.

On bulk RNA-seq this class is 22 % of real receptor fragments. Adding these scaffolds took
agreement on it from 0.0606 to 0.9844 and overall recall from 0.78 to 0.986 (16 datasets,
5,273 real fragments).

Per-allele junction anchors
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``cdr3_anchors.tsv`` records, for each of 1,240 human alleles, the conserved Cys104 (V) or
Phe/Trp118 (J) position, the residues that allele templates into the junction, and a status —
``ok``, ``truncated`` or ``no_anchor``. An allele with no derivable anchor is flagged rather than
given a guessed one.

These anchors, not the scaffold projection, are what bound the junction at runtime. The reason is
geometric: a scaffold has a 9 nt N pad where a real read has a 20–40 nt N-D-N region, so
projecting the scaffold's own coordinates would collapse exactly the window the D lives in.

.. note::

   A conserved-motif test is not a substitute for the anchor table. ``TRAJ35*01``'s anchor codon
   decodes Cys (TGC) rather than Phe or Trp, and it is a functional IMGT ``F`` gene — a
   ``[FW]GXG`` motif check silently deletes it.

The runtime path
----------------

1. Search
~~~~~~~~~

MMseqs2 searches the query against the scaffold reference and returns the best hit with its score,
E-value and aligned spans. Nucleotide input is searched on **both strands** by default;
reverse-complement reads are re-oriented and flagged ``rev_comp=T``. All loci are in one reference,
so a mixed bulk RNA-seq file is annotated across every locus in a single search.

Ties are kept as ties. A read aligned over ``[germline_start, germline_end]`` is explained exactly
as well by any germline carrying that same stretch, so ``v_call`` is a comma-joined set whenever
the alignment cannot separate them.

2. Projection
~~~~~~~~~~~~~

The C++ ``transfer_regions`` walks the hit's CIGAR and maps each of the scaffold's known region
boundaries into query coordinates, handling insertions, deletions, truncation at either end,
mid-codon alignment starts and the reverse strand. This is the hot path and the reason the
extension exists.

3. The junction and the D call
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The V..J interior is bounded by the per-allele anchors. A junction is emitted **only** when the
read actually reaches its Cys104 anchor — a read that lands mid-V gets an empty junction rather
than a shortened one.

Inside that interval, a gapless C++ local alignment against the locus's D germlines produces
``d_call`` and, when a second D fits, ``d2_call``, plus the non-templated stretches ``np1``,
``np2``, ``np3``. D mapping also runs on amino-acid input, against each D germline's three
translated frames.

The D call is accepted on a Karlin–Altschul E-value, reported as ``d_support``, rather than on a
per-locus score floor. Before the statistics are consulted at all, germline geometry rules out
what cannot have happened:

* an ``/OR`` orphon sits outside the locus and cannot rearrange;
* TRBD2 lies 3' of the entire TRBJ1 cluster, so a TRBJ1 rearrangement can never carry it;
* a tandem D-D must run in genomic 5'→3' order, because deletional joining cannot produce any
  other.

MMseqs2 is not used for D. An 8–31 nt germline is below the length where a k-mer prefilter is
reliable, which is why this one alignment is done directly in C++. Details in :doc:`d_segments`.

4. Output
~~~~~~~~~

Regions are translated, productivity is evaluated over the annotated span, and the row is written
as AIRR. Out-of-frame junctions are reported with an N-bridge (``_``) so FR4 still reads.

From reads to clonotypes
------------------------

The three pipeline modes run three stages. The split matters because only the first is per read.

.. list-table::
   :header-rows: 1
   :widths: 16 14 70

   * - Stage
     - Scope
     - What it does
   * - ``map``
     - per read
     - Stream FASTQ, keep the reads that map to a receptor scaffold, annotate them, write AIRR.
       With ``--reconstruct``, overlapping mate pairs are merged into one fragment first, resolving
       overlap mismatches by the higher-Phred base. In paired mode the isotype of a CDR3-bearing
       read is recovered from its constant-region mate.
   * - ``assemble``
     - global
     - Reconstruct clonotypes whose CDR3 is too long for any single 100–150 bp read to span
       (V(DD)J ultralong, roughly 20–40 aa), by greedy overlap extension anchored on Stage 1's
       per-read ``cdr3_start``.
   * - ``correct``
     - global
     - Aggregate reads into clonotypes and collapse sequencing-error variants of a junction onto
       their parent.

Only Stage 1 shards. Error correction compares a clonotype against its neighbours by abundance, so
running it per shard would ask the question against a fraction of the evidence. ``arda cluster``
shards Stage 1 and runs Stages 2–3 once over the merged output, which is what makes a sharded run
byte-identical to a single-node one.

Single cell takes a different route: ``arda cells`` assembles each cell's contigs with **no
germline reference** — a cell's UMI-consensus molecules tile the transcript even though no single
molecule spans it — and annotates the contigs afterwards. See :doc:`singlecell`.

Determinism
-----------

A sharded run, a Nextflow run and a single-node run over the same input produce byte-identical
output. That is a requirement rather than a coincidence: polars' ``group_by`` is a multithreaded
hash aggregation, so Stage 2 sorts its groups on a total key and sorts the reads inside each group.
Without those sorts, three runs of the same input produced three different tables with the same
row count.

.. tip::

   When comparing two runs, compare a **call digest** over
   ``(locus, v_call, j_call, junction, duplicate_count)`` rather than a row count. A count-equal,
   digest-different pair of runs is invisible to any count.

The C++ extensions
------------------

Four small nanobind modules, kept separate so an experimental one can be opt-in while a shipped
one stays load-bearing.

.. list-table::
   :header-rows: 1
   :widths: 16 20 64

   * - Module
     - Availability
     - What it holds
   * - ``_markup``
     - always on
     - The hot path: ``transfer_regions`` and ``project_region`` (CIGAR walk and coordinate
       projection), AIRR row formatting, ``translate``, ``reverse_complement``, the D local
       aligner, CIGAR utilities.
   * - ``_prefilter``
     - ``--prefilter``
     - An exact 16-mer screen that rejects reads which cannot align, before MMseqs2 is invoked.
   * - ``_segmap``
     - ``--fast-segments``
     - Structure-aware chained seed-and-extend: the best V and best J per read with no homology
       search.
   * - ``_denoise``
     - used by ``correct``, ``resolve-ties``
     - The per-clonotype string sweeps: mean Phred, fraction below a threshold, substitution
       counting, and the substring test behind tie resolution.

Binary wheels ship all four. If they fail to build, arda falls back to a pure-Python markup path —
correct, and much slower — which is why ``setup.sh`` verifies the imports rather than trusting that
``import arda`` succeeded.

Sequence primitives (``translate``, ``detect_coding_frame``, ``reverse_complement``,
``back_translate``) are re-exported from ``arda.refbuild.translate`` with a mirpy-compatible API,
so mirpy can import arda and reuse them.

Reference build
---------------

``arda build-db`` regenerates ``database/vdj/<organism>/`` from IMGT. It is offline, reproducible,
and needed only to change the reference — annotation uses the committed one, and a plain
``pip install`` fetches it. The steps and every output file are documented in
:doc:`reference_build`; provenance is in ``SOURCES.md``.

Indexes are precompiled and shipped under ``database/vdj/<organism>/mmseqs/``, and used when the
local MMseqs2 version matches. An index is only reusable by the release that compiled it, so a
mismatch makes arda rebuild a private cache rather than accept it. ``segments.fasta``, the
collapsed 924-target reference the two-pass path screens against, is generated on demand under a
build lock.

See also
--------

* :doc:`usage` — which mode and configuration a given library needs.
* :doc:`outputs` — every column the pipeline above produces.
* :doc:`benchmarks` — what it costs and how accurate it is, with methods.
* :doc:`reference_build`, :doc:`reference_export` — the offline half, and how to get it out again.
