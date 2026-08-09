Somatic hypermutation
=====================

A B-cell repertoire is not a set of germline calls. After activation, AID mutates the rearranged V
region at ~10⁻³ per base per division, and the resulting lineage — who descends from whom — is the
signal most BCR analyses are actually after. arda records that per read, in the one coordinate
frame a downstream tool can use directly.

What arda records
-----------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - field
     - what it is
   * - ``v_mutations`` / ``j_mutations``
     - The read's substitutions against its called germline: ``G45A,C112T`` — germline base,
       1-based position **in that segment's own germline allele**, read base. Empty when the read
       is germline over that segment.
   * - ``v_identity``
     - The same information as a fraction (identity over the V alignment).
   * - ``v_cigar`` / ``j_cigar`` / ``c_cigar``
     - Per-segment AIRR CIGARs. **Indels live here**, as ``I``/``D``.
   * - ``sequence_alignment`` / ``germline_alignment``
     - The full aligned strings. Everything above is derived from them — see the warning below
       about deriving it yourself.
   * - ``v_germline_start`` / ``v_germline_end``, ``j_germline_*``
     - The germline window the read actually covered. Mutation positions lie inside it.

The mutation lists come out of the same single walk of the alignment that builds the CIGARs, so
they are close to free: **+36 ms per 100,000 mapped reads** (measured on 35,825 real bulk IG
alignments, 26.9 → 39.7 ms, CIGARs byte-identical before and after) and **+2.25 %** on the output
TSV (21,461,770 → 21,944,646 bytes). Amplicon wall clock does not move.

Why not just diff the two alignment strings
-------------------------------------------

The information was never missing — ``sequence_alignment`` and ``germline_alignment`` carry every
column, and the germline arda reports matches the shipped allele on **28,365 of 28,365** mapped
reads of a real bulk IG library (66,526 V mismatches, zero disagreements). What it was not, was
*usable*.

.. warning::

   arda aligns to a ``V + N-pad + J [+ C]`` **scaffold**, not to a germline. A consumer that does
   the obvious AIRR thing — diff the two alignment strings — gets 100,091 mismatches on that
   library of which **20,140 (20.1 %) are N-pad or constant-region columns**. It attributes
   junction positions to a germline. Recovering it correctly needs arda's scaffold geometry
   (``mmseqs2_tstart``, ``mmseqs2_t_vend``, ``mmseqs2_t_jstart``, ``mmseqs2_t_vjend``); the two
   strings are also 30.7 % of the TSV.

   Use ``v_mutations`` / ``j_mutations`` in preference to the alignment strings — but read the
   defect below first.

.. danger::

   **⛔ KNOWN DEFECT (open): ``v_mutations`` / ``j_mutations`` DO include junction-internal
   positions.** This page previously claimed the opposite — that the N-pad is not a segment, so a
   junction position "cannot enter the list by any code path". That claim is false and is
   retracted. The pad is excluded, but the **V germline's 3' tail and the J germline's 5' head lie
   inside the junction**, and mutations are scoped by segment (``t <= t_vend``,
   ``[t_jstart, t_vjend]``), not by the junction boundary. So exonuclease chew-back and
   non-templated N/P bases are emitted as substitutions against a germline that does not template
   them.

   Measured on a **TRA amplicon** (SRR5233636, 500,000 reads) — T-cell receptors do **not**
   somatically hypermutate, so every entry there is spurious by construction:

   * **1.046 V and 1.658 J "mutations" per read**;
   * **86.2 % of J entries sit at J germline position <= 10**, i.e. the 5' head, inside the
     junction;
   * splitting the reported load by the frequency of each ``(allele, position, alt)`` across reads
     carrying that allele: **13.0 % of J entries are below frequency 0.01** (sequencing error),
     **80.8 % fall between 0.01 and 0.5** (junction diversity), and **6.2 % sit at 0.5 or above**
     (a genuinely wrong *allele* call, e.g. ``TRAV8-6*01`` positions 281/282 at 0.88, ``TRAJ8*01``
     position 1 at 0.67).

   Until this is fixed, treat the lists as *three superimposed populations* and separate them by
   frequency: ``>= 0.5`` is an allele artifact, ``< 0.01`` is sequencing error, and the remainder
   on TR loci — or at J positions below the ``[FW]118`` anchor on any locus — is junction
   diversity, not SHM. ``v_identity`` has the same scope defect: it runs to ``t_vend``, so it is
   depressed by junction diversity rather than by mutation load.

.. important::

   **A mutation inside the junction is not attributable to any germline**, and arda should not
   claim one.

   V(D)J recombination trims the V, D and J ends by a variable amount and inserts non-templated
   N/P bases in their place, so the *V-end / NDN / J-start* partition of a junction frequently is
   not identifiable from the sequence at all. A "mutation" placed inside the NDN is a statement
   about a boundary that has no ground truth — see :doc:`d_segments` for the same rule on the D
   side. ⚠ That is the design intent; the defect above is that the implementation does not yet
   honour it.

Two more properties worth knowing:

* **Substitutions only.** An SHM indel is one event of unbounded length, not a per-position call;
  it stays in the CIGAR as ``I``/``D``. Germline coordinates on the far side of an indel are still
  correct, because the walk tracks the target position across gap columns.
* **A read** ``N`` **is a no-call, not a mutation.** A column contributes only when both bases are
  unambiguous ACGT.
* **Positions are on the coding strand.** For ``rev_comp = T`` the ``sequence`` field is the read
  as submitted, so a read-side lookup (a Phred quality, say) must be made against the reverse
  complement of it.

Accuracy against IgBLAST
------------------------

3,000 hypermutated IGH reads (``v_identity < .97``), IgBLAST truth at ``v_score >= 70``, per-mutation
precision and recall over the whole mutation set:

.. list-table::
   :header-rows: 1
   :widths: 34 14 16 16 20

   * - arm
     - reads
     - precision
     - recall
     - identical mutation set
   * - shares a V allele with IgBLAST
     - 1,947
     - .98991
     - .98926
     - 1,880 / 1,947 = .9656
   * - **\+ germline frame verified**
     - 1,882
     - **.99902**
     - **.99834**
     - **1,876 / 1,882 = .9968**

The 65-read gap between the arms is a **scoring artifact, not an error**: IgBLAST answers with an
allele tie list and aligns to one member of it, tied members differ in 5′ length, so a shared
*gene* can be two different coordinate frames (``33N`` vs ``105N`` on IGHV4-30-2). The frame arm
additionally requires IgBLAST's own germline string at its own coordinate to *be* arda's allele
there. The 6 residual frame-verified disagreements are gap placement — 251 of the 3,000 reads carry
a V indel.

Building a lineage tree
-----------------------

A tree builder needs three things: a **root** (the germline), **member sequences in one coordinate
frame**, and the **clone membership** that says which reads belong together. arda gives all three,
but they come from two different stages.

.. code-block:: bash

   arda rnaseq map     --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv
   arda rnaseq correct -i mapped.airr.tsv -o clones.tsv --read-map read_map.tsv

``--read-map`` writes ``sequence_id -> corrected junction``, which is the join key: it links each
read's ``v_mutations`` back to the error-corrected clonotype it was assigned to. Replay the
mutations onto the clone's common V germline window and you have a gapless alignment whose root is
the germline — the input format ``dnapars``, ``gctree`` and ``IQ-TREE`` read directly.

On a real bulk IG library that yields **66 IGH clones with ≥ 5 reads, 735 member sequences, 70
parsimony-informative sites, 0 malformed alignments**, the largest clone 88 reads over 37 distinct
mutation sets.

.. warning::

   **Clone membership requires a complete junction, and most hypermutated reads do not have one.**
   On that library, **8,345 of 9,755 mutated IGH reads (85.5 %) carry SHM but join no clonotype** —
   they are V-only fragments that never reach the Cys104 anchor, so ``--read-map`` cannot see them
   by construction. That is exactly the population with the most SHM signal per read.

   If the question is the *mutation spectrum* (targeting, selection, replacement/silent ratios),
   score ``v_mutations`` over **all** mapped reads and do not go through the clonotype table. If the
   question is *lineage*, accept the coverage cost — a tree needs the junction to know what a clone
   is.

Two more cautions when consuming this:

* An **allele tie list** shares one reference target named with the comma-joined list. Look the full
  ``v_call`` string up first; matching on the first member silently changes the coordinate frame.
* Count over **clonotypes, never reads**, when inferring a donor's germline variants from apparent
  mutations. One expanded clone otherwise makes its own SHM look germline.

See also
--------

* :doc:`usage` — the full AIRR column list.
* :doc:`error_correction` — a multi-base indel survives correction on purpose; it is the SHM
  signature, not an instrument error.
* :doc:`d_segments` — the same non-identifiability rule, on the D side.
