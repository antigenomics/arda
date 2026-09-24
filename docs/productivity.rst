Productivity: ``productive``, ``stop_codon``, ``vj_in_frame``
=============================================================

Three AIRR columns say whether a rearrangement could make a receptor. They are the most
misread columns arda writes, because each one is scoped differently and none of them means
"the junction looks like germline". This page states the rule arda applies and then works
through real reads where the rule is the only thing that settles the answer.

The rule
--------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - column
     - what arda puts in it
   * - ``vj_in_frame``
     - ``T`` when the V's reading frame, carried through the junction, arrives at FR4 on a
       codon boundary — so the J's FR4 and the constant exon spliced to it are translated in
       their own frame. Computed as ``phase = (fwr4_start - v_coding_start) % 3``, ``T`` iff
       ``phase == 0``.
   * - ``stop_codon``
     - ``T`` when any annotated region — FR1–FR3, CDR1–CDR2, the junction, **or FR4** —
       contains ``*``. Scoped to the annotated span, not to the whole read.
   * - ``productive``
     - ``T`` when ``vj_in_frame`` is ``T`` **and** no annotated region carries a stop. It is
       the conjunction, not an independent judgement.
   * - all three
     - **Empty** when the read never reached both the junction and FR4. That is not
       "non-productive", it is *unevaluable* — on real bulk RNA-seq it is about 72 % of mapped
       reads, and writing ``F`` there would look like a repertoire of broken rearrangements.

``vj_in_frame`` is a statement about **where the frame lands at FR4**, which is what decides
whether the constant region translates. It is not a statement about the junction's length, and
it is not a statement about whether the junction's 3′ residues match the called J's germline
residues. Those are different questions and they give different answers — see
`The two tests that disagree`_.

Flags, never filters
--------------------

arda writes these columns and does not act on them. No stage drops a non-productive read, and
no threshold is applied anywhere. A non-functional rearrangement is real biology: one allele of
a T cell's two is usually out of frame, and a sample with *no* out-of-frame reads is usually a
sample that has been filtered upstream.

Worked examples
---------------

Every sequence below is a GenBank cDNA entry annotated by ``arda annotate``. They are shown
translated from Cys104 in the frame arda used.

An ordinary productive read
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``MH918759.1`` — human TRA, ``TRAJ31*01``. The entry runs 66 nt past FR4, so the constant
region is in the read itself and needs no reconstruction:

.. code-block:: text

   CVVNIRGNARLMF GDGTQLVVKP NIQNPDPAVYQLRDSKSSDKSV
   └─ junction ─┘└── FR4 ──┘└──────── TRAC ────────┘

   vj_in_frame=T   stop_codon=F   productive=T

The frame runs from the V through the junction, through the canonical FR4 ``FGDGTQLVVKP``, and
on into TRAC without a stop. That is ``vj_in_frame=T`` by definition.

An indel inside the J
~~~~~~~~~~~~~~~~~~~~~

``JX277388.1`` — mouse TRB, ``TRBJ2-5*01``. This read is the instructive one, because a
germline-frame test and a functional test give opposite answers.

.. code-block:: text

   germline TRBJ2-5*01   AACCAAGACACCCAGTACTTT··GGGCCAGGCACTCGGCTCCTCGTGTTAG
   in its own frame       N  Q  D  T  Q  Y  F      G  P  G  T  R  L  L  V  L

   the read               ..CCAAGACACCCAGTACTTTTTGGGCCAGGCACTCGGCTCCTC..
                                                ^^
                                    two nucleotides germline does not have

The read carries two nucleotides between the anchor codon and FR4 that ``TRBJ2-5*01`` does not
have. So the J's 5′ germline residues and its own FR4 **cannot both** be read in the germline
frame — whichever one you put in frame, the other is shifted.

arda's frame puts FR4 in frame:

.. code-block:: text

   CTCSAGGPRHPVLF GPGTR
   └─ junction ──┘└ FR4 ┘     FR4 matches germline GPGTR exactly

   vj_in_frame=T   stop_codon=F   productive=T

Splice ``TRBC2*01`` onto it and translate: **0 stop codons over 375 nt**, yielding the real
mouse TRBC2 protein.

.. code-block:: text

   CTCSAGGPRHPVLFGPGTRLLVL EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSG...
                           └──────────────────── TRBC2 ────────────────────────────

The receptor translates. ``vj_in_frame=T`` is correct. A tool that reads the J's frame from a
lookup table keyed on the J's 5′ alignment start reports ``F`` for this read, because on that
side of the extra two nucleotides the germline residues are indeed shifted.

``JX277345.1`` — mouse TRB, ``TRBJ1-3*01`` — is the same shape with one extra nucleotide
instead of two:

.. code-block:: text

   CASSLLRQGEVEIRSIF GEGSR      FR4 matches germline GEGSR exactly
   + TRBC2*01  ->  CASSLLRQGEVEIRSIFGEGSRLIVVEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARG...

   0 stop codons.   vj_in_frame=T   stop_codon=F   productive=T

A truncated J answers nothing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``X60894.1`` — mouse TRA, ``TRAJ9*01``. The entry is 45 nt and stops three nucleotides into
FR4:

.. code-block:: text

   germline TRAJ9*01   AGGAACATGGGCTACAAACTTACCTTCGGGACAGGAACAAGC...
   the read              ..GAACATGGGCTACAAACTTAC·TTCGGG            <- read ends here
                                                ^ one nucleotide deleted

   arda:  CALSDGEHGLQTYF G
                         └ the only FR4 residue the read carries

The single FR4 codon present (``GGG``) is germline FR4's own first codon and sits on arda's
codon boundary, and completing the J from germline then splicing ``TRAC*01`` gives 0 stops over
261 nt. But **one residue of FR4 is not evidence**. With the J truncated this far, no tool can
answer the question from this read alone; treat the call as provisional rather than as a
disagreement between annotators.

.. _The two tests that disagree:

The two tests that disagree
---------------------------

There are two natural ways to ask "is the J in frame", and they are not the same test:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - test
     - what it actually measures
   * - **Does FR4 land on a codon boundary?**
     - Whether the constant exon translates. This is what ``vj_in_frame`` means and what arda
       computes.
   * - **Do the junction's 3′ residues match the called J's germline residues?**
     - Whether the read's J end resembles germline. This is *not* a frame test.

The second test is tempting because it needs nothing but the junction string, and it is wrong
often enough to be dangerous. Measured across arda's fixtures it fires on **51.8 % of IGL,
51.4 % of IGK and 39.8 % of IGH** rows (mean ``v_identity`` .955–.968) against **3.4 % of TRB
and 7.1 % of TRA** (mean ``v_identity`` .992). The gradient follows hypermutation, not frame:
it is detecting SHM that has reached the J, which is biology, not a defect. Do not use it as a
productivity check.

How often does this matter
--------------------------

Over 3,601 fixture rows where arda and IgBLAST both produced an evaluable call and agreed on
``junction_aa``, the two tools disagree on ``vj_in_frame`` for **8 rows (0.22 %)**. Every one
was examined read by read; each is a read carrying an indel inside the J, or a J truncated too
far to call. The column is stable — but when it does disagree, the disagreement is about which
of the two tests above is being run.

Checking a call yourself
------------------------

The self-contained check, for any read where the answer matters:

.. code-block:: python

   import polars as pl

   df = pl.read_csv("sample.airr.tsv", separator="\t", quote_char=None)
   row = df.filter(pl.col("sequence_id") == "MH918759.1").row(0, named=True)

   # arda's frame runs from Cys104; the junction and FR4 are contiguous in it
   print(row["junction_aa"], row["fwr4_aa"], row["vj_in_frame"], row["productive"])

If FR4 reproduces the called J allele's germline FR4 residues, the frame is the germline J's
own frame, and the constant exon spliced to that J is in frame too. That is the whole of
``vj_in_frame``.

.. seealso::

   :doc:`shm` for ``v_mutations`` / ``v_identity``, which is what the second test above is
   really measuring, and :doc:`qc` for the per-sample counts of out-of-frame and stop-codon
   reads.
