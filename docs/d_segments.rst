D segments and tandem D-D
=========================

The D segment is the hardest call in a rearrangement. It is 10–37 nt of germline to begin with,
both of its ends are chewed back by exonuclease before joining, and what survives sits inside the
junction surrounded by non-templated N/P bases. arda therefore does not report a D as a fact: it
reports a D **with the E-value that accepted it**, and gives you the knob to move the gate.

.. code-block:: bash

   arda rnaseq map -i reads.fq -o mapped.airr.tsv                    # shipped gate, E <= 0.2
   arda rnaseq map -i reads.fq -o mapped.airr.tsv --d-max-evalue 0.01  # the strict band
   arda rnaseq correct -i mapped.airr.tsv -o clones.tsv --d-max-evalue 0.01
   arda annotate  -i mrna.fasta -o out.airr.tsv --d-max-evalue 0.05
   arda rnaseq run --r1 R1.fq.gz --r2 R2.fq.gz -p S -d out/ --d-max-evalue 0.01

``--d-max-evalue`` is accepted by ``annotate`` and by every ``rnaseq`` stage that maps D
(``map``, ``correct``, ``assemble``, ``run``, ``reduce``), and by
:func:`arda.annotate.dmap.map_d_junction` in the library.

How the call is made
--------------------

D mapping runs only on **VDJ loci** (IGH, TRB, TRD) and only inside the **V..J interior** — the
stretch of the junction the V and J germlines do not template. The interior is aligned gaplessly
against every allowed D germline of the locus, and the best segment is accepted when its
Karlin–Altschul E-value clears the gate. A second, non-overlapping segment is then sought in what
is left; if it also clears the gate *and* passes the orientation constraint below, the record
carries a tandem **D-D**.

The columns:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - field
     - meaning
   * - ``d_call`` / ``d2_call``
     - Allele call(s), comma-joined when the germlines tie byte-for-byte. ``d2_call`` is the
       **3′** segment of a D-D fusion, not "the runner-up" — the two are sorted by position on
       the read, not by score.
   * - ``d_support`` / ``d2_support``
     - The E-value that accepted that call. **Shipped so you can re-threshold offline** without
       re-running: it ranks correctness (see the bands below).
   * - ``np1`` / ``np2`` / ``np3``
     - Non-templated stretches: V→D1, D1→(D2 or J), and D2→J. ``np3`` is empty without a D-D.
   * - ``d_sequence_start/end``, ``d2_sequence_start/end``
     - 1-based closed coordinates of each D. In **read** space in the per-read AIRR; in
       **junction** space in the clonotype table, where there is no read (``-1`` = not located).
   * - ``d_germline_start/end``, ``d_cigar``
     - Where in the D allele the surviving piece came from. Empty on ``--seqtype aa`` on purpose:
       the alignment offsets there index a reading frame, not the D germline.

``d_support`` ranks correctness: the E-value bands
--------------------------------------------------

Measured against an IgBLAST truth at gene level (IgBLAST at ``v_score >= 70``), sweeping the gate
on two real libraries. The trade is call rate against agreement, and it is monotone in both.

**TRB amplicon** (SRR5233641; 45,604 reads with a projected V..J interior, 31,608 IgBLAST D calls):

.. list-table::
   :header-rows: 1
   :widths: 18 14 12 16 16 12

   * - ``--d-max-evalue``
     - D called
     - call rate
     - judged
     - gene agreement
     - tandem D-D
   * - **0.20** (shipped)
     - 18,362
     - .4026
     - 18,106
     - .9765
     - 5
   * - 0.10
     - 14,690
     - .3221
     - 14,571
     - .9842
     - 2
   * - 0.05
     - 10,749
     - .2357
     - 10,707
     - .9911
     - 2
   * - **0.01**
     - 5,355
     - .1174
     - 5,344
     - **.9985**
     - 0

**Bulk IGH** (SRR5233639, full 1.32 M pairs; 1,795 reads with an interior, 1,056 IgBLAST D calls):

.. list-table::
   :header-rows: 1
   :widths: 18 14 12 16 16 12

   * - ``--d-max-evalue``
     - D called
     - call rate
     - judged
     - gene agreement
     - tandem D-D
   * - **0.20** (shipped)
     - 648
     - .3610
     - 206
     - .9417
     - 0
   * - 0.10
     - 557
     - .3103
     - 158
     - .9494
     - 0
   * - 0.05
     - 461
     - .2568
     - 118
     - .9746
     - 0
   * - **0.01**
     - 347
     - .1933
     - 77
     - **1.0000**
     - 0

Read this as a **recall/precision dial, and pick it from the question**:

* **Repertoire-level D usage, VDJ-model fitting, anything summing over many reads** — keep 0.2.
  It is the highest-recall setting, and its errors are not systematic enough at scale to move a
  usage histogram much.
* **Per-clonotype D annotation you will act on, or a tandem D-D you intend to report** — use
  ``0.01``. That is the band where the call agrees with IgBLAST on essentially everything it
  makes, at roughly a third of the calls.
* The shipped 0.2 is deliberately the **weakest** band in both libraries. It ships as the default
  because dropping two thirds of the D calls is the wrong default for a repertoire tool, not
  because it is the accurate one.

.. note::

   Zero D calls come out of the VJ loci (TRA, TRG, IGK, IGL) — they have no D germlines and the
   locus gate is exact, not statistical.

.. warning::

   A shuffle null puts the **false-D rate at .0884 in TRB** against a real call rate of .4025 at
   the shipped gate, so roughly a fifth of TRB single-D calls at ``E <= 0.2`` are reproducible
   from composition alone. IGH is far cleaner (.0095) because its D germlines are longer. This is
   the same statement as the table above, from the other side: a TRB D at 0.2 is a ranked
   hypothesis, not an identification.

Germline geometry: two constraints, applied before the statistics
-----------------------------------------------------------------

Both of these refuse products that **recombination cannot make**. Neither is a score threshold,
and neither can be tuned away with ``--d-max-evalue``.

Orphons cannot rearrange
~~~~~~~~~~~~~~~~~~~~~~~~

IMGT ships ``/OR`` D genes — ``IGHD.../OR15-...`` lie on **chromosome 15**, outside the IGH locus.
They are not disfavoured; they are not producible, and they are excluded unconditionally
(``transfer._allowed_d``). This is not academic: on a real bulk IGH library, **11 of 11 tandem D-D
calls named** ``IGHD2/OR15-2a*01,IGHD2/OR15-2b*01`` **as their second D**, so the library's entire
tandem-D signal was this one vocabulary artifact. Excluding them leaves 639 of 650 single-D calls
untouched, moves 9 to a rearrangeable gene, loses 2, and takes IGH tandem D-D from 11 to 2.

TRBD2 cannot join a TRBJ1
~~~~~~~~~~~~~~~~~~~~~~~~~

The TRB locus runs ``TRBD1 – TRBJ1 cluster – TRBC1 – TRBD2 – TRBJ2 cluster – TRBC2`` and V(D)J
joining deletes the DNA between the joined segments, so TRBD2 sits downstream of the entire J1
cluster and can never reach it. Without this rule, TRBD2 (16 nt) simply outscores TRBD1 (12 nt) on
noise: 17 % of real human TRB J1-cluster records were assigned an impossible TRBD2, at a median
E-value of 0.096 — the chance band — against 0.014 for genuinely producible TRBJ2 × TRBD2. An
ambiguous J spanning both clusters excludes nothing.

A tandem D-D must run in genomic order
--------------------------------------

D-D fusion is a rearrangement like any other: the upstream D joins the downstream D and everything
between is deleted, so the fused product carries them in **genomic 5′→3′ order**. A read whose 5′ D
lies 3′ of its 3′ D in the germline locus names a product deletional joining cannot make, and so
does the same gene twice (that needs two germline copies).

``transfer._dd_orientation_ok`` enforces it over ``_D_GENOMIC_ORDER`` — ``TRBD1 < TRBD2`` and
``TRDD1 < TRDD2 < TRDD3``, the architectures pinned independently of species. It is applied **after**
the two segments are sorted by position on the read (the rule is about order on the read, not about
which scored higher), and a refused pair collapses to the single higher-scoring D — arda does not
reach further down the score list to manufacture a producible partner. Where either side is an
allele ambiguity list, one producible assignment is enough.

On the TRB amplicon this takes tandem calls from **15 to 5**:

.. list-table::
   :header-rows: 1
   :widths: 30 18 18 20

   * - pair
     - before
     - after
     - producible?
   * - TRBD1 → TRBD2
     - 5
     - **5**
     - yes
   * - TRBD2 → TRBD1
     - 7
     - **0**
     - no
   * - TRBD2 → TRBD2
     - 3
     - **0**
     - no

It removes **none** of the producible calls (5 of 5 kept) and the single-D call count is identical
either way (18,362) — it deletes only the second call, never the first.

.. warning::

   **Fixing the orientation does not make TRB tandem D-D real.** Under a composition-preserving
   shuffle (100 permutations, conditioned on a real first D, D1's span held fixed, only the
   flanking non-templated bases permuted): with the gate on, **5 observed against 2.71 expected**,
   per-permutation range 0–7, Poisson ``p(X >= obs) = .139``. The pre-fix ``p = .0031`` (15 observed,
   6.53 expected) is not evidence either — the "excess" it measured *was* the 10 genomically
   impossible calls, which a flank-only shuffle cannot generate and therefore under-counts. What
   the gate buys is that the residual signal is composed only of producible pairs; it does not buy
   significance. On bulk IGH the null is empty on both sides (0 observed, 0.15 expected).

   Report a TRB tandem D-D only from the strict band (``--d-max-evalue 0.01``, where the count on
   this library is 0), or from long reads / assembled contigs where the D-D is actually spanned.

.. note::

   **IGH is deliberately absent from** ``_D_GENOMIC_ORDER``. In *human* IMGT the second number of
   ``IGHD<family>-<position>`` is the genomic position, but in *mouse* it is a family-member index
   with no locus meaning — and the two vocabularies collide on real gene names (``IGHD1-1``,
   ``IGHD2-15``, ``IGHD5-5``, ``IGHD5-12``, ``IGHD6-6`` exist in both). The D mapper is handed
   sequences, not an organism, so a name-parsed IGH rank would silently mis-order mouse. A gene
   absent from the table imposes no constraint.

Using the output: D-D markup
----------------------------

Naming a second D and giving no way to find it is not markup. Every stage now emits the
coordinates and the non-templated stretches alongside the calls, and the clonotype table
(``arda rnaseq correct``) carries them in **junction space**, 1-based closed, with ``-1`` for
"not located":

``v_sequence_end``, ``d_sequence_start``/``d_sequence_end``, ``d2_sequence_start``/
``d2_sequence_end``, ``j_sequence_start``, ``np1``, ``np2``, ``np3``.

The partition closes exactly::

   np1 + D1 + np2 + D2 + np3 == junction[v_sequence_end : j_sequence_start - 1]

verified on **every** record carrying a ``d2_call`` in both local libraries (5 of 5 at read level
and 4 of 4 at clonotype level on the TRB amplicon; bulk IGH has none).

From the library, :meth:`arda.annotate.dmap.DCall.markup` hands you the cut directly — labelled
parts that concatenate back to the junction byte-for-byte:

.. code-block:: python

   from arda.annotate.dmap import map_d_junction

   jn = "TGTGCTCTTGGGCCCCGGCCTTCCTACAGCGAGGAGTTGGGGGATACCCATCGGGCCGATAAACTCATCTTT"
   call = map_d_junction(jn, "TRDV1*01", "TRDJ1*01", "human", d_max_evalue=0.01)
   call.d2_call                      # 'TRDD3*01' — a tandem D-D
   call.markup(jn)
   # [('V', ...), ('np1', ...), ('TRDD2*01', ...), ('np2', ...),
   #  ('TRDD3*01', ...), ('np3', ...), ('J', ...)]
   "".join(s for _, s in call.markup(jn)) == jn        # True, always

Without a D the markup is ``[("V", …), ("N", …), ("J", …)]``; with a single D the ``np3`` and second
D entries are absent. It returns ``[]`` when the V/J split could not be located at all.

.. important::

   **The boundaries inside the junction are one consistent reading, not ground truth.**

   V(D)J recombination is probabilistic: exonuclease chew-back removes a variable number of bases
   from the V, D and J ends and N/P addition inserts untemplated bases in their place, so the
   *V-end / np1 / D / np2 / J-start* partition frequently is **not identifiable from the sequence
   at all** — several partitions explain the same nucleotides equally well. This binds hardest on
   D, which is trimmed at **both** ends inside the NDN.

   So a boundary here disagreeing with another tool's is not an error in either, and neither is a
   ``d_sequence_start`` that moves by a base or two between arda versions. Do not score it, and do
   not pin one in a test. What *is* checkable, and what the tables above score: the **gene/allele
   call**, the junction's **outer bounds** (Cys104 … [FW]118), whether the partition **closes**,
   and whether a call is **producible** at all.

See also
--------

* :doc:`error_correction` — D is mapped once per clonotype, into the already-corrected junction,
  rather than voted over reads.
* :doc:`usage` — the full AIRR column list and the aa-input behaviour.
* ``examples/dd.airr.tsv`` — the two real human reads (of 7,341, across five organisms) that carry
  a tandem D-D, regenerated by ``python examples/regenerate.py``.
