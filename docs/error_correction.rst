Error correction
================

``arda correct`` turns a table of mapped reads into a table of clonotypes. Two things
happen there: reads are aggregated on the clonotype key, and CDR3 variants that are better
explained as sequencing error than as biology are collapsed onto the clone they came from.

.. code-block:: bash

   arda correct -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv

The clonotype key is ``(locus, v_call, j_call, junction)``. Abundance is the AIRR
``duplicate_count`` — every read that encompasses the junction, spanning or partial, assigned by
alignment — and ``consensus_count`` carries the distinct-fragment count beside it. Each
clonotype's D is mapped once, into its already-corrected junction, rather than voted over its
reads: D is a function of the junction, so calling it per read and taking a majority only adds
noise.

.. important::

   Pass ``--extra-airr`` whenever ``assemble`` ran. Without it the assembled long-CDR3 contigs
   are discarded and the clones that no single read spans never appear in the output.

The abundance test
------------------

A clonotype ``C`` is an **error child** of a more abundant neighbour ``P`` — differing by
``n_subs`` substitutions and ``n_indel`` inserted or deleted bases — when the number of misread
``P`` reads the error model *expects* is at least as large as the number of ``C`` reads actually
seen:

.. code-block:: text

   count[P] * p_sub**n_subs * p_ind**n_indel  >=  count[C]

with the per-event probabilities scaled by junction length ``L``:

.. code-block:: text

   p_sub = min(0.5, error_rate * L)
   p_ind = min(0.5, indel_rate * L)

If the inequality holds, ``C`` routes to ``P``; chains collapse to their ultimate ancestor.
Because ``p_err < 1``, counts strictly increase along parent pointers, so the parent graph
cannot contain a cycle. Where several neighbours qualify, the most abundant one wins, with a
deterministic tie-break on junction sequence.

Neighbours are found with :mod:`seqtree`, an edit-bounded index, so the cost is proportional to
the number of *near* neighbours rather than to the square of the clonotype count.

Why the rate is per base and the probability is not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--error-rate`` is a **per-base** rate (0.001 ≈ Phred 30). The collapse probability is
``error_rate * L``, because a single mismatch anywhere along a longer junction sheds
proportionally more error mass — there are more positions at which it could have happened. The
familiar vdjtools threshold of ~1/20 is calibrated for a 45 nt (15 aa) junction;
``error_rate = 0.001`` reproduces exactly that there (``0.001 * 45 ≈ 1/22``) and scales
correctly for junctions that are shorter or longer, which a flat 1/20 does not.

A multi-base indel costs ``p_ind ** len``, so a 3–9 bp in-frame indel — the somatic
hypermutation signature — is vanishingly unlikely as an instrument error and survives as a real
clonotype, while a 1 bp indel collapses. That asymmetry is deliberate and is why the indel term
is exponentiated in the length of the indel rather than counted as one event.

Counts are spanning reads: "2/2, not 2/200"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The counts entering the test are the **spanning** read depths — reads that fully observe the
junction — so the comparison is made only over the reads that actually saw the discriminating
position. A clonotype seen twice out of the two reads that reached its position is 100 % present
at that position, i.e. real; it is not a 1 % error of an abundant neighbour whose 200 reads never
covered that base. Testing it against the neighbour's total depth would delete it.

For very low coverage, ``--error-method binom`` or ``betabinom`` instead builds a per-position
read-depth pileup that also counts partial reads, buying depth at the cost of a distributional
assumption. ``simple`` (the default) is the spanning-count test above.

.. _quality gate:

Quality: the evidence abundance does not have
---------------------------------------------

The test above sees abundance and nothing else, so a sequencing miscall and a real low-frequency
variant are **the same object** to it — both are "a rare neighbour of an abundant clonotype" — and
``--error-rate`` can only trade them off globally. Phred separates them, because it is a different
measurement: a miscall is a detector artifact and reads low-Q, while a real base (a true variant,
or a template error made before the UMI) reads high-Q.

Measured at the discriminating base over 310,559 real MIGEC spike-in windows:

.. list-table::
   :header-rows: 1
   :widths: 44 20 18 18

   * - population
     - n bases
     - median Q
     - % below Q30
   * - parent clone, all 48 positions
     - 14,079,696
     - **38**
     - **5.1 %**
   * - EHEB-V1 (published, real)
     - 1,094
     - 35
     - 17.6 %
   * - EHEB-V2 (published, real)
     - 42
     - 34
     - 16.7 %
   * - 1-substitution error cloud
     - 12,281
     - **24**
     - **54.3 %**
   * - 2-substitution error cloud
     - 5,094
     - **6**
     - **91.1 %**

``--min-junction-q Q`` drops a read whose junction differs from its putative parent at **any base
below Q**.

.. code-block:: bash

   arda map     --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --junction-quality
   arda correct -i mapped.airr.tsv -o clones.tsv --error-rate 1e-5 --min-junction-q 20

Three properties of the gate, each deliberate:

* **Only the mismatching bases are evidence.** A junction that agrees with its parent everywhere
  else says nothing about whether the one differing base is real, so gating on the junction's
  *minimum* or *mean* quality asks the wrong question and mostly measures read length.
* **A clonotype with no more-abundant neighbour is never gated.** There is no hypothesis "this is a
  misread of X" to test, so nothing is dropped. A read whose quality string is missing or the wrong
  length is **kept** — absent evidence is not evidence of error.
* **It needs Stage 1 to have carried the quality, and raises if it did not.** ``map`` reads FASTQ
  quality only for ``--reconstruct``'s tie-break and otherwise discards it, so
  ``--junction-quality`` is what makes the gate possible. Without the ``junction_quality`` column
  ``correct`` **errors out** rather than silently not gating — an unapplied gate is
  indistinguishable from a gate that found nothing.

``map --junction-quality`` emits ``junction_quality``: the read's Phred+33 string over exactly the
bases of ``junction``, in the same orientation, so position *i* of one indexes position *i* of the
other. It is off by default (it appends a non-schema column, so the default output does not move),
costs **+2.2 % wall and +4.4 % bytes**, and is refused with ``--reconstruct`` — a merged fragment's
bases come from two reads, so no input quality string describes it.

.. note::

   For a ``rev_comp`` hit the quality belongs to the read as submitted while every coordinate is on
   the coding strand, so it is reversed and then **verified against the junction it claims to
   describe** before being written: a same-length slice off the wrong strand is a corruption
   nothing downstream could detect. Verified on 4,370 junction-bearing reads (4,365 of them
   reverse-strand) against an independent FASTQ extraction: 4,370 exact, 0 wrong.

What the gate does before any error model runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read counts on the MIGEC library. ``Err1``/``Err2`` are the most abundant error clonotype at the
same edit distance from the parent — the paper's own discriminating statistic.

.. list-table::
   :header-rows: 1
   :widths: 16 12 10 12 10 16 24

   * - ``--min-junction-q``
     - V1
     - V2
     - Err1
     - Err2
     - **V1/Err1**
     - distinct error clonotypes
   * - **0** (off)
     - **1,094**
     - **21**
     - 811
     - 76
     - 1.349
     - 2,690
   * - **20**
     - 941
     - 19
     - 446
     - 21
     - **2.110**
     - **300**
   * - 30
     - 901
     - 16
     - **421**
     - 21
     - **2.140**
     - 237
   * - 35
     - 612
     - 2
     - 342
     - 20
     - 1.789
     - 213

The effect is **flat over Q20–32 and degrades by Q35**, where it starts eating the real variant
(1,094 → 612 reads). Q20 keeps 86 % of V1's reads while removing 89 % of the distinct error
clonotypes.

Where the abundance model stops working: the ladder and the cliff
-----------------------------------------------------------------

The abundance test asks whether the parent could have produced this many misreads by chance. That
question is **only answerable for near neighbours**, and it is possible to measure exactly where it
stops being answerable.

If a clonotype *k* substitutions from its parent is accumulated independent error, the 1- and
2-substitution intermediates must **also** be observed — the same process generates them at far
higher rate. So the intermediates are a testable prediction. On a monoclonal T line (Jurkat,
dominant TRB clone 9,932 reads over 48 nt):

.. list-table::
   :header-rows: 1
   :widths: 8 14 24 16 38

   * - *k*
     - clonotypes
     - expected at ``1e-3``
     - median mean-Q
     - has an observed *(k−1)* neighbour
   * - 1
     - 108
     - 476.74
     - 31.4
     - —
   * - 2
     - 82
     - 11.20
     - 25.2
     - **82 / 82**
   * - 3
     - 28
     - 0.17
     - 24.2
     - 3 / 28
   * - 4
     - 13
     - 0.0019
     - 24.0
     - **0 / 13**
   * - 5
     - 18
     - 0.0000
     - 20.1
     - **0 / 18**
   * - ≥ 6
     - 14
     - 0.0000
     - 16.5–18.2
     - **0 / 14**

**k ≤ 2 is a ladder.** Every one of the 82 two-substitution variants has an observed
one-substitution neighbour on its path to the parent, and the observed counts are within an order of
magnitude of the binomial prediction. The abundance model is valid here, and chain collapse walks
it — which is why ``--max-subs`` defaults to 3 (2 plus headroom).

**k ≥ 4 is a cliff.** Zero intermediates at every *k*, and the model predicts 0.0019 clonotypes at
*k* = 4 where 13 are observed. These are not accumulated substitutions; they are single reads whose
whole junction window is unreliable, and the quality says so independently — median mean junction
Phred falls monotonically from 31.4 to 16.5, and the *k* ≥ 5 class is 100 % sub-Q30 against the
dominant clone's 5.9 %.

.. warning::

   Widening ``--max-subs`` to 10 *does* clean that library up (53 → 11 clonotypes, TRB purity
   .99990) — **for the wrong reason**. The abundance test it applies there has probability 0 to
   every printed digit. ``--error-rate`` is likewise inert from 1e-3 to 1e-1 on that class: the
   model has no neighbour to work with. If the class is to be collapsed, it must be collapsed on
   the evidence that actually distinguishes it, which is read quality — see the rescue below.

The quality-directed rescue
---------------------------

``--ec-mode amplicon`` and ``--ec-mode rnaseq`` add a second tier that reaches the cliff class.
Only clonotypes whose reads are **measurably bad** (median of the per-read mean junction Phred
below ``lowq_mean_q``) are considered, each may only join a neighbour at least ``lowq_min_ratio``
times more abundant within ``lowq_max_subs`` substitutions, and — the part that is not negotiable —
a candidate with **no** such neighbour keeps its reads.

.. important::

   **The rescue ignores the V and J calls on purpose, and honours the locus on purpose.**

   The abundance model above defaults to ``--require-vj``: a true sequencing error keeps the
   germline V/J call, which holds for the 1–3 substitution neighbours it collapses. The rescue
   targets the opposite class — the cliff, where the *whole junction window* is unreliable (median
   mean Phred 16.5–20.1). **A read that bad has an unreliable V/J call for the same reason its
   junction is unreliable**: the call came from aligning that sequence. Requiring the calls to
   agree would filter on the corrupted evidence. Measured on a full-depth TRA amplicon
   (SRR5233636), of 9,025 rescues under ``amplicon`` **4,593 (50.9 %) cross the V call** and 908
   the J; under ``rnaseq`` it is 24 and 2 of 215. What protects a genuine clone is not the call but
   the two gates that *are* trustworthy — its reads must be measurably bad, and the parent must be
   ``lowq_min_ratio`` times more abundant.

   The **locus** is not ignored: it is fixed by the whole read (V, J *and* C together), not by
   junction bases, so a locus flip is not a plausible consequence of miscalls, and a rearrangement
   of another locus is not a sequencing error of this one. Before the search was partitioned by
   locus, 3 of those 9,025 were 1-read **TRB** clonotypes absorbed into abundant **TRA** clonotypes
   at 11–12 substitutions — misattributed expression, not lost reads.

.. important::

   **Nothing in this framework discards a read.** A read that reached a complete junction came off
   a real rearrangement of that locus; deciding its junction carries a miscall is a statement about
   bases, not about whether the molecule existed. The rescue returns *parent assignments*, and the
   sum of ``duplicate_count`` is invariant — pinned by
   ``tests/unit/test_denoise.py::test_no_regime_or_key_ever_loses_a_read`` over every regime and
   both clonotype keys.

That is also why this is a rescue *radius* and not a quality *filter*. A whole-junction mean-Q
floor looks better on every cell line, but measured on a polyclonal hypermutated repertoire
(IGH_repertoire, 31,943 clonotypes) a floor at Q30 strands **3.70 %** of all junction-bearing reads
with no parent to inherit them — 47 % of everything it removes — against 0.148 % at a floor of 20.

Modes: ``--ec-mode``
--------------------

Four presets over the knobs above. An explicitly passed ``--error-method`` or ``--min-junction-q``
always wins over the mode.

.. list-table::
   :header-rows: 1
   :widths: 18 30 52

   * - ``--ec-mode``
     - is
     - use when
   * - ``fast`` (default)
     - ``--error-method simple``, no quality gate
     - Anything abundance-driven: repertoire overlap, diversity, clonal expansion. Byte-identical
       to arda's historical behaviour.
   * - ``accurate``
     - ``--error-method simple --min-junction-q 20``
     - Low-frequency variants matter: spike-ins, MRD-style tracking, a monoclonal control whose
       purity you are measuring. Needs ``map --junction-quality``.
   * - ``amplicon``
     - ``accurate`` + rescue at mean-Q < 25, radius 12 subs, ratio 50×
     - A targeted library, where a real clonotype is **deep**, so a 1-read neighbour of an abundant
       clone is almost always error and the rescue can be wide.
   * - ``rnaseq``
     - ``accurate`` + rescue at mean-Q < 20, radius 6 subs, ratio 200×
     - Bulk RNA-seq, where the receptor fraction is 0.02–3 % and **singletons are the norm and
       mostly real**, so the rescue stays narrow and demands a much stronger abundance ratio.

The two regimes differ because their clonotype-size distributions differ, not by taste. Measured on
Jurkat (14,531 reads, monoclonal), with ``--clonotype-key junction``:

.. list-table::
   :header-rows: 1
   :widths: 16 14 12 14 14 14

   * - ``--ec-mode``
     - clonotypes
     - reads
     - TRB purity
     - rescued
     - orphans
   * - ``fast``
     - 89
     - 14,531
     - .99562
     - 0
     - 0
   * - ``accurate``
     - 53
     - 14,531
     - .99696
     - 0
     - 0
   * - ``amplicon``
     - **10**
     - **14,531**
     - **.99990**
     - 43
     - 2
   * - ``rnaseq``
     - 40
     - 14,531
     - .99800
     - 13
     - 16

The read total does not move. On the MIGEC spike-ins all three published clonotypes survive every
mode (error clonotypes 1,630 → 108 under ``amplicon``) at exactly 310,559 reads.

The clonotype key: ``--clonotype-key``
--------------------------------------

The default key is ``(locus, v_call, j_call, junction)``. ``--clonotype-key junction`` canonicalises
V/J to the junction's majority first, so **call splits collapse**: a junction byte-identical to an
abundant clone's under a different V or J call. No error model can see that class — identical
junctions have no discriminating base — and on Jurkat it is the largest error class *by reads*
(130 of 14,531, including an allele-level ``TRGJ1*01`` / ``TRGJ1*02`` split).

Measured cost on a **polyclonal** TRA amplicon: 132 of 19,956 clonotypes merge (0.66 %), and in
every ambiguous case inspected the minority call carried **one** read against 4–10 for the majority,
on a short 30–39 nt junction where little V sequence sits inside the junction to call V from — a
call error on a low-abundance read, not a second clone.

Q20 is the **low end of the measured plateau**, not the optimum on any one library — on the MIGEC
data Q25–30 is marginally better on every axis. Shipping the conservative end is deliberate:
``accurate`` exists to protect rare real variants, and one library is not enough to tune a default
with.

.. warning::

   ``binom``/``betabinom`` are deliberately **not** in a mode, and "accurate" does not mean them.
   Measured on the same 302 k-read MIGEC library at ``--error-rate 1e-5``:

   .. list-table::
      :header-rows: 1
      :widths: 22 22 20 20 16

      * - arm
        - ``--error-method``
        - ``--min-junction-q``
        - wall
        - clonotypes
      * - MIGEC
        - ``simple``
        - 0
        - **3.32 s**
        - 1,633
      * - MIGEC
        - ``betabinom``
        - 0
        - 325.42 s
        - 1,633 (identical, 98×)
      * - MIGEC
        - ``binom``
        - 0
        - 238.68 s
        - 123
      * - MIGEC
        - **``simple``**
        - **20**
        - **3.79 s**
        - **127**
      * - MIGEC
        - ``binom``
        - 20
        - 216.33 s
        - 23
      * - Jurkat (monoclonal)
        - all three
        - 0
        - 0.34 / 0.38 / 0.39 s
        - **301 / 301 / 301**

   ``betabinom`` returns byte-for-byte the same answer as ``simple`` at 98× the wall — a pure loss.
   ``binom`` is a real gain on the deep library but at 57×, and it is an exact **no-op** on the
   monoclonal precision arm: every point of Jurkat purity came from the quality gate, none from the
   depth model. They stay explicit flags, for the very-low-coverage case they were written for.

What the options actually do
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - option
     - default
     - meaning
   * - ``--error-rate``
     - ``0.001``
     - Per-**base** substitution rate. Enters the test as ``error_rate * junction_len``.
       Lowering it collapses **less**; raising it collapses more. This is the knob to calibrate.
   * - ``--indel-rate``
     - ``0.001``
     - Per-**base** indel rate, length-scaled the same way, then raised to the power of the
       indel length.
   * - ``--max-subs``
     - ``3``
     - A search **radius**, not a threshold. It bounds how far ``seqtree`` looks for candidate
       parents; the length-scaled probability model still decides every collapse. Saturates at
       3 — beyond that the abundance test rejects the pairs anyway, so raising it only costs
       time.
   * - ``--max-indel``
     - ``0``
     - Indel bases searched. 1–2 bp indel errors are frameshifts that ``--complete-only``
       already removes, so this only matters with ``--all-junctions``; multi-bp SHM indels are
       kept regardless.
   * - ``--require-vj``
     - on
     - Only collapse neighbours that share V **and** J. A true sequencing error does not change
       the germline call, so a neighbour with a different V is a different clone.
   * - ``--complete-only``
     - on
     - Keep only junctions spanning Cys104 to [FW]118, in frame, no stop. A read that stops
       short of the anchor yields a *prefix of* a junction, not a clonotype.
   * - ``--error-method``
     - ``simple``
     - ``simple`` = spanning-read counts. ``binom``/``betabinom`` = per-position pileup for very
       low coverage — and ~270× slower, see the warning above.
   * - ``--ec-mode``
     - ``fast``
     - Preset over ``--error-method`` + ``--min-junction-q``. ``accurate`` = ``simple`` + Q20.
       An explicit knob overrides it.
   * - ``--min-junction-q``
     - off
     - Drop a read whose junction differs from its putative parent at any base below this Phred
       score. Needs ``map --junction-quality``; **raises** without it.
   * - ``--map-d``
     - on
     - Add the D columns per clonotype — ``d_call``/``d2_call``/``d_support``/``d2_support``, the
       D and V/J spans, and ``np1``–``np3`` — mapped into the corrected junction. Coordinates are
       1-based closed in **junction** space (``-1`` = not located).
   * - ``--d-max-evalue``
     - calibrated
     - E-value gate on those D calls. Lower is stricter; ``0.01`` is the band where D agrees .9985
       with IgBLAST on a TRB amplicon. See :doc:`d_segments`.

The one to reach for is ``--error-rate``. ``--max-subs`` looks like a stringency knob and is not
one: it changes which pairs are *considered*, never which are *accepted*.

Calibration: the MIGEC spike-ins
--------------------------------

The published MIGEC spike-in experiment (PRJNA239303) puts two known variants into a library at
known abundance, which makes it the natural test of whether an abundance-based corrector erases
real low-frequency clones.

**At the default ``--error-rate 1e-3``, arda erases both spike-in variants.** That is the
measurement, and it is worth stating precisely rather than filing as a bug, because the same
data says no abundance-based method can do better.

On the paper's own metric, computed on the raw reads before any correction:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - variant
     - ratio to its worst error neighbour
     - separable by abundance?
   * - V1
     - V1/Err1 = **1.35**
     - marginal
   * - V2
     - V2/Err2 = **0.28**
     - no

``V2`` is **less abundant than the worst 2-substitution PCR error in the same library**. No
abundance threshold — arda's, vdjtools', or anyone's — can separate a real variant from an error
that is more common than it is. This is the regime that UMI consensus exists for: collapsing
reads onto their molecular identifier first moves ``V1/Err1`` from 1.35 to **26–76**, and the
abundance test then works trivially. It is a property of signal-to-noise ~1 in the input, not of
the corrector.

What to do about it
~~~~~~~~~~~~~~~~~~~

Two measured settings, both actionable:

* ``--error-rate 1e-5`` recovers **both** spike-in variants exactly. Use it when low-frequency
  variant recovery is the goal and you can tolerate uncollapsed PCR error.
* ``--error-rate 1e-4`` — measured on an independent error cloud — kept both variants while
  still removing **72 %** of the real PCR errors. This is the useful middle.

``--error-rate`` is a **per-library calibration**, not a constant. It depends on the polymerase,
the cycle count, the platform and the read depth, and the default of ``1e-3`` is a Phred-30
sequencer assumption that says nothing about the PCR in front of it. If low-frequency clones
matter to the question, calibrate on a spike-in or on a known monoclonal control and pin the
value in the pipeline.

``--error-rate`` has a physical reading
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is a per-base error rate, so pick it from **what the input actually is**:

* **~1e-3 for raw reads.** Phred 30 — the sequencer's own substitution rate, which is what
  dominates when every read is an independent observation.
* **~1e-5 for UMI-consensus input.** Consensus over a UMI family removes the sequencing error
  almost entirely; what is left is the error that was already in the molecule.

Recovering the variants used to cost purity; with the gate it does not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reason ``1e-5`` was not simply the answer is that it also stops collapsing real PCR error. Two
arms, one for recall (MIGEC: both published variants must survive) and one for precision (a
monoclonal Jurkat TRB library: everything except the one true clone is spurious):

.. list-table::
   :header-rows: 1
   :widths: 34 16 20 16 14

   * - configuration
     - MIGEC variants
     - MIGEC error clonotypes
     - Jurkat TRB purity
     - Jurkat junctions
   * - ``1e-3`` (shipped default)
     - **0 / 2**
     - 0
     - .99540
     - 86
   * - ``1e-5``, no gate
     - 2 / 2
     - 1,630
     - **.96034**
     - 297
   * - ``1e-5 --ec-mode accurate`` (Q20)
     - **2 / 2**
     - **124**
     - **.99530**
     - **62**
   * - ``1e-5 --min-junction-q 35``
     - 2 / 2
     - 64
     - **.99600**
     - 59

Keeping both published variants used to cost 3.5 points of Jurkat purity. With the gate it costs
nothing: purity is back at or above the default's, while the default keeps **neither** variant.
Cross-lineage false positives are 0 throughout. The gate does not overrule ``--error-rate`` — at
``1e-3`` the variants are gone before it runs — so the two are set together.

.. code-block:: bash

   # low-frequency variant recovery, without paying for it in spurious clonotypes
   arda map     --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --junction-quality
   arda correct -i mapped.airr.tsv -o clones.tsv --error-rate 1e-5 --ec-mode accurate

The template-error floor: why V2 stays unrecoverable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Neither the gate nor UMI consensus rescues the second spike-in, and the reason is physical rather
than algorithmic.

Reverse transcription and the first few linear cycles of PCR happen **before the UMI is attached**.
An error made there is copied into every molecule of that UMI family, so it is in the consensus by
construction: UMI collapsing cannot remove it, and its Phred score is *high* — it is a real base,
faithfully read. At a template error rate of ~1e-5–1e-6 per base, a 48 nt junction carries a
spurious substitution in about **0.048 %** of molecules.

That is the floor. Against it:

.. list-table::
   :header-rows: 1
   :widths: 26 26 48

   * - variant
     - abundance vs parent
     - vs the ~0.048 % floor
   * - EHEB-V1
     - 0.373 %
     - ~8× above it — recoverable, and the gate recovers it
   * - EHEB-V2
     - 0.0072 %
     - **below it** — indistinguishable from template error by any method

So V2 is not a corrector failure and not a threshold that needs lowering. It is a variant present
at a frequency the library's own chemistry manufactures noise at. The way past it is experimental
— a higher-fidelity RT/polymerase, fewer pre-UMI cycles, or duplex tagging — not a parameter.

.. important::

   Do not chase a variant below the template-error floor by lowering ``--error-rate`` further. The
   1e-6 row of the same sweep keeps exactly the same 2 variants as 1e-5 while taking Jurkat TRB
   purity from .96034 to .94672 — it buys nothing and costs precision.

.. _what a junction disagreement means:

What a junction disagreement means
----------------------------------

.. important::

   **A V/J boundary disagreement inside the junction is not an error.**

   V(D)J recombination is probabilistic. Exonuclease chew-back removes a variable number of
   bases from the V, D and J ends, and N/P-nucleotide addition inserts untemplated bases in
   their place. The consequence is that for a large fraction of real junctions the partition
   into *V-end / NDN / J-start* **is not identifiable from the sequence alone** — several
   partitions explain the same nucleotides equally well, and the true one is unknowable. The
   ground truth does not exist to be recovered.

   So overlapping or disagreeing V / J / NDN assignments are **acceptable**, and two tools
   placing the V-end at different positions inside the same junction are not one right and one
   wrong. Never write such a disagreement up as an error rate.

What *is* checkable, and what arda is measured on:

* **The junction's outer bounds** — that it starts at Cys104 and ends at the [FW]118 anchor.
  These are conserved positions in the germline, not products of recombination.
* **The gene and allele calls** — ``v_call`` / ``j_call``, scored as tie lists (see
  :ref:`choosing a mode`).
* **Whether a tool invents a junction it has no anchor for.** If the Cys104 codon is not in the
  read, there is nothing to anchor on, and emitting a junction anyway is a defect. arda
  suppresses the junction in that case rather than extrapolating past its own alignment.

Those three are the ones to hold a corrector to. The internal partition is not.
