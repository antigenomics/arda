Error correction
================

``arda rnaseq correct`` turns a table of mapped reads into a table of clonotypes. Two things
happen there: reads are aggregated on the clonotype key, and CDR3 variants that are better
explained as sequencing error than as biology are collapsed onto the clone they came from.

.. code-block:: bash

   arda rnaseq correct -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv

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
       low coverage.
   * - ``--map-d``
     - on
     - Add ``d_call``/``d2_call``/``d_support`` per clonotype, mapped into the corrected
       junction.

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
