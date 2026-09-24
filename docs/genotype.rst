Personalized germline: genotype and call restriction
====================================================

A germline reference is a catalogue of every allele anyone carries. No donor carries all of them,
and no donor carries more than two per gene — so a call naming a third is wrong before any sequence
is looked at, and a read that "could not choose" between two alleles the donor does not have was
never ambiguous at all. Immcantation's TIgGER measured the size of this on full-length BCR:
restricting V calls to an inferred genotype took ambiguous assignments from 11.2 % to 1.5 %.

arda splits this into the two halves that have very different confidence:

``arda resolve-ties --genotype``
    **Apply** an allele set you already trust. Solid, cheap, and useful on its own.
``arda genotype``
    **Infer** the allele set from mapped reads, by a likelihood ratio between diploid genotypes.
    Whether it can say anything depends on read length — see below.

.. _genotype-restriction:

Applying a genotype
-------------------

The simplest possible genotype file is a one-column TSV of allele names:

.. code-block:: text

   allele
   TRBV19*01
   TRBV19*03
   TRBV20-1*01

.. code-block:: bash

   arda resolve-ties -i sample.airr.tsv -o sample.genotyped.tsv --genotype genotype.tsv

This adds one column, ``v_call_genotyped``, and leaves ``v_call`` byte-identical — the evidence the
genotype was applied to stays in the same file, so a reader can re-restrict at a different genotype
without re-running anything.

**It never re-aligns and never rebuilds a reference.** Given the span a read already aligned over,
the restriction is a set intersection: the germlines a read cannot rule out, intersected with the
ones its donor carries. The same reasoning as TIgGER's ``reassignAlleles``, plus three arda-specific
reasons a per-donor reference would be a trap — scaffold ids are positional so changing the allele
set renumbers every scaffold in the locus; ``arda build-db`` needs IgBLAST and IMGT network access
that a ``pip install`` does not have; and the mmseqs freshness contract is an mtime with no
allele-set identity recorded, so a donor-specific FASTA beside the shared one would silently
invalidate the index for every concurrent process.

``v_call_genotyped`` has exactly three value classes, and they are deliberately distinguishable:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - value
     - meaning
   * - narrower than ``v_call``
     - the usual case, and the point
   * - equal to ``v_call``
     - either the tie machinery declined (span under 30 nt, allele not in the reference, tie list
       past the cap) **or** no candidate's gene was genotyped. A refusal to answer is not a
       contradiction; treating the two alike would delete the call of every short read.
   * - empty
     - every candidate belongs to a genotyped gene and none is carried — this read contradicts the
       genotype

An allele survives unless **its own gene** was genotyped and did not name it. Tie lists routinely
span genes, so a gene-blind test empties every read whose tie list merely brushes a genotyped gene.

Any source works: an OGRDB set, a library MiXCR inferred, a list typed by hand. Every named allele
is validated against the reference and an unknown one **raises** — 884 human V alleles are
functional in IMGT and 801 reach a scaffold, so naming one of the others is an easy mistake that
would otherwise restrict nothing, silently.

Inferring a genotype
--------------------

.. code-block:: bash

   arda genotype -i sample.airr.tsv -o sample.genotype.tsv --loci TRB

**The junction is de facto a UMI.** V(D)J junctional diversity makes a nucleotide junction
essentially unique to one rearrangement, so grouping reads by ``(locus, gene, J gene, junction)``
groups them by molecule. Two things follow:

* **The clonotype, not the read, is the unit of observation** — each distinct junction is one
  independent draw from the donor's two chromosomes, however many reads carry it. Read depth is
  reported beside the clonotype count, never substituted for it.
* **Disagreement within a junction is error, not allele.** Every read of one rearrangement carries
  the same V allele by construction, so a read naming a different one is hypermutation or
  sequencing error. That is where the error rate comes from — measured on the library rather than
  picked. ⚠ It is measured over *all* reads, not the germline-exact subset used for assignment:
  those were selected for carrying no mismatch, so they never disagree and the estimate collapses
  silently onto its floor (0 discordant of 77,345 reads on the amplicon below).

The call is then a **likelihood ratio between diploid genotypes**. Every single allele and every
pair is scored by its multinomial likelihood under that error rate, and the winner is called only
if it beats the runner-up by ``--min-log10-bf`` (default 1.0 — ten times more likely). Homozygous
and heterozygous are the same formula at genotype size 1 and 2, so neither is special-cased.

.. note::

   **Why a likelihood and not a coverage rule.** The first version of this was TIgGER's frequency
   rule — the fewest alleles explaining 7/8 of the calls — which has no error model and so cannot
   tell overwhelming evidence from none. On a real TRB amplicon it called ``TRBV11-2`` off **754 of
   757** clonotypes that singled the allele out and ``TRBV20-1`` off **1 of 2,544**, reported both
   as ``explained = 1.0000, ok``, and returned 43 of 53 genes with not one heterozygous. A ratio of
   counts is not a confidence.

Ambiguous clonotypes still constrain the answer. Discarding them is not the conservative option, it
is wrong: if a donor is ``*01/*02`` and only ``*02`` is separable at the read length, every
unambiguous observation is ``*02``, the likelihood says homozygous, and the restriction then
empties every ``*01`` read — systematically, on every sample with that genotype. A genotype that
fails to intersect ``--fraction-to-explain`` of the compatible sets is rejected however well it
fits the counts.

Read length is the binding constraint
-------------------------------------

Separating a gene's alleles needs, measured from the 3' end of the V germline:

.. list-table::
   :header-rows: 1

   * -
     - multi-allele genes
     - median nt needed
     - separable within 100 nt of the 3' end
   * - TRBV
     - 44 / 56
     - 150
     - 16 / 44
   * - TRAV
     - 30 / 45
     - 175
     - 8 / 30
   * - IGHV
     - 57 / 75
     - 230
     - 9 / 57

Measured end to end on ``SRR5233641`` (human TRB amplicon, 151 nt paired, 99,839 mapped reads →
**33,440 clonotypes / 77,345 reads**, error rate **5.25 × 10⁻⁴**): reads cover a median of **72 nt**
of V germline — **56 nt** after clipping at the Cys104 anchor — and **not one** reaches the ~150 nt
TRBV needs. Only **33.1 %** of clonotypes can be assigned an allele at all. The result is **17 of 53
genes called**, of which 11 are ``single_allele`` (one catalogued allele, carried without being
inferred) and **6 genuinely inferred**, with log₁₀ Bayes factors of 223 (``TRBV11-2``, 754
clonotypes), 251 (``TRBV5-6``, 846) and 11.7 (``TRBV5-8``, 39). The other 36 are refused —
including ``TRBV10-3`` with 1,037 clonotypes, none of which can separate its alleles.

That is the honest output for that library, and the refusals are the point: a full-length, 5'RACE
or ``arda cells`` library (contig N50 536 nt) has the resolution this one does not.

The output has one row per carried allele, and a row with an empty ``allele`` for every gene that
could not be called. ``clonotypes`` and ``reads`` are that allele's coverage at both levels;
``log10_bf`` is the margin over the runner-up genotype.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``note``
     - meaning
   * - ``ok``
     - called, beating the runner-up by at least ``--min-log10-bf``
   * - ``single_allele``
     - the reference has one allele for this gene; carried without being inferred
   * - ``low_support``
     - fewer than ``--min-clonotypes`` clonotypes could be assigned an allele at all
   * - ``undetermined``
     - the best genotype does not beat the runner-up by ``--min-log10-bf`` — the data cannot
       separate them
   * - ``unexplained``
     - the best genotype fits the counts but leaves too many ambiguous clonotypes unexplained

What this is not
----------------

**A claim about the reference, never about the repertoire.** It answers which germline alleles a
donor carries — annotation, the same kind of question as which J a read uses. Clonal composition,
diversity and overlap are ``vdjtools``'.

**Not** :doc:`qc`'s ``allele_candidate`` scope, which stays a shortlist of recurrent high-quality V
mutations to look at and is deliberately never a call.

Downstream, ``vdjtools`` already ingests a user germline through
``model.reference.read_germline_fasta`` and ``model.extend_alleles``, whose own docstring names the
use case: *a newer IMGT release, a population-specific library, your own genotyped alleles*.
