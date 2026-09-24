Recombination scenarios
=======================

``arda scenarios`` reads nucleotide junctions and estimates the generative model of V(D)J
recombination they came from — the trimming and insertion distributions, and which D goes with
which J.

Why it exists
-------------

:mod:`arda.dpost` places and identifies a D from an *amino-acid* junction by marginalising a
generative model. That model ships as ``database/vdj/<org>/d_prior.tsv``, and every number in it
is **borrowed from OLGA / vdjrearm**. This command is how arda estimates its own, from its own
output:

.. code-block:: bash

   arda rnaseq --r1 R1.fq.gz --r2 R2.fq.gz -o results -p sample
   arda scenarios -i results/sample.clones.tsv -o my_prior.tsv --organism human

The output is the same long ``locus / kind / key / value`` table the shipped file uses, with the
same key grammar, so it is a **drop-in** for it.

.. list-table::
   :header-rows: 1
   :widths: 16 26 58

   * - kind
     - key
     - what it is
   * - ``insVD``
     - ``<n>``
     - P(non-templated nt between V and D)
   * - ``insDJ``
     - ``<n>``
     - P(non-templated nt between D and J)
   * - ``dlen``
     - ``<d_allele>:<n>``
     - P(surviving D length | allele)
   * - ``d_marginal``
     - ``<d_allele>``
     - P(D)
   * - ``d_given_j``
     - ``<d_allele>|<j_allele>``
     - P(D | J)
   * - ``delV`` / ``delJ``
     - ``<allele>:<n>``
     - P(germline nt deleted | allele) — **new**, nothing shipped these

A scenario is not identifiable, and that is the whole design
------------------------------------------------------------

Given a junction and its V/J calls, a *scenario* is the tuple that reproduces it exactly:

.. code-block:: text

   junction == Vg[:len(Vg)-delV] + N1 + Dg[delDl:len(Dg)-delDr] + N2 + Jg[delJ:]
                                   insVD                           insDJ

Several tuples reproduce the same junction. A first non-templated base that happens to match
germline is indistinguishable from one less nucleotide of trimming — the same ambiguity that
makes a V/J boundary inside a junction not a defect to be fixed. On one real human TRB junction
there are **4,346** admissible scenarios.

So the counts are **expected counts under the current model, summed over scenarios**, never the
counts of one reading. Counting a single MAP reading assigns weight 1 to one member of the
ambiguity set and 0 to the rest, which biases every trimming and insertion distribution toward
less trimming and shorter inserts.

Having a weight requires a model; a model plus expected counts plus renormalisation is EM. The
log-likelihood is echoed per pass:

.. code-block:: text

   scenarios: iteration 1/4 -- 2357 records, log-likelihood -67159.5
   scenarios: iteration 2/4 -- 2357 records, log-likelihood -61008.0
   scenarios: iteration 3/4 -- 2357 records, log-likelihood -60477.4
   scenarios: iteration 4/4 -- 2357 records, log-likelihood -60331.0

.. warning::

   **An insertion costs its own sequence, not just its length.** The junction's 5′ end reads
   either as templated V or as an insertion that happens to match V germline — and with a length
   term alone the second is *free*, so EM walks straight into the corner where everything is
   insertion. Measured before ``P(seq | len) = 0.25^len`` was added: three iterations on 503 real
   human TRB junctions moved ``insVD`` mass onto **10–11 nt**, with the log-likelihood rising
   monotonically the whole way. With the term, the same data gives ``insVD`` peaking at **4 nt**
   and ``dlen`` for ``TRBD1*01`` peaking at **4–5 surviving nt** — which independently reproduces
   the *"median surviving D is 5 nt for human TRB"* figure measured in :mod:`arda.dpost`.

What it reads, and how it is weighted
-------------------------------------

Any table with ``junction``, ``v_call`` and ``j_call`` — a ``.clones.tsv`` from ``correct`` or an
AIRR TSV from ``map``.

``--weight abundance`` (default) weights each record by ``duplicate_count``; ``--weight rows``
weights every row by 1. Neither is right for every question — a clonotype row is a clonotype, not
an observation of the recombination process, but a clonal expansion is one recombination event
seen many times — so the choice is written into the output header rather than defaulted silently.

Every locus contributes. VJ loci (TRA, TRG, IGK, IGL) have no D and no ``insDJ``, but they carry
most of the trimming evidence, so they are not skipped.

.. note::

   ``dlen == 0`` — the D trimmed away entirely — is enumerated and is a common outcome, not a
   failure. The shipped table's single largest entry is ``IGHD1-1*01:0 = 0.8759``. Dropping
   records where no D survived would truncate ``dlen`` at 1 and inflate every insertion
   distribution by the length of the D that was really there.

Limits
------

* **Exact matching, so unmutated receptors.** A hypermutated IGH junction breaks the premise that
  germline segments appear verbatim. TR is the intended input; IG works where SHM has not reached
  the junction.
* **Generating a prior is not adopting one.** ``database/vdj/<org>/d_prior.tsv`` is unchanged;
  swapping in an estimate is a measurement and a release decision, not a side effect of running
  this.
* **The insertion composition is uniform** (0.25/base) rather than a fitted first-order Markov
  chain. The per-base cost is what breaks the degeneracy; composition is a refinement.

Scoring a junction: ``arda.hmm``
--------------------------------

The same recursion, used the other way round. :func:`arda.scenarios.lattice` **is** the
forward-backward pass of a semi-Markov model of V → N1 → D → N2 → J; an E-step and a posterior
differ only in what you do with the same term weights, so there is one implementation.

.. code-block:: python

   from arda.hmm import model_for, log_likelihood, posterior_d, best_scenario

   m = model_for("human")                       # or model_for("human", prior="my_prior.tsv")
   log_likelihood(junction, "TRBV20-1*01", "TRBJ2-1*01", m)   # -40.03
   post = posterior_d(junction, "TRBV20-1*01", "TRBJ2-1*01", m)
   post.probabilities   # {'TRBD1*01': 0.001, 'TRBD2*01': 0.010, 'TRBD2*02': 0.989}
   best_scenario(junction, "TRBV20-1*01", "TRBJ2-1*01", m)    # the Viterbi path, for reading

**Semi-Markov, not Markov**, because germline deletions and insertion lengths have explicit,
tabulated, non-geometric durations — a plain HMM would impose geometric ones on quantities that
are measured not to be. **Conditioned on the V and J call** that mmseqs already made, which is
what keeps the state space to ``(delV, insVD, D, delDl, delDr, insDJ, delJ)`` and makes it cheap
per clonotype after ``correct``, never per read.

``log_likelihood`` is the marginal, not the best path: a junction explainable many mediocre ways
is more probable than one explainable a single slightly-better way, and a Viterbi score cannot
say so. ``model_for(prior=...)`` reads a table written by ``arda scenarios``, so a model
estimated on one cohort can score another.

.. warning::

   This does **not** replace :mod:`arda.dpost`, and it does not gate anything. ``dpost`` answers
   the *amino-acid* question — a record with no nucleotides, where the D is often invisible in
   the translated junction. Different input, both ship. Nothing in the annotation path calls
   ``arda.hmm``: two measured negatives (recorded in ``ROADMAP.md``) say that re-ranking
   nucleotide D candidates by a scenario likelihood changes nothing, and that replacing the
   E-value gate with a Bayes factor would need a *per-locus* shipped threshold.

Inspecting one junction
-----------------------

:func:`arda.scenarios.enumerate_scenarios` returns the whole set for one junction, unweighted —
which is what inspection wants, and what the estimator sums in factorised form:

.. code-block:: python

   from arda.scenarios import enumerate_scenarios

   found = enumerate_scenarios("TGCAGTGCTAGAGATTTAGCGGGAGGGACTGAAGCTTTCTTT",
                               "TRBV20-1*01", "TRBJ2-1*01", "human")
   len(found)                     # 4346
   sorted({s.d_call for s in found})
   # ['TRBD1*01', 'TRBD2*01', 'TRBD2*02']

.. seealso::

   :doc:`d_segments` for how a D is called on a read, and ``project/design-scenarios.md`` for the
   design record.
