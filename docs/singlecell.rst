Single cell
===========

Per-cell contigs with no germline reference, chain pairing, and the doublets called.

``arda cells`` takes one UMI consensus per molecule with the cell barcode in the record name,
assembles each cell's contigs **with no germline reference**, annotates them, pairs the chains and
writes the diagnostics a droplet run is read through.

.. code-block:: bash

   # 1. upstream: demultiplex, correct barcodes, call cells, collapse molecules  (migec)
   migec checkout R1.fq.gz R2.fq.gz -o co/ --bc-pattern '^X{16}N{10}' --barcodes sheet.tsv
   migec refine co/PBMC_R2.fq.gz -o ref/ --cell-whitelist 737K-august-2016.txt
   migec assemble ref/PBMC.fq.gz -o asm/ --contig --min-reads 1

   # 2. arda: per-cell assembly, annotation, pairing, QC
   arda cells asm/PBMC.consensus.fq.gz -p out/PBMC \
       --cells ref/PBMC.cells.tsv --plot svg

arda does **no** barcode demultiplexing, **no** barcode correction, **no** UMI consensus and
**no** cell calling. The upstream tool does all four and leaves the barcode in the record name,
which :mod:`arda.cell` reads. ``--cells`` is how the called set gets in; without it every empty
droplet is assembled too.

Why a cell's molecules are enough
---------------------------------

Reads sharing a ``(cell, UMI)`` are co-terminal in 5' droplet chemistry, so a molecule's consensus
is a pile at one position: it covers one window of the transcript however deep it is. At
``--min-reads 30`` the mean migec consensus on ``sc5p_v2_hs_PBMC_1k`` is 204 nt against a ~508 nt
amplicon, and depth does not extend it.

What makes the full-length contig recoverable is that *different molecules start at different
positions*, so a cell's molecules tile the transcript. Scored directly, before any assembly runs:
**99.90%** of the 25-mers of Cell Ranger's 943 filtered contigs are already present in their own
cell's raw molecules, and 942 of 943 CDR3 nucleotide sequences appear verbatim in one of them.
The contig is in the data; assembling it is an engineering problem, not an information one.

The method
----------

Per cell, in :mod:`arda.singlecell`:

1. **Trim** the read-through adapter and any homopolymer tail off each molecule.
2. **Seed** on every 25-mer; wherever two molecules share one, propose the offset it implies.
3. **Verify** the implied overlap (``--min-overlap``, ``--min-identity``) before joining, and
   union-find the survivors into components.
4. **Phase** each component into haplotypes on columns that carry a real minor base, so two chains
   of one locus do not average into a third.
5. **Consense** each group by weighted column majority, the weight being the molecule's read
   depth.

Contigs are named ``<cell>_contig_<n>``, which is the Cell Ranger dialect, so any later
``arda`` run over the FASTA recovers ``cell_id`` for free.

Measured against Cell Ranger
----------------------------

``sc5p_v2_hs_PBMC_1k`` VDJ-T, both lanes, 479 Cell Ranger cells, 249,635 migec molecules.
Reproduce with ``migec/scripts/compare_cellranger_contigs.py``; the table is
``migec/assets/cellranger_contigs.tsv``.

.. list-table::
   :header-rows: 1

   * - variant
     - k-mer coverage
     - contigs at >=90%
     - CDR3 exact
     - contig N50
     - chain recall
     - doublets
   * - molecules, no assembly (ceiling)
     - 0.9990
     - 943 / 943
     - 942 / 943
     - --
     - --
     - --
   * - **arda cells, default**
     - **0.9759**
     - **892 / 943**
     - **933 / 943 (0.9894)**
     - 536
     - **0.9777**
     - 17
   * - no phasing (``--no-split``)
     - 0.9663
     - 868 / 943
     - 926 / 943 (0.9820)
     - 460
     - 0.9714
     - 12
   * - no adapter trim
     - 0.9525
     - 879 / 943
     - 907 / 943 (0.9618)
     - 580
     - 0.9491
     - 16

Per chain, at the default: **TRA 454/464 (0.9784), TRB 479/479 (1.0000)**.

.. warning::

   **Contig N50 is a description, never a score.** The no-adapter-trim row has the highest N50 in
   the table, 580 against 536, and the lowest accuracy in every other column. A contig built
   across an adapter is a longer contig and a wronger one.

Reading a cell: chains and doublets
-----------------------------------

``<prefix>.chains.tsv`` is one row per ``(cell, locus, junction)``, ranked within the cell by
supporting molecules. Rank 1 is ``primary``. A lower-ranked chain is believed only if it clears
all three gates, and then:

``doublet_candidate``
   a second **heavy** chain (TRB, TRD, IGH). Two rearranged heavy chains in one droplet is what a
   doublet looks like.
``secondary``
   a second **light** chain (TRA, TRG, IGK, IGL). Allelic inclusion at TRA and at the light-chain
   loci is a known, real population, not an artifact.
``extra``
   everything filtered out, kept in the table with the columns that say why -- ``productive``,
   ``molecules``, ``fraction_of_top``.

.. list-table:: the three gates
   :header-rows: 1

   * - option
     - default
     - what it is for
   * - ``--require-productive``
     - on
     - **the discriminator.** Of the extra chains Cell Ranger agrees with, 60/60 are productive;
       of those it does not, 20/128 are.
   * - ``--min-extra-molecules``
     - 3
     - absolute support. An extra chain on exactly one molecule is contamination 96-97% of the
       time and a real second TRB in none of the observed cases.
   * - ``--min-extra-fraction``
     - 0.2
     - support relative to the cell's top chain at that locus. Real second chains run at a median
       0.63 of the top; the ambient ones at 0.22.

Choosing them on your own data
------------------------------

Pass ``--reference`` with a call table you are willing to score against -- Cell Ranger's
``filtered_contig_annotations.csv``, a hashtag demultiplex, a simulation's truth -- and
``<prefix>.sweep.tsv`` reports precision and recall of the extra-chain filter per locus, with the
productivity gate off and on. That is the table the defaults were chosen from, and it is emitted
rather than baked in because the right setting is a property of your library's ambient load.

Before the productivity gate, on this library:

.. list-table::
   :header-rows: 1

   * - locus
     - min molecules
     - precision
     - recall
   * - TRA
     - 1
     - 0.360
     - 1.000
   * - TRA
     - 3
     - 1.000
     - 0.917
   * - TRB
     - 1
     - 0.034
     - 1.000
   * - TRB
     - 5
     - 0.889
     - 1.000

The distribution is bimodal, not a tail: the second chain is one molecule or it is a cell.

Agreement as a clustering, not as a list
----------------------------------------

Recall answers "did we find this chain". It cannot see a clone split in two or two clones merged
into one -- both leave every individual chain call correct. With ``--reference``,
``<prefix>.partition.tsv`` scores our per-cell clonotype partition against the reference's with
the homogeneity-parsimony trade-off (:mod:`arda.partition`): three families of two scores each,
one punishing merging and one punishing splitting.

.. warning::

   Read ``classes_singleton`` before the scores. An unexpanded repertoire is 476 clonotypes over
   479 cells, so there is almost no clustering to agree about and every score is decided by a
   handful of cells. These metrics earn their keep on an expanded repertoire or a simulation.

Where the curve breaks
----------------------

``arda cells`` reports the knee of the molecules-per-barcode curve in its JSON report and marks
the row in ``.cell_rank.tsv``. It **does not call cells** with it -- that is upstream's job -- but
it is what tells you whether the barcode set you handed it is a cell set or a cell set plus its
ambient tail.

The rule is Kneedle (`Satopaa 2011 <https://raghavan.usc.edu/papers/kneedle-simplex11.pdf>`_)
taken at the **global maximum** of its difference curve. On a decreasing curve, distance to the
chord joining the endpoints and Kneedle's normalised ``|y_n - (1 - x_n)|`` are the same quantity
up to a constant; only the selection rule differs.

.. warning::

   **Kneedle's published rule is degenerate without the smoothing spline the paper specifies.**
   It walks the *local* maxima and stops at the first that then falls by a sensitivity step.
   Measured on a 136,032-barcode 10x library:

   ==========================================  =======  ===============
   rule                                        rank     molecules
   ==========================================  =======  ===============
   Kneedle, unsmoothed (287 local maxima)      15       1,057
   Kneedle, smoothed over 0.01 log-rank        358      314
   Kneedle, smoothed over 0.10 log-rank        180      435
   **global maximum (what is implemented)**    **376**  **308**
   ==========================================  =======  ===============

   The unsmoothed answer would call fifteen cells. Properly smoothed it lands on the global
   maximum to within 5%, and over-smoothing drifts it off again -- so the global maximum reaches
   the same place with no smoothing parameter to tune.

.. warning::

   **A global maximum always exists, so the knee needs a floor.** An ambient-only library returns
   a rank too, and the threshold it implies calls every barcode a cell. ``find_knee`` refuses a
   knee below **ten times the mean molecules per barcode** -- one order of magnitude, the unit a
   log-log curve is read in. The mean alone was tried and is too weak: an ambient-only library of
   1-3 molecules per barcode is a step curve whose corner sits at 3 against a mean of 2.0 and
   clears a mean-only test at 1.5x. On the real library the knee is 308 against a mean of 3.38, a
   factor of 91. The two cases are two orders apart, so the boundary is not a tuned one.

   When the guard fires, ``knee_rank`` is 0 and the run says so. That is a real answer: it is what
   you get when the barcodes handed in were already filtered to called cells, because a cell set
   with its ambient tail removed has no knee left in it.

migec's ``refine`` computes the same knee in C++ beside its OrdMag call, and the two agree
exactly: rank 376, 308 molecules, on the same 136,032-barcode curve.

What is written
---------------

===========================  ================================================================
``.contigs.fasta``           the assembled per-cell contigs, ``<cell>_contig_<n>``
``.contigs.airr.tsv``        their AIRR annotation, with ``cell_id``, ``molecules``, ``reads``
``.chains.tsv``              one row per (cell, locus, junction), ranked, with a status
``.cells.tsv``               one row per cell, including the doublet scatter's two columns
``.cell_rank.tsv``           the knee curve: molecules per cell, ranked
``.contig_lengths.tsv``      the contig length histogram
``.sweep.tsv``               the filter's precision and recall (needs ``--reference``)
``.partition.tsv``           clustering agreement (needs ``--reference``)
``.arda.json``               the run report
``.stats.tsv`` / ``.json``   the QC table, in the same scopes a bulk run writes
===========================  ================================================================

The QC table
------------

``arda cells`` writes ``<prefix>.stats.tsv`` and ``.stats.json`` like every other mode, and in the
**same scopes** — ``sample``, ``chain`` per locus, ``v_gene``/``j_gene`` and ``junction_aa_len``.
"Which loci did this sample yield, and at what junction lengths" is the same question whether the
reads came one per cell or in bulk, so a cohort mixing single-cell and bulk samples has one table
to join on. What is genuinely single-cell — cells, molecules, the knee, the contig N50 — arrives
through the ``run`` scope from the report above, where it cannot be mistaken for a read count.

Three ratios are derived for you, because they are what a batch gets compared on:

``pairing_rate``
   ``cells_paired / cells`` — the AIRR Community's chain-pairing QC.
``doublet_rate``
   ``cells_doublet_candidate / cells``.
``molecules_placed_fraction``
   ``molecules_placed / molecules_in`` — how much of the input reached a contig.

A chain the extra-chain gate marked ``extra`` is **not counted as yield** anywhere in the table.
That gate exists because an extra chain supported by one molecule is ambient 96–97 % of the time,
and counting it would report the contamination as signal.

Several samples roll up with :doc:`arda qc batch <qc>` exactly as bulk ones do.

QC figures
----------

``--plot svg`` (or ``png``, ``pdf``) writes a gnuplot script per table and renders it if gnuplot
is on PATH. gnuplot is **not** a dependency: the ``.gp`` scripts and their tables are written
either way, so a figure can always be redrawn from the table beside it.

``cell_rank``
   the barcode rank plot, with the knee marked. **Molecules, never reads** -- one over-amplified
   molecule otherwise puts an empty droplet high on the curve, which is the artifact the plot
   exists to show. The knee row is marked in the table, not recomputed by the plot, so the figure
   and the number in the report cannot drift apart.
``doublet_scatter``
   the second heavy chain's support against the first, doublet candidates picked out. Only cells
   that *have* a second heavy chain are drawn; a cell with one has y = 0, which a log axis cannot
   show, and substituting a floor would invent the datum the plot is read for.
``chain_support``
   molecules per chain, rank 1 against rank 2+. The two are different distributions, which is the
   whole basis of the filter.
``contig_lengths``
   the assembled length histogram.
``filter_sweep``
   precision and recall against ``--min-extra-molecules``, with the productivity gate off and on.

Every panel addresses its columns by **name**, not by index. These tables gain columns, and a
positional ``$9`` then draws the wrong series with no error and a plausible figure.

Worked example
--------------

``notebooks/singlecell_qc.py`` is a marimo notebook that reads the tables above and reproduces
every panel plus the filter sweep interactively.

.. code-block:: bash

   uvx marimo edit notebooks/singlecell_qc.py

Limits
------

* **One UMI consensus per molecule is required.** Raw reads will assemble into something, but the
  column weights then mean read depth of an uncollapsed pile and the phasing gate sees PCR error
  as an allele.
* **Cell calling is upstream, and an unfiltered library is refused.** ``arda cells`` does not
  call cells. Without ``--cells`` it checks the molecules-per-barcode curve for a knee, and
  refuses when there is one well below the barcode count -- on ``sc5p_v2_hs_PBMC_1k`` the knee is
  at rank 376 of 136,032 observed barcodes, so 99.7% of what would be assembled is ambient RNA
  and it would land in ``.cells.tsv`` and the chain table looking exactly like cells. A plate or
  combinatorial library has no ambient tail, reports no knee, and is never refused.
* **One sample per file.** Droplet barcodes come from a fixed whitelist, so two samples share
  them by design; the grouping key here is the barcode alone, so a concatenated pair of consensus
  FASTQs would merge two cells into one and the result would read as an ordinary doublet. Two
  ``migec`` sample ids in one file are refused by name.
* **A saturated barcode breaks the premise.** If two molecules routinely share a ``(cell, UMI)``,
  a "molecule consensus" is a mixture and the layout is assembling chimeras. migec reports
  ``expected_molecules_per_group`` and warns.
* **Assembly is single-threaded Python**, about 54 ms per cell (479 cells in 26 s).
  ``--threads`` applies to the annotation pass. A process pool over cells is the obvious upgrade
  and was tried and reverted: on macOS the default ``spawn`` start method re-imports the caller's
  ``__main__`` in every worker, so a library function that opens a pool crashes any script or
  notebook without an ``if __name__ == "__main__":`` guard, and ``forkserver`` -- which does not
  touch ``__main__`` -- fails to start here at all. The upgrade path is the assembler in C++.
