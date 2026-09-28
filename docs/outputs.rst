Output files and columns
========================

What a run writes, and what is in it. This is the reference to reach for when a column's meaning
matters.

.. contents::
   :local:
   :depth: 2

Reading arda's files
--------------------

All of arda's tabular output is tab-separated, uncompressed, with a single header line and
**quoting turned off**.

.. important::

   Read it with quoting off too. A reader using the usual default turns a blank field into the
   two-character value ``""``, which silently makes every clonotype look chimeric and every empty
   flag look set.

   .. code-block:: python

      import polars as pl
      df = pl.read_csv("SAMPLE.clones.tsv", separator="\t", quote_char=None)

   .. code-block:: python

      import pandas as pd
      df = pd.read_csv("SAMPLE.clones.tsv", sep="\t", quoting=3)   # csv.QUOTE_NONE

Coordinates are **1-based and closed**, in query space, as both AIRR and GFF3 are. A region's
``*_start`` and ``*_end`` both name bases that belong to the region, so
``sequence[start-1:end]`` is that region's sequence. Anything 0-based half-open must convert.

An empty field means **not evaluated**, never zero and never false. A locus no read spanned has no
minimum junction length rather than a minimum of zero; a read that never reached both the junction
and FR4 has an empty ``productive`` rather than ``F``.

Files a run writes
------------------

``arda rnaseq`` and ``arda amplicon`` write these under ``--out-dir``, named from
``--out-prefix`` (or from the sample id when one is given):

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - File
     - Contents
   * - ``<prefix>.airr.tsv``
     - One AIRR Rearrangement row per **mapped read**. Stage 1's output. The largest file, and the
       one to keep if you want to re-derive anything.
   * - ``<prefix>.clones.tsv``
     - One row per **clonotype**. Stage 2's output, with Stage 3's assembled clones folded in.
       This is the table most analyses start from.
   * - ``<prefix>.assembled.airr.tsv``
     - Stage 3's assembled long-CDR3 contigs, as AIRR rows. Written when ``--assemble`` is on.
   * - ``<prefix>.arda.json``
     - The **run report**: what was read, what mapped, why the rest did not, wall time, peak RSS,
       and the versions of arda, MMseqs2 and the reference.
   * - ``<prefix>.stats.tsv`` / ``.stats.json``
     - **Run QC** in long format. The numbers that decide whether a sample is usable, without
       re-reading the FASTQ.

``arda cells`` writes a different set, described under :ref:`outputs-singlecell`.

.. _outputs-airr:

``<prefix>.airr.tsv`` — one row per mapped read
-----------------------------------------------

A spec-valid `AIRR Rearrangement <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`_
file; it passes ``airr.schema`` validation. The same schema is what ``arda annotate`` writes.

Identity and locus
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Meaning
   * - ``sequence_id``
     - The read or record name, verbatim from the input.
   * - ``sequence``
     - The query sequence, re-oriented to the plus strand if it mapped reversed.
   * - ``locus``
     - ``TRA``, ``TRB``, ``TRG``, ``TRD``, ``IGH``, ``IGK`` or ``IGL``.
   * - ``rev_comp``
     - ``T`` when the read mapped on the minus strand and was re-oriented.
   * - ``cell_id``
     - The cell barcode, present when ``--cell-from`` or ``--cell-regex`` was given.

Gene calls
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Meaning
   * - ``v_call``, ``j_call``
     - IMGT allele names. A comma-joined **tie list** when the alignment cannot separate several
       germlines — standard AIRR practice, and an exact statement of what the read supports.
   * - ``d_call``, ``d2_call``
     - The D germline, and the second D of a tandem D-D. VDJ loci only; empty on a VJ locus,
       which is not a miss.
   * - ``d_support``, ``d2_support``
     - The Karlin–Altschul E-value the D call was accepted on. Lower is better. Re-thresholdable
       with ``--d-max-evalue``. See :doc:`d_segments`.
   * - ``c_call``
     - The constant-region gene, from a hit on a :term:`J+C scaffold` or from a read's
       constant-region mate.
   * - ``c_class``
     - The IGH isotype **class** — ``IGHG``, ``IGHM``, ``IGHA`` — never the noise-prone subclass.
   * - ``v_call_genotyped``
     - Added by ``resolve-ties --genotype``. ``v_call`` is left byte-identical. See
       :doc:`genotype`.

Regions
~~~~~~~

For each region ``r`` in ``fwr1``, ``cdr1``, ``fwr2``, ``cdr2``, ``fwr3``, ``cdr3``, ``fwr4``:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Meaning
   * - ``<r>_start``, ``<r>_end``
     - 1-based closed coordinates in query space.
   * - ``<r>``
     - The region's nucleotide sequence.
   * - ``<r>_aa``
     - Its translation.

A region absent from the read is empty rather than guessed: a V-only query yields ``fwr1``
through ``fwr3`` and nothing else, and a J-only query yields ``fwr4`` alone. There is no coverage
filter.

The junction
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Meaning
   * - ``junction``, ``junction_aa``
     - Cys104 through Phe/Trp118 **inclusive of both anchors**. This is AIRR's ``junction`` and
       MiXCR's ``CDR3``.
   * - ``cdr3``, ``cdr3_aa``
     - The IMGT CDR3: the same stretch **excluding** both anchors, two residues shorter.
   * - ``np1``, ``np2``, ``np3``
     - The non-templated stretches: V→D, D→D2, D→J.
   * - ``v_anchor_nt``, ``j_anchor_nt``
     - The anchor bases the called germlines template into each end of the junction.
   * - ``junction_completed_nt``
     - How many bases were imputed from the germline to finish a truncated junction, under
       ``--complete-junctions``. Zero or empty when nothing was imputed.

.. warning::

   ``junction`` and ``cdr3`` are not interchangeable. VDJdb's ``cdr3`` column holds what arda
   calls ``junction``. Conflating them shifts every coordinate by one residue at each end.

A junction is emitted only when the read actually reaches its Cys104 anchor. A read that lands
mid-V and has no anchor gets an empty junction rather than a shortened one.

Alignment and mutation
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Meaning
   * - ``sequence_alignment``, ``germline_alignment``
     - The aligned query and germline strings.
   * - ``v_cigar``, ``d_cigar``, ``d2_cigar``, ``j_cigar``, ``c_cigar``
     - Per-segment CIGAR against the germline.
   * - ``v_sequence_start``, ``v_sequence_end``, ``j_sequence_start``,
       ``d_sequence_start``, ``d_sequence_end``, ``d2_sequence_start``, ``d2_sequence_end``
     - Segment boundaries in query space.
   * - ``v_germline_start``, ``v_germline_end``, ``j_germline_start``, ``j_germline_end``,
       ``d_germline_start``, ``d_germline_end``, ``d2_germline_start``, ``d2_germline_end``
     - The same boundaries in germline space. Filter on the V pair when you need to know how much
       V a read actually covered.
   * - ``v_identity``
     - Fraction identity to the called V germline, scoped **outside the junction** by default.
       ``--shm both`` adds the junction-inclusive ``v_identity_full``.
   * - ``v_mutations``, ``j_mutations``
     - Per-segment substitutions in germline coordinates, comma-joined. See :doc:`shm`.
   * - ``v_mutation_quality``, ``j_mutation_quality``
     - The Phred behind each mutation entry, comma-joined **as integers**, one for one with the
       mutation list. Requires ``map --mutation-quality``.
   * - ``junction_quality``
     - Phred over exactly the bases of ``junction``, as **raw Phred+33 characters**. Requires
       ``map --junction-quality``.

.. note::

   The two quality columns use different encodings on purpose: ``junction_quality`` is a Phred+33
   string, ``v_mutation_quality`` is comma-joined integers. Reading one as the other yields
   plausible numbers off by 33.

Productivity
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Column
     - Meaning
   * - ``vj_in_frame``
     - ``T`` when the V's frame, carried through the junction, arrives at FR4 on a codon boundary.
   * - ``stop_codon``
     - ``T`` when any annotated region — including FR4 — carries a ``*``. Scoped to the annotated
       span, not the whole read.
   * - ``productive``
     - The conjunction of the two.

All three are **empty** when the read never reached both the junction and FR4 — *unevaluable*
rather than non-productive. On real bulk RNA-seq that is around 72 % of mapped reads. They are
flags: no stage drops a non-productive read. Worked examples in :doc:`productivity`.

Search diagnostics
~~~~~~~~~~~~~~~~~~

Columns prefixed ``mmseqs2_`` carry the raw search result — score, E-value, identity, and query
and target spans. They are diagnostics rather than AIRR fields; ignore them unless you are
debugging a call.

.. _outputs-clones:

``<prefix>.clones.tsv`` — one row per clonotype
-----------------------------------------------

Twenty-one columns, in this order, plus ``umi_count`` when the input supports it:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Column
     - Meaning
   * - ``junction``, ``junction_aa``
     - The clonotype's junction, anchors included.
   * - ``v_call``, ``j_call``, ``c_call``, ``locus``
     - The calls the reads agreed on. ``c_call`` is voted per fragment, then per clonotype.
   * - ``duplicate_count``
     - **Reads encompassing the junction.** arda's abundance estimate. Invariant across every
       ``--ec-mode``: correction moves reads onto a parent, it never discards them.
   * - ``consensus_count``
     - **Distinct fragments.** Use this when you need molecules rather than abundance.
   * - ``umi_count``
     - Distinct UMIs. Present only under ``correct --cell-from``, and omitted — never 0 or 1 —
       when the dialect names no UMI, so a default run's table has 21 columns and not this one.
   * - ``d_call``, ``d2_call``, ``d_support``, ``d2_support``
     - As in the AIRR table.
   * - ``np1``, ``np2``, ``np3``
     - The non-templated stretches.
   * - ``v_sequence_end``, ``d_sequence_start``, ``d_sequence_end``,
       ``d2_sequence_start``, ``d2_sequence_end``, ``j_sequence_start``
     - Segment boundaries **within the junction**.

.. note::

   A V/J boundary disagreement *inside* the junction is not an error. Exonuclease chew-back and
   N/P addition mean the V-end / N-D-N / J-start partition is often not identifiable from sequence
   alone, so the ground truth for it is unknown. What is checkable is the junction's outer bounds,
   the gene calls, and whether a junction was invented without an anchor.

Grouping is controlled by ``--clonotype-key`` (``full`` = locus, V, J, junction; or ``junction``
alone) and ``--call-level`` (``allele`` or ``gene``).

.. _outputs-report:

``<prefix>.arda.json`` — the run report
---------------------------------------

Top level: ``arda_version``, ``mmseqs_version``, ``reference`` (path, size, mtime),
``wall_seconds``, and one object per stage.

.. code-block:: json

   {
     "arda_version": "2.30.1",
     "mmseqs_version": "18-8cc5c",
     "wall_seconds": 3.394,
     "map": {
       "total_reads": 1320,
       "mapped_reads": 453,
       "mapped_fraction": 0.3431818181818182,
       "per_locus": {"IGH": 104, "IGK": 107, "IGL": 104, "TRB": 88, "TRA": 45, "TRG": 5},
       "read_length_mean": 100.0,
       "paired": true,
       "threads": 4,
       "peak_rss_mb": 352.1,
       "unmapped": {
         "prefilter_rejected": 788,
         "no_hit": 8,
         "constant_only": 54,
         "below_min_score": 17,
         "accounted": 1320
       }
     },
     "assemble": {"seeds": 58, "contigs": 25, "contigs_complete": 10, "reads_rescued": 6},
     "correct": {"clonotypes_in": 49, "clonotypes_out": 49, "reads": 57, "collapsed": 0}
   }

Two fields deserve naming.

``unmapped`` **accounts for every read that produced no row**, so a low mapped fraction has a
cause: rejected by the k-mer prefilter, no search hit, a constant-region-only hit, or a hit below
``--min-score``. ``accounted`` equals ``total_reads``.

``fast_fraction`` appears under ``map`` when the two-pass search ran. It is the fraction of reads
that hit both a V and a J segment, and it predicts whether the amplicon speed configuration pays
better than the library's name does.

.. _outputs-stats:

``<prefix>.stats.tsv`` — run QC
-------------------------------

Four columns — ``scope``, ``key``, ``metric``, ``value`` — one value per cell, so a metric can be
grepped, ``join``-ed across samples or plotted without reshaping.

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - ``scope``
     - ``key``
     - What it holds
   * - ``run``
     - ``map`` / ``correct`` / ``assemble``
     - The run report verbatim: total and mapped reads, FASTQ bytes, read length, paired,
       threads, wall time, peak RSS.
   * - ``sample``
     - —
     - Library-wide totals, junction lengths and quality, SHM rate, V/J gene coverage.
   * - ``chain``
     - ``TRB``, ``IGH``, …
     - Per locus, reads **and** clonotypes: functional, non-functional, stop codons, truncated
       junctions, junction length min/max/mean, junction quality, SHM rate, chimeras.
   * - ``v_gene`` / ``j_gene``
     - ``TRBV19``
     - Reads and clonotypes per germline gene.
   * - ``allele_candidate``
     - ``TRBV19*01:G45A``
     - A recurrent, high-quality V mutation, with its frequency and mean Phred. A shortlist, not a
       genotype call.
   * - ``junction_aa_len``, ``read_len``, ``clone_size``, ``isotype``, ``chain_support``
     - ``IGH:17``
     - The distributions, keyed ``locus:bucket``.

.. code-block:: text

   $ awk -F'\t' '$1=="chain" && $2=="IGH"' SAMPLE.stats.tsv
   chain  IGH  reads                       104
   chain  IGH  reads_truncated_junction    1
   chain  IGH  junction_nt_mean            48.75
   chain  IGH  junction_quality_mean       35.6718
   chain  IGH  shm_rate                    0.0413
   chain  IGH  clonotypes_chimeric         2

The distribution scopes exist because a mean cannot show a **shape**, and shape is most of the
diagnosis. A bimodal junction length is two primer sets in one tube; a read-length cliff is an
adapter left on; a clone-size distribution with no singletons is a library amplified before it was
sequenced. Each of those reads as an unremarkable mean.

Cohort roll-up and the dashboard are in :doc:`qc`.

.. _outputs-singlecell:

``arda cells`` output
---------------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - File
     - Contents
   * - ``<prefix>.contigs.fasta``
     - The assembled per-cell contigs.
   * - ``<prefix>.contigs.airr.tsv``
     - Their AIRR annotation, with ``cell_id``.
   * - ``<prefix>.chains.tsv``
     - One row per ``(cell, locus, junction)``, ranked, with a status.
   * - ``<prefix>.cells.tsv``
     - One row per cell, including the two columns the doublet scatter uses.
   * - ``<prefix>.stats.tsv`` / ``.stats.json``
     - The same QC scopes as the bulk modes, so one cohort table holds every mode.

Column-level detail and the doublet criteria are in :doc:`singlecell`.

Reference exports
-----------------

``arda export-ref`` writes the reference itself, in three kinds × four formats — scaffolds,
collapsed per-allele segments, or per-allele CDR3 anchors, as TSV, FASTA, GFF3 or AIRR. Because
coordinates are already 1-based closed, GFF3 passes through unchanged. See
:doc:`reference_export`.
