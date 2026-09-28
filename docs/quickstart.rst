Quickstart
==========

This page takes you from an empty environment to a clonotype table you can read. It uses a small
FASTQ pair committed to the arda repository, so it runs in a few seconds and needs no data of
your own. Every number shown below is that run's actual output.

.. contents::
   :local:
   :depth: 1

1. Install
----------

.. code-block:: bash

   pip install arda-mapper

That is the whole installation. arda imports as ``arda``, and the binary wheels carry its C++
extensions. Two things it needs at runtime — the MMseqs2 search binary and the curated germline
reference — are fetched automatically into ``~/.cache/arda`` the first time you run a command
that needs them.

Check what arda resolved:

.. code-block:: bash

   arda info

If you are setting up a development checkout instead, see :doc:`installation`.

2. Run a library
----------------

Clone the repository for the test data, or substitute your own FASTQ pair:

.. code-block:: bash

   git clone https://github.com/antigenomics/arda
   cd arda

   arda rnaseq \
       --r1 tests/data/rnaseq_real/reads_1.fq.gz \
       --r2 tests/data/rnaseq_real/reads_2.fq.gz \
       -p DEMO -d out/ --threads 4

Progress goes to stderr, and the paths of the files written go to stdout, one per line:

.. code-block:: text

   [arda] map: reads_1.fq.gz + reads_2.fq.gz -> out/DEMO.airr.tsv | 4 threads, chunk 400000, min_score 75
   [arda] map: 453/1,320 reads mapped (34.32 %) in 1.9 s (681 reads/s), peak 352 MB;
          loci={'IGH': 104, 'IGK': 107, 'IGL': 104, 'TRB': 88, 'TRA': 45, 'TRG': 5}
   [arda] assemble: 10/25 complete contigs from 58 seeds; rescued 6 reads
   [arda] correct: 49 -> 49 clonotypes (0 collapsed) over 57 reads
   [arda] stats: 645 rows -> out/DEMO.stats.tsv

.. note::

   ``arda rnaseq`` is for whole-transcriptome libraries. For a targeted RepSeq or 5'RACE library
   run ``arda amplicon`` with the same arguments, and for single cell run ``arda cells``. The
   mode name carries that library's speed configuration — see :doc:`usage`.

3. Read the clonotype table
---------------------------

``out/DEMO.clones.tsv`` is one row per clonotype. Sorted by abundance, its head looks like this
(the first seven of its 21 columns):

.. list-table::
   :header-rows: 1
   :widths: 28 10 18 20 10 7 7

   * - ``junction_aa``
     - ``locus``
     - ``v_call``
     - ``j_call``
     - ``c_call``
     - ``duplicate_count``
     - ``consensus_count``
   * - ``CMQALQTFTF``
     - IGK
     - ``IGKV2D-28*02``
     - ``IGKJ4*01``
     - ``IGKC``
     - 3
     - 2
   * - ``CMQATQFPPLTF``
     - IGK
     - ``IGKV2-24*01``
     - ``IGKJ5*01``
     - ``IGKC``
     - 3
     - 2
   * - ``CARQDMIMFGGIRGYYGMDVW``
     - IGH
     - ``IGHV4-39*08``
     - ``IGHJ6*04``
     - ``IGH``
     - 2
     - 1
   * - ``CASSLRGSYEQYF``
     - TRB
     - ``TRBV7-8*01``
     - ``TRBJ2-7*01``
     - ``TRB``
     - 2
     - 1
   * - ``CQVWDSSSDHVVF``
     - IGL
     - ``IGLV3-21*04``
     - ``IGLJ2*01,IGLJ3*01``
     - ``IGLC2``
     - 2
     - 2

Three things are worth reading off it immediately.

``junction_aa`` **includes both anchors.** It runs from Cys104 to Phe/Trp118 inclusive, which is
why every string above opens with ``C``. IMGT CDR3 excludes them and is two residues shorter;
arda writes that separately as ``cdr3_aa``. See :term:`junction`.

``j_call`` **can name more than one gene.** ``IGLJ2*01,IGLJ3*01`` is a tie list: over the span
this read aligned, the two germlines are indistinguishable. Naming one would be a claim the data
does not support.

``duplicate_count`` **and** ``consensus_count`` **answer different questions.** The first counts
reads spanning the junction — abundance. The second counts distinct fragments — molecules. Pick
the one your analysis needs.

The full column list for every file is in :doc:`outputs`.

4. Check the run
----------------

``out/DEMO.arda.json`` is the run report: what was read, what mapped, and what it cost.

.. code-block:: bash

   jq '.map | {total_reads, mapped_reads, mapped_fraction, peak_rss_mb, unmapped}' out/DEMO.arda.json

.. code-block:: json

   {
     "total_reads": 1320,
     "mapped_reads": 453,
     "mapped_fraction": 0.3431818181818182,
     "peak_rss_mb": 352.1,
     "unmapped": {
       "prefilter_rejected": 788,
       "no_hit": 8,
       "constant_only": 54,
       "below_min_score": 17,
       "accounted": 1320
     }
   }

The ``unmapped`` breakdown accounts for every read that did not produce a row, so an unexpectedly
low mapped fraction has a cause rather than a shrug. This fixture is receptor-enriched; a genuine
bulk RNA-seq library maps 0.02–3 % of its reads, and that is normal rather than a failure.

``out/DEMO.stats.tsv`` is the same run as QC in long format — four columns, one value per cell:

.. code-block:: bash

   awk -F'\t' '$1=="chain" && $2=="IGH"' out/DEMO.stats.tsv

Per-sample QC, cohort roll-up and the self-contained HTML dashboard are covered in :doc:`qc`.

5. Analyse it
-------------

arda writes TSV with quoting **off**. Read it with quoting off too, or a blank field arrives as
the two-character string ``""``:

.. code-block:: python

   import polars as pl

   clones = pl.read_csv("out/DEMO.clones.tsv", separator="\t", quote_char=None)

   top = (
       clones
       .filter(pl.col("locus") == "TRB")
       .with_columns(
           (pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias("frequency")
       )
       .sort("duplicate_count", descending=True)
       .select("junction_aa", "v_call", "j_call", "duplicate_count", "frequency")
   )

More of these — clonal fractions, gene usage, repertoire overlap, gating a run on its own report,
cohort QC — are in :doc:`examples`.

Annotating sequences you already have
-------------------------------------

If your input is assembled sequences rather than reads, skip the pipeline and annotate directly.
FASTA or FASTQ, nucleotide or amino acid, all loci at once:

.. code-block:: bash

   arda annotate -i examples/example.fasta -o example.airr.tsv --organism human

Or from Python:

.. code-block:: python

   import arda

   records = arda.annotate_sequences(
       ["GACGTGCAG...", ("clone7", "CAGGTG...")],   # strings, or (id, sequence) pairs
       seqtype="nt",
       organism="human",
   )

Each record is a dict of AIRR fields. The TSV form passes ``airr.schema`` validation.

Next steps
----------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Question
     - Page
   * - Which mode does my library need?
     - :doc:`usage`
   * - What does this column mean?
     - :doc:`outputs`, :doc:`glossary`
   * - What does this command do?
     - :doc:`cli`
   * - Is this sample usable?
     - :doc:`qc`
   * - How do I run a whole cohort?
     - :doc:`samples`, :doc:`cluster`, :doc:`pipeline_integration`
   * - How is the annotation actually computed?
     - :doc:`how_it_works`
