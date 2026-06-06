Usage
=====

Command line
------------

.. code-block:: bash

   arda info
   arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
   arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa

The output is an AIRR rearrangement TSV with 1-based, closed region coordinates
(``fwr1_start``/``fwr1_end`` … ``cdr3_start``/``cdr3_end``), region nucleotide and
amino-acid sequences, ``v_call``/``j_call``, ``junction``, and ``productive``.

Python library
--------------

.. code-block:: python

   import arda

   records = arda.annotate_sequences(
       ["GACGTGCAG...", ("clone7", "CAGGTG...")],
       seqtype="nt",
       organism="human",
   )

Each record is a dict keyed by the AIRR fields above.

Scaling
-------

MMseqs2 runs multi-threaded (``--threads``); inputs may be FASTA or FASTQ, plain
or gzipped. For cluster runs, shard the input across SLURM array tasks and
concatenate the per-shard AIRR TSVs.

Supported organisms
-------------------

* **human, mouse** — full IG and TR loci.
* **rat, rabbit, rhesus_monkey** — IG only (IgBLAST ships no TR internal
  annotation for these organisms).
