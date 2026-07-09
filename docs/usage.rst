Usage
=====

Command line
------------

.. code-block:: bash

   arda info
   arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
   arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa

The output is a **spec-valid AIRR Rearrangement** TSV (it passes ``airr.schema``
validation) with 1-based, closed region coordinates (``fwr1_start``/``fwr1_end`` …
``cdr3_start``/``cdr3_end``), region nucleotide and amino-acid sequences,
``v_call``/``d_call``/``d2_call``/``j_call``, the constant-region ``c_call``/``c_class``
(isotype), per-segment CIGARs (``v_cigar``/``d_cigar``/``j_cigar``/``c_cigar``) and the
matching V/J/D germline coordinates, ``sequence_alignment`` / ``germline_alignment``,
``v_identity``, ``stop_codon``, ``vj_in_frame``, ``junction``, and ``productive``.
On a score tie ``d_call``/``d2_call`` are comma-separated allele ambiguity lists, and
``d2_call`` is the second (3′) segment of a D-D fusion — called in every D locus (IGH,
TRB, TRD). The ``sequence`` field holds the read **as submitted**; ``rev_comp`` = ``T``
signals that the other output fields describe its reverse complement (per the AIRR spec).

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

Bulk RNA-seq mode
-----------------

``arda rnaseq`` extracts the receptor repertoire from bulk RNA-seq, where only
1–5% of reads are receptor-derived:

.. code-block:: bash

   arda rnaseq map --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --report run.json
   arda rnaseq correct -i mapped.airr.tsv -o clones.tsv

``map`` streams paired FASTQ and writes only the reads that map to a receptor
scaffold. The reference includes ``J + C`` constant-region scaffolds, so a read
spanning the J→C splice still maps and carries a ``c_call`` (CH1 exon) and a
``c_class`` isotype (``IGHG``/``IGHM``/``IGHA`` … — the class, never the subclass).
With ``--reconstruct``, overlapping mates are merged into one fragment; a mismatch in
the overlap is resolved base-by-base in favour of the higher-Phred call (the default,
un-merged path is unchanged). ``correct`` collapses CDR3 sequencing errors into
clonotypes by a parent:child count ratio; a clonotype is keyed by
``(locus, v_call, j_call, junction)`` and its ``count`` is the number of distinct
fragments (paired mates counted once). ``arda igblast -i reads.fastq -o truth.airr.tsv``
runs IgBLAST across all loci as a gold-standard reference for benchmarking.

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
