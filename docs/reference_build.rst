Reference database build
========================

The reference build is an offline, reproducible process (``arda build-db``). It
is only needed to regenerate ``database/vdj/<organism>/``; normal annotation uses
the committed references.

Pipeline
--------

#. **Download** the IMGT/V-QUEST germline reference directory.
#. **Ungap** each gene file with IgBLAST's ``edit_imgt_file.pl``.
#. **Enumerate** deduplicated in-frame V·J scaffolds. The D segment is *not*
   enumerated — it only affects the CDR3 interior, which is query-specific at
   runtime — but VDJ loci get a short frame-neutral N spacer where D would sit.
#. **Annotate** the scaffolds with ``igblastn -outfmt 19`` (AIRR) and extract the
   FR/CDR coordinates with polars.
#. **Translate** each scaffold and derive protein markup.

Outputs (per organism)
----------------------

* ``alleles.fasta`` / ``alleles.aa.fasta`` — scaffold nucleotide / protein seqs.
* ``markup.tsv`` / ``markup.aa.tsv`` — region coordinates + sequences.
* ``combinations.tsv`` — scaffold → contributing (V, J) allele pairs.
* ``build.log`` — per-locus counts and dropped/incomplete summaries.

Reading frames
--------------

V germline is normalized to its coding frame by stop-free frame detection (the
IMGT "codon start" header field is unreliable for 5'-partial alleles). The J
coding frame comes from the IgBLAST auxiliary file; junction N-padding keeps the
J frame aligned to V so the conserved FR4 ``[FW]GXG`` motif translates correctly.
