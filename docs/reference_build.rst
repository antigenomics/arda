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
* ``d_germlines.fasta`` — per-locus D germlines (VDJ loci only), for runtime D mapping.
* ``cdr3_anchors.tsv`` — **per allele**, the conserved Cys104 (V) or Phe/Trp118 (J)
  position, the residues that allele templates into the junction, and a ``status``
  (``ok`` / ``truncated`` / ``no_anchor`` — flagged, never guessed). This is what lets
  :mod:`arda.cdr3fix` mark up a junction with no read behind it, and what pins the
  V..J interior for D mapping.
* ``d_prior.tsv`` — generative-model summaries used by :mod:`arda.dpost`
  (``insVD``/``insDJ`` insert lengths, surviving-D length, ``P(D | J)``). Derived, not
  measured; shipped only for the (organism, locus) pairs with a published model.
  See ``SOURCES.md`` and ``scripts/build_d_priors.py``.
* ``loci_manifest.tsv`` — one row per defined locus: V·J and J+C scaffold
  counts, D-germline count, unreachable-D count, and ``ok`` / ``EMPTY`` status.
* ``build.log`` — per-locus counts and dropped/incomplete summaries.

The build warns at the end on any ``EMPTY`` locus (no scaffolds — reads from it
fall through unannotated) and on any D germlines shipped for a locus with no
scaffolds (unreachable, since runtime D lookup is keyed by a hit scaffold's
locus). IMGT carries no TCR reference for rat, rabbit, or rhesus, so their TR
loci (TRA, TRB, TRG, TRD) build ``EMPTY`` — now recorded in the manifest and
warned, not silently skipped. Only human and mouse ship a full TCR reference.

Reading frames
--------------

V germline is normalized to its coding frame by stop-free frame detection (the
IMGT "codon start" header field is unreliable for 5'-partial alleles). The J
coding frame comes from the IgBLAST auxiliary file; junction N-padding keeps the
J frame aligned to V so the conserved FR4 ``[FW]GXG`` motif translates correctly.
