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
* ``d_prior.tsv`` — generative-model summaries read by :func:`arda.scenarios.load_prior_table`,
  :func:`arda.hmm.model_for`, the SHM-aware scorer
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
J coding sequence in frame with V. Auxiliary coordinates also cover
noncanonical J anchors such as the Cys118 of TRAJ35*01.

Opt-in TRDV–TRAJ scaffolds
--------------------------

``arda build-db --allow-chimeras`` enables the existing optional TRDV × TRAJ
cross-product. It does not change the default pairing policy. Build into an
isolated reference checkout when evaluating it; do not replace the shipped
reference merely to test a candidate.

The flag's legacy name does not imply that every such pairing is chimeric.
`IMGT documents TRDV1*01 rearranged to TRAJ54*01
<https://www.imgt.org/IMGTrepertoire/Regulation/cDNA/human/TRD/TRDV/Hu_TRDVcdna.html>`_
and another TRDV1–TRAJ example. A TRDV V-gene name alone does not establish a
TRD chain; J/D/C evidence determines the locus. Ordinary TRG/TRD annotation
is already supported by the default reference. This option is not a separate
γδ mode, and the broader cross-product remains a reference-scope choice.

IgBLAST can identify the exact TRDV germline while omitting J markup. For these
constructed scaffolds, the builder can recover J coordinates from the known
``V + N-padding + J`` sequence and the same germline anchor authority used by
``cdr3_anchors.tsv``. Auxiliary coordinates take precedence, including the
noncanonical Cys118 of TRAJ35*01. The existing motif fallback applies only to
alleles without auxiliary metadata. J padding follows the anchor frame when
the auxiliary coding-frame field disagrees with its own anchor.

The opt-in build fails if any expected TRDV × TRAJ scaffold lacks complete
markup. It appends new targets while preserving existing target IDs, order
and padding, so vocabulary expansion does not reshuffle equal-score legacy
alignments. Check every region's sequence/coordinate round trip before packaging a
candidate. Complete scaffold coverage establishes reference construction;
measure recovered junctions and changed calls on real reads separately.
