Glossary
========

Terms used throughout these pages, in the sense arda uses them.

.. glossary::
   :sorted:

   AIRR Rearrangement
      The `AIRR Community <https://docs.airr-community.org/>`_ standard schema for one annotated
      receptor sequence: gene calls, region coordinates, the junction and its translation, with
      defined column names and value formats. Every per-sequence TSV arda writes is a valid AIRR
      Rearrangement file. See :doc:`outputs`.

   allele
      The most specific level of IMGT gene nomenclature, written with a ``*`` suffix:
      ``TRBV19*01``. Separating a gene's alleles usually needs substantially more read coverage
      than naming the gene — see :doc:`genotype`.

   amplicon
      A targeted library in which nearly every read is receptor-derived and reads span V into J,
      produced by multiplex V primers or 5'RACE. Sometimes called RepSeq. Handled by
      :ref:`arda amplicon <cli-amplicon>`.

   anchor
      The conserved residue that bounds the junction: Cys104 on the V side, Phe or Trp 118 on the
      J side, in IMGT numbering. arda stores the anchor position and templated residues per
      allele in ``cdr3_anchors.tsv`` and refuses to emit a junction whose V anchor the read never
      reached.

   bulk RNA-seq
      A whole-transcriptome library in which 0.02–3 % of reads are receptor-derived and a read
      lands anywhere in a transcript rather than spanning V into J. Handled by
      :ref:`arda rnaseq <cli-rnaseq>`.

   CDR3
      Complementarity-determining region 3, in the IMGT definition: the junction **excluding**
      both conserved anchors. See :term:`junction` — the two differ by two residues and
      conflating them corrupts every downstream coordinate.

   chimera
      A clonotype whose V and J calls come from loci that cannot rearrange together. arda counts
      them in :doc:`qc` and flags them; it never removes them.

   clonotype
      A group of reads collapsed onto one rearrangement. The grouping key is set by
      ``--clonotype-key``: ``full`` (locus, V, J, junction) by default, or ``junction`` alone.
      One row of ``<prefix>.clones.tsv``.

   consensus_count
      AIRR column: the number of distinct fragments supporting a clonotype. Use it when you need
      molecules. Compare :term:`duplicate_count`.

   duplicate_count
      AIRR column: the number of reads encompassing the junction of a clonotype. arda's
      expression estimate. Use it when you need abundance. Invariant across every
      ``--ec-mode`` — error correction moves reads onto a parent, it never discards them.

   ec-mode
      The denoising preset given to ``arda correct`` (``--ec-mode``). Chosen from what the answer
      will be used for, not from the library type. See :doc:`error_correction`.

   fast_fraction
      Reported in ``<prefix>.arda.json``: the fraction of reads that hit **both** a V and a J
      segment in the two-pass search. It predicts whether the amplicon speed configuration will
      pay, better than the library's name does. Roughly .85 on a primer-anchored TCR amplicon,
      .50 on IGH RepSeq, .05 on bulk.

   IMGT
      The `international ImMunoGeneTics information system <https://www.imgt.org/>`_: the source
      of arda's germline sequences and of the region and numbering conventions it follows.

   junction
      The nucleotide (or amino acid) stretch from Cys104 through Phe/Trp118 **inclusive of both
      anchors**. This is what AIRR calls ``junction`` and what MiXCR labels ``CDR3``. It is two
      residues longer than IMGT :term:`CDR3`, which arda writes separately as ``cdr3``.

   J+C scaffold
      A reference entry consisting of a J germline with the CH1 constant exon spliced onto it,
      and no V. It gives reads that span the J→C splice — which end in the constant region and
      have no V to anchor on — somewhere to map, and is what lets arda report a ``c_call`` for
      them. 345 of them ship for human.

   locus
      The receptor chain a sequence belongs to: ``TRA``, ``TRB``, ``TRG``, ``TRD``, ``IGH``,
      ``IGK``, ``IGL``. arda searches all of them at once and reports the one that won.

   mode
      One of ``arda rnaseq``, ``arda amplicon``, ``arda cells``. The mode name carries the speed
      configuration for its library type, because the configurations do not compose. See
      :doc:`usage`.

   read group
      One ``(R1, R2)`` FASTQ pair belonging to a sample. A sample delivered as several Illumina
      lanes is one repertoire made of several read groups; arda maps each and concatenates before
      clonotype calling. See :doc:`samples`.

   regime
      The library class a speed configuration is tuned for — bulk or amplicon. Named once in
      ``arda.cluster.regime_flags`` so the CLI, SLURM, Nextflow and Snakemake cannot disagree.

   scaffold
      A reference entry made of one V allele and one J allele joined in frame, with a short
      frame-neutral N spacer where D would sit, carrying IgBLAST-derived FR1–FR4 / CDR1–CDR3
      coordinates. Queries are mapped to scaffolds and the markup is projected back onto them.
      15,069 V·J scaffolds ship for human.

   segment
      A single germline allele (V, J or C) in the collapsed 924-target reference that the
      two-pass search screens against, as opposed to a full V·J :term:`scaffold`.

   SHM
      Somatic hypermutation: the AID-driven mutation of rearranged IG loci in germinal centres.
      arda reports it as ``v_identity``, ``v_mutations`` and ``j_mutations``, scoped outside the
      junction by default. See :doc:`shm`.

   stage
      One of the three steps a mode runs: ``map`` (Stage 1, per read), ``correct`` (Stage 2,
      global), ``assemble`` (Stage 3, global). Each is also a standalone command, which is what
      makes sharding possible.

   tie list
      A comma-joined set of gene or allele names that the alignment cannot separate, written into
      a single ``v_call`` / ``j_call`` / ``d_call`` field. Standard AIRR practice; scoring a tie
      list as a miss is a scoring artifact rather than a calling error.

   UMI
      Unique molecular identifier: a random barcode attached to a cDNA molecule before
      amplification, used to collapse PCR duplicates into one consensus. arda does not collapse
      UMIs; ``arda cells`` consumes the consensus an upstream tool produced.
