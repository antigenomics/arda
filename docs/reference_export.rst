Exporting the reference
=======================

``arda export-ref`` writes arda's reference — every germline scaffold or segment together with
its FR1–4 / CDR1–3 markup — to stdout or to a file. The reference is arda's most valuable
offline artifact: an in-frame V×J germline set carrying IgBLAST-quality region coordinates. Until
this command existed it was reachable only by hand-joining the build's TSVs against its FASTAs.

.. code-block:: bash

   arda export-ref --kind scaffolds --format tsv --locus TRB -o trb_scaffolds.tsv

.. note::

   **Coordinates are 1-based CLOSED** — the AIRR convention, and also GFF3's, so GFF3 output
   passes coordinates through unchanged. ``fwr1_start = 1`` is the first base and ``fwr1_end``
   is the last base *of* FWR1, not one past it. An empty region pair means the record does not
   carry that region at all, which is not the same as a zero-length one.

The three kinds
---------------

.. list-table::
   :header-rows: 1
   :widths: 16 16 68

   * - ``--kind``
     - records
     - what it is
   * - ``scaffolds``
     - 15,414
     - The V×J (and J+C) reference the mapper aligns reads against. One record per germline
       combination, with a padded, in-frame V→J junction region.
   * - ``segments``
     - 924
     - The collapsed per-allele V / J / C reference the two-pass search nominates from. 775 V,
       124 J, 25 C targets.
   * - ``anchors``
     - 1,240
     - The per-allele CDR3 anchor table: ``anchor_nt``, ``templated_aa``, ``germline_nt``,
       ``status``.

Counts are for human, ``--seqtype nt``. Mouse gives 19,010 scaffolds; ``--seqtype aa`` gives
15,069 for human.

The four formats
----------------

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - ``--format``
     - output
   * - ``tsv``
     - Sequence plus every region as its own start / end / seq column triple. A leading ``#``
       comment line records the seqtype and the coordinate convention.
   * - ``fasta``
     - One record per scaffold or segment, with locus and gene calls in the description line.
   * - ``gff3``
     - Regions as features on each sequence. GFF3 is 1-based closed like arda, so coordinates
       are unchanged.
   * - ``airr``
     - The same rows shaped as an AIRR Rearrangement, so a scaffold can be fed straight into
       anything that already reads arda's own output.

.. warning::

   ``--kind anchors`` supports ``--format tsv`` **only**. It is a per-allele table with no
   sequence of its own, so the other three shapes have nothing to describe; asking for one
   raises rather than writing an empty file.

Other options: ``--organism`` (default ``human``), ``--seqtype`` (``nt`` or ``aa``), ``--locus``
(comma-separated, e.g. ``TRB,IGH``; default all loci), and ``-o``/``--output`` (default stdout).
The ``[arda] exported N ... record(s)`` progress line goes to **stderr**, so piping stdout stays
clean.

Worked examples
---------------

Scaffolds as TSV
~~~~~~~~~~~~~~~~

.. code-block:: bash

   arda export-ref --kind scaffolds --format tsv --locus TRB

.. code-block:: text

   # arda reference export (nt); coordinates are 1-based CLOSED (AIRR/GFF3 convention), an empty region means the record does not carry it
   sequence_id  locus  v_call       j_call       c_call  productive  junction  junction_aa       sequence_length  ...
   TRB_0        TRB    TRBV5-1*01   TRBJ2-5*01           T                     CASSLXXXXQETQYF   343              ...

Then ``fwr1_start``/``fwr1_end``/``fwr1_seq`` … ``fwr4_start``/``fwr4_end``/``fwr4_seq``. For
``TRB_0`` the region coordinates are FWR1 1–78, CDR1 79–93, FWR2 94–144, CDR2 145–162, FWR3
163–273, CDR3 274–312, FWR4 313–342.

The ``X`` runs in ``junction_aa`` are the padded V→J span: a scaffold is a germline
*combination*, not an observed rearrangement, so the untemplated N/P region has no sequence to
carry.

Scaffolds as FASTA
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   arda export-ref --kind scaffolds --format fasta --locus TRB

.. code-block:: text

   >TRB_0 locus=TRB TRBV5-1*01|TRBJ2-5*01
   AAGGCTGGAGTCACTCAAACTCCAAGATATCTGATCAAAACGAGAGGACAGCAAGTGACA
   CTGAGCTGCTCCCCTATCTCTGGGCATAGGAGTGTATCCTGGTACCAACAGACCCCAGGA
   CAGGGCCTTCAGTTCCTCTTTGAATACTTCAGTGAGACACAGAGAAACAAAGGAAACTTC

The description line carries the locus and the ``V|J`` pair, so the FASTA is self-describing
without the markup table beside it.

Scaffolds as GFF3
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   arda export-ref --kind scaffolds --format gff3 --locus TRB

.. code-block:: text

   ##gff-version 3
   ##sequence-region TRB_0 1 343
   TRB_0  arda  region  1    78   .  +  .  ID=TRB_0:fwr1;Name=FWR1;locus=TRB
   TRB_0  arda  region  79   93   .  +  .  ID=TRB_0:cdr1;Name=CDR1;locus=TRB
   TRB_0  arda  region  94   144  .  +  .  ID=TRB_0:fwr2;Name=FWR2;locus=TRB
   TRB_0  arda  region  145  162  .  +  .  ID=TRB_0:cdr2;Name=CDR2;locus=TRB
   TRB_0  arda  region  163  273  .  +  .  ID=TRB_0:fwr3;Name=FWR3;locus=TRB
   TRB_0  arda  region  274  312  .  +  .  ID=TRB_0:cdr3;Name=CDR3;locus=TRB
   TRB_0  arda  region  313  342  .  +  .  ID=TRB_0:fwr4;Name=FWR4;locus=TRB

Pair it with the FASTA export of the same ``--kind``/``--locus`` and any genome browser or
``bedtools``-style toolchain can work against the germline directly.

Scaffolds as AIRR
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   arda export-ref --kind scaffolds --format airr --locus TRB -o trb_ref.airr.tsv

The rows carry the AIRR Rearrangement fields arda emits for real reads — ``sequence_id``,
``sequence``, ``locus``, ``v_call``/``d_call``/``d2_call``/``j_call``/``c_call``, ``c_class``,
``rev_comp``, ``productive``, ``stop_codon``, ``vj_in_frame``, ``v_identity``,
``sequence_alignment``, ``germline_alignment``, the per-segment CIGARs, and the ``mmseqs2_*``
alignment columns. A germline scaffold can therefore be fed into any consumer of arda's own
output — useful as a positive control, since a scaffold aligned against itself should annotate
perfectly.

Segments
~~~~~~~~

.. code-block:: bash

   arda export-ref --kind segments --format tsv --locus IGH

.. code-block:: text

   sequence_id     locus  v_call        j_call  c_call  ...  sequence_length
   V|IGHV1-18*01   IGH    IGHV1-18*01                   ...  296
   V|IGHV1-18*03   IGH    IGHV1-18*03                   ...  296

Segment ids are prefixed by kind: ``V|``, ``J|``, ``C|``. A V segment carries FWR1–FWR3 in full
plus the CDR3 *prefix* it templates — for ``V|IGHV1-18*01`` that is FWR1 1–75, CDR1 76–99, FWR2
100–150, CDR2 151–174, FWR3 175–288, CDR3 289–296. The eight CDR3 bases are the germline
contribution up to Cys104; everything 3′ of that is recombination and belongs to no segment.

Anchors
~~~~~~~

.. code-block:: bash

   arda export-ref --kind anchors --format tsv --locus TRA

.. code-block:: text

   locus  segment  allele        functionality  anchor_nt  partial_nt  templated_aa  germline_nt                       status     source
   TRA    V        TRAV1-1*01    F              261        2           CAVR          TGCGCTGTGAGAGA                    ok         ndm
   TRA    V        TRAV1-2*02    F              -1         0                                                           no_anchor  no_anchor
   TRA    J        TRAJ10*01     F              30         0           ILTGGGNKLTF   ATACTCACGGGAGGAGGAAACAAACTCACCTTT  ok         aux

``anchor_nt`` is the 0-based offset of the anchor codon within the germline segment, and
``partial_nt`` how many bases of a split codon precede it. ``status = no_anchor`` (with
``anchor_nt = -1``) marks an allele whose anchor could not be placed — it is reported rather
than silently dropped, because an unanchorable allele is a real gap in the reference vocabulary
and a consumer may need to know which genes it affects.

.. note::

   Read ``anchor_nt`` from this table rather than re-deriving the anchor with a ``[FW]GXG``
   motif search. A motif check is not an anchor: ``TRAJ35*01`` is a functional gene whose anchor
   codon decodes **Cys (TGC), not [FW]**, and a motif-based scan silently deletes it.

Reproducibility
---------------

The export is deterministic. Where IMGT ships two accessions under one allele name, the anchor
table can carry two rows for the same key with different ``templated_aa`` and ``germline_nt``.
These are resolved by an explicit rule — prefer ``status = ok``, then the row that templates
*more* junction, then the sequence itself as a final tie-break — and the choice is logged to
stderr rather than made silently:

.. code-block:: text

   cdr3_anchors.tsv has conflicting rows for V IGKV10-96*01 in mouse; keeping status=ok
   germline_nt=TGCCAACAGGGTAGTACGCTTCCTCC (IMGT ships two accessions under one allele name)

Three mouse alleles are affected (``IGKV10-96*01``, ``IGLV2*01``, ``IGLV3*01``). The loader used
to be last-wins, so the winner depended on row order and two builds of the same reference could
disagree about which germline the Cys104 gate scored against.
