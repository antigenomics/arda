API reference
=============

arda is usable as a library as well as a CLI. The single entry point most callers need is
:func:`arda.annotate_sequences`, which takes a batch of sequences and returns AIRR record dicts. It invokes MMseqs2 and
manages temporary search files internally:

.. code-block:: python

   import arda

   records = arda.annotate_sequences(
       [("q1", "TGTGCCAGCAGCTTAGCGGGAGGGAACACCGGGGAGCTGTTTTTTGGA")],
       seqtype="nt",          # "nt" or "aa"
       organism="human",      # human | mouse | rat | rabbit | rhesus_monkey
       threads=2,             # 0 (default) uses all CPUs allocated to this process/job
       map_d=True,            # D segments for VDJ loci; works on amino-acid input too
   )

Every input produces a record. Unmapped inputs retain ``sequence_id`` and ``sequence`` with
empty annotation fields. Submit one batch instead of calling the mapper once per sequence. The batch is held in memory; use the streaming CLI for a whole FASTQ.

The fields each record carries are documented in :doc:`outputs`. The modules below are the rest of
the public surface, grouped as the package is.

Library entry point
-------------------

.. automodule:: arda.adapter
   :members:
   :undoc-members:
   :show-inheritance:

Reference provenance
--------------------

``arda.reference_index_identity(organism, seqtype="aa")`` resolves and validates
an index without annotating queries. Its JSON-compatible result contains the
resolved ``source_fasta`` and ``database_path``, ``source_fasta_sha256``, and a
``files_sha256`` mapping from MMseqs payload basenames to SHA256 digests. The
mapping includes split sequence/header files and indexes. ``schema_version``,
``dbtype``, ``mmseqs_version``, ``organism`` and ``seqtype`` describe the identity.

Only ARDA's own source module triggers automatic source-root discovery; explicit
``ARDA_HOME`` remains authoritative. Index reuse checks content, not timestamps.
Legacy indexes remain untouched. New generations publish atomically and remain
immutable, so concurrent readers retain coherent sequence/header/index files.
Capture the API result before and after a batch and require equality to detect a
reference change during that batch.

Runtime annotation
------------------

.. automodule:: arda.annotate.mapper
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.transfer
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.reference
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.io
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.cigar
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.dmap
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.contig
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.annotate.shortlist
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.prefilter
   :members:
   :undoc-members:
   :show-inheritance:

Junction markup and repair
--------------------------

Working from a bare ``(junction_aa, v_call, j_call, species)`` record — a VDJdb row, with no
read to align — rather than from a sequenced fragment.

.. automodule:: arda.cdr3fix
   :members:
   :undoc-members:
   :show-inheritance:

.. note::

   ``arda.dpost`` was **removed in 2.33.0**. A posterior over the D gene is a
   recombination-model question, not a germline-reference one, so it lives in vdjtools as
   ``vdjtools.model.annotate_junctions``
   (`arda#144 <https://github.com/antigenomics/arda/issues/144>`_). The module was ported, not
   rewritten: its answer is identical field-for-field, and it gained the batch entry point
   `arda#142 <https://github.com/antigenomics/arda/issues/142>`_ asked for. arda still ships the
   prior table it reads (``database/vdj/<org>/d_prior.tsv``) and still fits one
   (:func:`arda.scenarios.estimate`).

.. automodule:: arda.rnaseq.map
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.rnaseq.assemble
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.rnaseq.correct
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.rnaseq.pipeline
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.rnaseq._res
   :members:
   :undoc-members:
   :show-inheritance:

Single cell
-----------

.. automodule:: arda.cell
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.singlecell
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.partition
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.scplot
   :members:
   :undoc-members:
   :show-inheritance:

Run QC and logging
------------------

.. automodule:: arda.stats
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda._log
   :members:
   :undoc-members:
   :show-inheritance:

Cluster sharding
----------------

.. automodule:: arda.cluster
   :members:
   :undoc-members:
   :show-inheritance:

Reference build
---------------

.. automodule:: arda.refbuild.build
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.refbuild.combinations
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.refbuild.segments
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.refbuild.translate
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.refbuild.imgt
   :members:
   :undoc-members:
   :show-inheritance:

External tool wrappers
----------------------

.. automodule:: arda.mmseqs
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda.igblast
   :members:
   :undoc-members:
   :show-inheritance:

Cache layout and reference fetch
--------------------------------

Where arda keeps its reference and how a plain ``pip install`` acquires one.

.. automodule:: arda.paths
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: arda._database_fetch
   :members:
   :undoc-members:
   :show-inheritance:

Preliminary segment evidence
----------------------------

.. autofunction:: arda.alignment_candidates
