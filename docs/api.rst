API reference
=============

arda is usable as a library as well as a CLI. The single entry point most callers need is
:func:`arda.annotate_sequences`, which takes sequences and returns AIRR record dicts with no
subprocess and no temporary files:

.. code-block:: python

   import arda

   records = arda.annotate_sequences(
       ["GACGTGCAG...", ("clone7", "CAGGTG...")],   # strings or (id, sequence) pairs
       seqtype="nt",          # "nt" or "aa"
       organism="human",      # human | mouse | rat | rabbit | rhesus_monkey
       map_d=True,            # D segments for VDJ loci; works on amino-acid input too
   )

The fields each record carries are documented in :doc:`outputs`. The modules below are the rest of
the public surface, grouped as the package is.

Library entry point
-------------------

.. automodule:: arda.adapter
   :members:
   :undoc-members:
   :show-inheritance:

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
