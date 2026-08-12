API reference
=============

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

.. automodule:: arda.dpost
   :members:
   :undoc-members:
   :show-inheritance:

Bulk RNA-seq
------------

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
