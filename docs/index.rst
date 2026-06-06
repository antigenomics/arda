arda — Antigen Receptor Domain Annotation
=========================================

Fast FR/CDR region annotation for TCR and BCR nucleotide and amino-acid
sequences. ``arda`` builds a pre-aligned IgBLAST reference database once, then
maps queries with MMseqs2 and transfers the region markup through the alignment
in a small C++ hot path — producing AIRR-formatted output that matches IgBLAST.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   reference_build
   api

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
