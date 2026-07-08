"""RNA-seq mode: recall-first filtering + mapping of receptor reads from bulk RNA-seq.

Bulk RNA-seq is mostly (95-99 %) non-receptor. This package reuses the streaming
annotator (``arda.annotate.mapper``) — whose MMseqs2 k-mer prefilter already rejects
non-receptor reads cheaply — but, unlike ``arda annotate``, keeps only the reads that
map to a V/J scaffold, keyed by read id.

Stages (see the ``rna-seq`` plan):

* ``map``     — filter + map (this module: :mod:`arda.rnaseq.map`).
* ``correct`` — CDR3 error correction (:mod:`arda.rnaseq.correct`).
* assembly    — deferred (see ``ROADMAP.md``).
"""

from __future__ import annotations

from .map import map_rnaseq, read_pairs, RnaseqReport

__all__ = ["map_rnaseq", "read_pairs", "RnaseqReport"]
