arda
====

.. raw:: html

   <div class="proj-intro">
     <div>
       <p class="proj-intro__eyebrow">ANTIGEN RECEPTOR DOMAIN ANNOTATION</p>
       <p class="proj-intro__lead">arda annotates T&nbsp;and B&nbsp;cell receptor sequences.
       Give it sequencing reads &mdash; bulk RNA-seq, targeted amplicon or single cell &mdash;
       and it returns AIRR Rearrangement records with V/D/J/C gene calls, FR1&ndash;FR4 and
       CDR1&ndash;CDR3 coordinates and the junction, plus a clonotype table. One command,
       no workflow engine, no container.</p>
       <p class="proj-intro__links">
         <a href="quickstart.html">Quickstart</a>
         <span>&middot;</span>
         <a href="installation.html">Install</a>
         <span>&middot;</span>
         <a href="cli.html">Commands</a>
         <span>&middot;</span>
         <a href="outputs.html">Output files</a>
       </p>
     </div>
   </div>

.. code-block:: bash

   pip install arda-mapper
   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="quickstart.html">
       <h3>Quickstart</h3>
       <p>Install, run a library, read the table that comes out. Ten minutes.</p>
     </a>
     <a class="proj-card" href="usage.html">
       <h3>Choosing a mode</h3>
       <p>Bulk RNA-seq, amplicon or single cell &mdash; and what each one writes.</p>
     </a>
     <a class="proj-card" href="outputs.html">
       <h3>Output reference</h3>
       <p>Every file a run produces, and every column in it.</p>
     </a>
     <a class="proj-card" href="examples.html">
       <h3>Recipes</h3>
       <p>Copy-paste shell and polars snippets for the analyses that follow a run.</p>
     </a>
   </div>

What it does
------------

**Annotation.** A sequence in, an AIRR Rearrangement record out: ``v_call``, ``d_call``,
``j_call``, ``c_call``, the four framework regions, the three CDRs, the junction and its
translation, all in 1-based closed query coordinates. Nucleotide or amino acid input, FASTA or
FASTQ, all loci searched at once.

**Repertoires.** ``arda rnaseq``, ``arda amplicon`` and ``arda cells`` take raw reads to a
clonotype table in one call, running mapping, contig assembly and error correction for you.

**Records with no read behind them.** A CDR3 amino acid plus a V and a J call — a VDJdb row —
can be marked up, its junction repaired against the germline anchors, and its D gene inferred
from junction length.

Why arda
--------

.. raw:: html

   <div class="proj-feature-grid">
     <div class="proj-feature">
       <h3>IgBLAST-quality regions, without IgBLAST at runtime</h3>
       <p>Region concordance with IgBLAST is 98&ndash;99.7&nbsp;% on ~7,300 real GenBank mRNA
       records across five organisms. IgBLAST itself is needed only to rebuild the reference,
       never to annotate.</p>
     </div>
     <div class="proj-feature">
       <h3>Cheap enough to run per sample on a shared node</h3>
       <p>Mapping holds a flat 300&ndash;650&nbsp;MB of RAM at any read depth. On a 100,000-read
       TRA amplicon arda uses 3.2&times; less peak memory than MiXCR and 1.4&times; less CPU.</p>
     </div>
     <div class="proj-feature">
       <h3>Precise on bulk RNA-seq</h3>
       <p>Across 16 bulk datasets and 5,273 real fragments, recall 0.986 [.982&ndash;.989] and a
       precision lower bound of 0.889 [.881&ndash;.897] &mdash; a precision bound above every
       compared tool's upper bound.</p>
     </div>
     <div class="proj-feature">
       <h3>It declines rather than guessing</h3>
       <p>A D call carries the E-value it was accepted on. A junction whose Cys104 anchor is not
       in the read is not emitted. A germline with no derivable anchor is flagged, not
       invented.</p>
     </div>
     <div class="proj-feature">
       <h3>Runs where your data already is</h3>
       <p>A plain CLI over named files: the same sample sheet drives the CLI, a SLURM array, the
       Nextflow module and the Snakemake workflow. A sharded run is byte-identical to a
       single-node one.</p>
     </div>
     <div class="proj-feature">
       <h3>Embeddable</h3>
       <p><code>import arda; arda.annotate_sequences(...)</code> returns AIRR record dicts. No
       subprocess, no temporary files, no server.</p>
     </div>
   </div>

Every performance and accuracy figure quoted in these pages is measured, and the method behind
each one is recorded on the :doc:`benchmarks` page.

Where to go next
----------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - If you want to…
     - Read
   * - run a library end to end for the first time
     - :doc:`quickstart`
   * - decide which mode your library needs
     - :doc:`usage`
   * - understand a column in the output
     - :doc:`outputs`
   * - look up a command or a flag
     - :doc:`cli`
   * - analyse a clonotype table
     - :doc:`examples`
   * - scale to a cohort or a cluster
     - :doc:`cluster`, :doc:`pipeline_integration`
   * - know how the annotation is computed
     - :doc:`how_it_works`

.. toctree::
   :hidden:
   :caption: Get started
   :maxdepth: 2

   quickstart
   installation
   usage

.. toctree::
   :hidden:
   :caption: Guides
   :maxdepth: 2

   examples
   use_cases
   singlecell
   qc
   samples
   cluster
   pipeline_integration

.. toctree::
   :hidden:
   :caption: Reference
   :maxdepth: 2

   cli
   outputs
   glossary
   benchmarks
   reference_export
   reference_build
   api

.. toctree::
   :hidden:
   :caption: Background
   :maxdepth: 2

   how_it_works
   productivity
   error_correction
   shm
   d_segments
   genotype
   scenarios

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
