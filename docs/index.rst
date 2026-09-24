arda
====

.. raw:: html

   <div class="proj-intro">
     <div>
       <p class="proj-intro__eyebrow">ANTIGEN RECEPTOR DOMAIN ANNOTATION</p>
       <p class="proj-intro__lead">arda does the expensive IgBLAST work <strong>once, offline</strong>
       &mdash; a reference of every in-frame V&middot;J germline scaffold with FR1&ndash;4 /
       CDR1&ndash;3 markup &mdash; then maps reads to it with MMseqs2 and transfers the markup
       through the alignment in C++. Out comes spec-valid AIRR Rearrangement, from FASTQ or from
       assembled sequences, for TCR and BCR.</p>
       <p class="proj-intro__links">
         <a href="installation.html">Install</a>
         <span>&middot;</span>
         <a href="usage.html">The three modes</a>
         <span>&middot;</span>
         <a href="examples.html">Recipes</a>
         <span>&middot;</span>
         <a href="api.html">API</a>
       </p>
     </div>
   </div>

   <div class="proj-card-grid">
     <a class="proj-card" href="installation.html">
       <h3>Installation</h3>
       <p>pip, conda or a source checkout; the reference fetches itself on first use.</p>
     </a>
     <a class="proj-card" href="usage.html">
       <h3>Bulk RNA-seq &amp; amplicon</h3>
       <p><code>arda rnaseq</code> and <code>arda amplicon</code>: which one your library wants,
       and what each writes.</p>
     </a>
     <a class="proj-card" href="singlecell.html">
       <h3>Single cell</h3>
       <p><code>arda cells</code>: reference-free per-cell contigs from UMI consensus, chain
       pairing, doublets and QC.</p>
     </a>
     <a class="proj-card" href="examples.html">
       <h3>Recipes</h3>
       <p>Copy-paste CLI and polars snippets for the analyses that follow a run.</p>
     </a>
   </div>

   <div class="proj-feature-grid">
     <div class="proj-feature">
       <h3>One command, FASTQ to clonotypes</h3>
       <p><code>arda rnaseq</code> and <code>arda amplicon</code> run
       map &rarr; assemble &rarr; correct in one call and write an AIRR table, a clonotype table,
       a run report and a QC table. The mode <em>name</em> picks the speed configuration, because
       the two speed levers do not compose and each is a loss in the other's regime.</p>
     </div>
     <div class="proj-feature">
       <h3>Precision, on real bulk RNA-seq</h3>
       <p>Over 16 datasets where every tool ran (5,273 real fragments, Wilson 95&nbsp;% CIs):
       recall <strong>0.986</strong> [.982&ndash;.989], precision lower bound
       <strong>0.889</strong> [.881&ndash;.897]. Recall ties TRUST4; precision does not &mdash;
       arda's lower bound sits above every competitor's upper bound.</p>
     </div>
     <div class="proj-feature">
       <h3>Cheap, not just fast</h3>
       <p>On a TRA amplicon of 100,000 reads arda is <strong>1.10&times;</strong> MiXCR's wall
       clock at <strong>3.6&times;</strong> less CPU and <strong>4.8&times;</strong> less RSS
       (631 MB against 3,027 MB). <code>map</code> is flat at 300&ndash;650 MB at any depth, which
       is what lets many samples share a node.</p>
     </div>
     <div class="proj-feature">
       <h3>Error correction you choose by question</h3>
       <p>Denoising is set by what you will do with the answer, not by the library: abundance
       collapse, quality-directed rescue at 12 subs / 50&times; for a deep targeted library,
       the same rescue kept narrow for bulk RNA-seq, or nothing at all.</p>
     </div>
     <div class="proj-feature">
       <h3>Single cell, with the doublets called</h3>
       <p>From one UMI consensus per molecule, <code>arda cells</code> assembles each cell's
       contigs with <em>no germline reference</em>, then annotates, pairs the chains and flags
       multiplets. Against Cell Ranger on <code>sc5p_v2_hs_PBMC_1k</code> VDJ-T,
       <strong>98.2&nbsp;%</strong> of its 943 CDR3s appear verbatim inside one of arda's
       contigs.</p>
     </div>
     <div class="proj-feature">
       <h3>It runs where your data is</h3>
       <p>The same sheet drives the CLI, a SLURM array, the Nextflow module and the Snakemake
       workflow, and all four ask <code>arda.cluster.regime_flags</code> for a regime's flags
       rather than restating them. A sample delivered as several lane FASTQs is handled in one
       call &mdash; see <a href="samples.html">samples split across files</a>.</p>
     </div>
   </div>

Every speed and accuracy figure on this site is measured, and the measurements live outside this
repository so they can be re-run: see the ``README`` benchmark tables and
``~/vcs/projects/2026-arda-benchmark``.

.. toctree::
   :hidden:
   :caption: Start here
   :maxdepth: 2

   installation
   usage
   singlecell
   examples

.. toctree::
   :hidden:
   :caption: Analysis
   :maxdepth: 2

   use_cases
   productivity
   error_correction
   shm
   d_segments

.. toctree::
   :hidden:
   :caption: At scale
   :maxdepth: 2

   cluster
   samples
   pipeline_integration
   qc

.. toctree::
   :hidden:
   :caption: Reference
   :maxdepth: 2

   reference_build
   reference_export
   api

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
