Run within a CPU and memory budget
==================================

Use this guide when running several samples on a workstation or scheduling a
sample on a cluster. For a first run, start with :doc:`quickstart`; for command
signatures, use :doc:`cli` and :doc:`api`.

Choose the mapping configuration
--------------------------------

Use ``arda rnaseq`` for bulk RNA-seq and ``arda amplicon`` for targeted RepSeq
or 5'RACE. Their presets target different read geometries; combining their speed
flags is not a general optimization. ``arda map`` exposes individual flags for
controlled comparisons. See :doc:`usage` and :doc:`benchmarks` for the recorded
measurements supporting each preset.

Keep both mates for paired libraries. A constant-anchored IGH amplicon can need
assembly to recover its junction even when nearly every read maps. Retain the
default assembly stage; mapped fraction alone does not measure junction recovery.

Assign CPUs once
----------------

``--threads 0`` selects the CPUs available to the process: the host count,
restricted by process affinity, visible Linux cgroup CPU quotas and
``SLURM_CPUS_PER_TASK`` when set. Fractional quotas select at least one CPU and
round down to avoid exceeding the quota with an extra worker. A positive value selects an explicit count within that allocation;
negative values and requests exceeding it raise an error. The Python batch API
accepts the same ``threads`` argument.

.. code-block:: bash

   arda rnaseq --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --threads 4

Streaming commands reserve one CPU for reading and use the rest for search and native
prefilter/segment batches in successive phases. With a budget of one, reads are parsed
synchronously. An embedding application whose existing table pool fills the
entire budget also reads synchronously, keeping native search within that same
budget. In-memory batches use the full selected budget. The CLI keeps its Polars
pool at one thread; MMseqs preparation and reduction phases also have explicit limits. Submit sequences in one batch; an outer worker pool calling a threaded
search per sequence duplicates work and can oversubscribe the machine.

The sample-sheet CLI processes samples sequentially. To process samples or read
groups concurrently, use the existing SLURM, Nextflow or Snakemake integration
and give each task its own allocation. See :doc:`samples`, :doc:`cluster` and
:doc:`pipeline_integration`. ARDA initializes its table pool serially when it owns the first initialization.
An existing Polars pool is preserved; mapping raises an actionable error if that
pool exceeds the requested CPU budget. Initialize the embedding process with an
appropriate ``POLARS_MAX_THREADS`` before importing Polars in that case.

Measure the stage that matters
------------------------------

Read ``SAMPLE.arda.json`` for stage wall time and resource reports. Reported RSS
is a high-water mark across the process lifetime and reaped children; it is not
a sampled simultaneous total of parent and child RSS. To attribute a high-water
mark to one stage, run ``map``, ``assemble`` and ``correct`` separately.

When comparing thread counts or code versions, use identical inputs, reference
identities, mode and output options. Record wall time, CPU time, peak simultaneous
process-tree RSS, live thread counts and repeated-run variation. Include initial
reference/index preparation in first-use timing, and label subsequent mapping
timings separately. Require equality of calls and sequence identifiers, rather
than only equal row counts, before accepting a speed improvement.

Use the small committed fixtures for functional checks. Measure throughput and
scaling on a representative library; a short fixture can be dominated by process
startup and does not establish a speedup for a full cohort. Heavy Aldan-3 checks
belong in scheduled compute jobs.
