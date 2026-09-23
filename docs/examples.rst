Recipes
=======

Runnable snippets for the things people actually do after (and around) a run. Every command and
every script on this page was executed against ``tests/data/rnaseq_real`` before it was written
down; the printed numbers are that fixture's, not illustrations.

:doc:`usage` is the reference for what each flag means, :doc:`samples` for multi-file samples,
:doc:`use_cases` for choosing a mode and a denoising preset.

.. contents::
   :local:
   :depth: 2

FASTQ to clonotypes
-------------------

One command per library type. The mode **name** picks the speed configuration — the two speed
levers do not compose, and each is a loss in the other's regime, so do not hand-assemble them.

.. code-block:: bash

   # bulk RNA-seq / WTS: a few receptor reads in a lot of transcriptome
   arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz --out-prefix PT01 -d results/ --threads 16

   # targeted RepSeq / 5'RACE amplicon: nearly every read is a receptor
   arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz --out-prefix PT01 -d results/ --threads 16

   # 10x-style single cell: per-cell contigs, chain pairing, doublets
   arda cells    --r1 R1.fq.gz --r2 R2.fq.gz --out-prefix PT01 -d results/ --threads 16

Four files come out, and the paths are echoed on **stdout** (one per line) while progress goes to
stderr — so ``$(arda rnaseq ... | head -1)`` is a usable idiom in a shell script:

.. code-block:: text

   results/PT01.airr.tsv     one AIRR Rearrangement row per mapped read
   results/PT01.clones.tsv   one row per clonotype, with duplicate_count
   results/PT01.arda.json    the run report: what was read, what mapped, wall time, peak RSS
   results/PT01.stats.tsv    run QC in long format (scope / key / metric / value)

On the 1,320-read test fixture that is ``453/1,320 reads mapped (34.32 %)`` and 49 clonotypes.

A sample delivered as several FASTQ pairs
-----------------------------------------

Lanes and in-house chunks are **read groups** of one sample, and one sample gives one clonotype
table. Repeat ``--r1``/``--r2`` and name the sample with ``--id`` — the grouping is declared,
never guessed from filenames:

.. code-block:: bash

   arda rnaseq -d results/ --threads 16 \
       --r1 PT01_S1_L001_R1_001.fastq.gz --r2 PT01_S1_L001_R2_001.fastq.gz --id PT01 \
       --r1 PT01_S1_L002_R1_001.fastq.gz --r2 PT01_S1_L002_R2_001.fastq.gz --id PT01

Do **not** ``cat`` them first: arda maps each read group and concatenates after Stage 1, which is
byte-identical to the same reads in one file and skips a full copy of the data. :doc:`samples` has
the sheet form, the ordering rule and the one-worker-per-read-group recipe.

A cohort from one sheet
-----------------------

.. code-block:: bash

   cat > sheet.tsv <<'SHEET'
   sample	fastq_1	fastq_2
   PT01	PT01_S1_L001_R1_001.fastq.gz	PT01_S1_L001_R2_001.fastq.gz
   PT01	PT01_S1_L002_R1_001.fastq.gz	PT01_S1_L002_R2_001.fastq.gz
   PT02	PT02_S2_L001_R1_001.fastq.gz	PT02_S2_L001_R2_001.fastq.gz
   SHEET

   arda rnaseq --samples sheet.tsv -d results/ --threads 16      # one box, samples in turn
   arda cluster submit-samples --samples sheet.tsv --work-dir work/ -d results/ \
       --regime rnaseq --threads 16 --partition medium --submit  # SLURM, one task per read group

The columns are nf-core's, so an existing nf-core samplesheet works unmodified. Repeated
``sample`` values merge in row order.

Annotate sequences you already have
-----------------------------------

No reads, just sequences — assembled contigs, Sanger, a synthesised panel, someone else's
consensus. ``arda annotate`` streams FASTA/FASTQ to AIRR and is memory-bounded:

.. code-block:: bash

   arda annotate -i contigs.fasta -o contigs.airr.tsv --organism human
   cut -f1,3,4,7,55 contigs.airr.tsv | head -3

.. code-block:: text

   sequence_id	locus	v_call	j_call	junction_aa
   PZ235980.1	IGH	IGHV3-9*01	IGHJ6*03	CARDIGAGGFGDNFYFFYYMDVW
   PV083657.1	IGK	IGKV1-33*01,IGKV1D-33*01	IGKJ5*01	CQQYDSLPYTF

Bare ``(junction_aa, V, J)`` records — a VDJdb-style table, a published supplement — go through
``arda markup``, which marks up the junction and repairs an anchor the record lost:

.. code-block:: bash

   arda markup -i junctions.tsv -o fixed.tsv --id-col id
   # 7 records -> fixed.tsv  (5 repaired, 1 failed)

``FLVGPQGSSASKIIF`` comes back as ``CLVGPQGSSASKIIF`` (the Cys104 anchor restored) and
``CAIRDDKII`` as ``CAIRDDKIIF``. Add ``--vdjdb`` to read VDJdb's own column names.

From Python
-----------

The library entry point takes sequences and returns AIRR record dicts — no files, no subprocess:

.. code-block:: python

   from arda import annotate_sequences

   rows = annotate_sequences(
       [("q1", "TGTGCCAGCAGCTTAGCGGGAGGGAACACCGGGGAGCTGTTTTTTGGA")], seqtype="nt"
   )
   print({k: rows[0][k] for k in ("locus", "v_call", "j_call", "junction_aa")})
   # {'locus': 'TRB', 'v_call': 'TRBV7-2*04', 'j_call': 'TRBJ2-2*01',
   #  'junction_aa': 'CASSLAGGNTGELFF'}

``seqtype="aa"`` annotates amino-acid input. For a whole FASTQ use the CLI or
:func:`arda.rnaseq.pipeline.run` — ``annotate_sequences`` holds its batch in memory.

.. _examples-analysis:

Analysing the clonotype table
-----------------------------

.. important::

   **Read arda's TSVs with quoting off.** Every arda writer uses ``quote_style="never"`` and every
   reader ``quote_char=None``. Read one with polars' default quoting and a blank
   ``chimera_parents`` (written by ``arda correct --flag-chimeras``) becomes the two-character value
   ``""`` — every clonotype then reads as chimeric.

.. code-block:: python

   import polars as pl

   READ = dict(separator="\t", quote_char=None)
   clones = pl.read_csv("results/PT01.clones.tsv", **READ)

**Top clonotypes within a locus, with their clonal fraction.** The fraction is per locus, because
a TRB frequency computed over a table that also holds IGK is not a TRB frequency:

.. code-block:: python

   top = (
       clones.filter(pl.col("locus") == "TRB")
       .with_columns(
           (pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias("frequency")
       )
       .sort("duplicate_count", descending=True)
       .select("junction_aa", "v_call", "j_call", "duplicate_count", "frequency")
       .head(5)
   )

.. code-block:: text

   ┌─────────────────┬─────────────┬────────────┬─────────────────┬───────────┐
   │ junction_aa     ┆ v_call      ┆ j_call     ┆ duplicate_count ┆ frequency │
   ╞═════════════════╪═════════════╪════════════╪═════════════════╪═══════════╡
   │ CSQSGGFGADTQYF  ┆ TRBV29-1*03 ┆ TRBJ2-3*01 ┆ 2               ┆ 0.181818  │
   │ CSATPPDSWTGELFF ┆ TRBV20-1*07 ┆ TRBJ2-2*01 ┆ 2               ┆ 0.181818  │
   │ CASSLRGSYEQYF   ┆ TRBV7-8*01  ┆ TRBJ2-7*01 ┆ 2               ┆ 0.181818  │
   └─────────────────┴─────────────┴────────────┴─────────────────┴───────────┘

**Which chains the library actually carries**, in clonotypes and in reads:

.. code-block:: python

   chains = (
       clones.group_by("locus")
       .agg(pl.len().alias("clonotypes"), pl.col("duplicate_count").sum().alias("reads"))
       .sort("reads", descending=True)
   )
   # IGL 17/23, IGK 14/19, TRB 8/11, TRA 6/8, IGH 4/6 on the test fixture

**V-gene usage as a frequency, gene not allele.** Split on ``*`` rather than trusting a gene
column: ``v_call`` is allele-level by default and can be a comma-separated tie:

.. code-block:: python

   usage = (
       clones.filter(pl.col("locus") == "IGH")
       .with_columns(pl.col("v_call").str.split("*").list.first().alias("v_gene"))
       .group_by("v_gene")
       .agg(pl.col("duplicate_count").sum().alias("reads"))
       .with_columns((pl.col("reads") / pl.col("reads").sum()).alias("frequency"))
       .sort("reads", descending=True)
   )

Run ``arda rnaseq --call-level gene`` if you want the collapse done upstream, at clonotype-key
level, instead.

**Overlap between two repertoires**, on the amino-acid junction within one locus:

.. code-block:: python

   a = pl.read_csv("results/PT01.clones.tsv", **READ).filter(pl.col("locus") == "TRB")
   b = pl.read_csv("results/PT02.clones.tsv", **READ).filter(pl.col("locus") == "TRB")
   shared = a.join(b.select("junction_aa"), on="junction_aa")
   f = shared.height / min(a.height, b.height)     # overlap as a fraction of the smaller set

.. note::

   ``junction_aa`` includes the Cys104 and [FW]118 anchors; IMGT ``cdr3_aa`` excludes them and is
   two residues shorter. Joining one table's ``junction_aa`` to another's ``cdr3_aa`` silently
   finds nothing. MiXCR's ``CDR3`` column is AIRR ``junction``.

**Somatic hypermutation per read**, from the AIRR table — ``v_mutations`` is a comma-separated
list of the mutations found outside the junction, so its length is the count:

.. code-block:: python

   airr = pl.read_csv("results/PT01.airr.tsv", **READ)
   shm = (
       airr.filter(pl.col("locus").str.starts_with("IG"))
       .with_columns(
           pl.when(pl.col("v_mutations").is_null()).then(0)
           .otherwise(pl.col("v_mutations").str.split(",").list.len())
           .alias("v_mutation_count")
       )
       .group_by("locus")
       .agg(pl.col("v_mutation_count").mean().round(2).alias("mean_v_mutations"))
   )

:doc:`shm` explains what is and is not counted (the junction is excluded on purpose) and how to
recount an existing AIRR file with ``arda shm``.

Gate a run on its own report
----------------------------

Do not parse the log. ``<prefix>.arda.json`` carries the same numbers as structured data, and a
pipeline should assert on them:

.. code-block:: python

   import json

   rep = json.loads(open("results/PT01.arda.json").read())
   m = rep["map"]
   print(f"mapped {m['mapped_fraction']:.3f} of {m['total_reads']} reads, "
         f"{rep['correct']['clonotypes_out']} clonotypes, {rep['wall_seconds']} s")

   assert m["mapped_fraction"] > 0.01, "nothing mapped — wrong organism or wrong reference"
   assert rep["correct"]["clonotypes_out"] > 0, "no clonotypes"

.. code-block:: text

   mapped 0.343 of 1320 reads, 49 clonotypes, 2.47 s

A sample that arrived as several read groups adds ``read_groups`` and reports
``wall_seconds_max`` / ``wall_seconds_sum`` instead of a single ``wall_seconds`` for the map
stage: summing forty array tasks' wall time and calling it "wall seconds" would be a lie.
``reads_per_second``, ``peak_rss_mb``-style memory and every count are present either way.

``<prefix>.stats.tsv`` is the same run in long format, plus per-chain and per-gene QC —
``scope``/``key``/``metric``/``value``, which pivots:

.. code-block:: python

   qc = pl.read_csv("results/PT01.stats.tsv", **READ)
   per_chain = qc.filter(pl.col("scope") == "chain").pivot(
       on="metric", index="key", values="value"
   )

:doc:`qc` lists every scope and metric it writes.

Cohort QC in one table
----------------------

The QC tables are long format precisely so a cohort concatenates without a schema decision:

.. code-block:: python

   import glob
   import polars as pl

   cohort = pl.concat([
       pl.read_csv(p, separator="\t", quote_char=None).with_columns(
           pl.lit(p.split("/")[-1].removesuffix(".stats.tsv")).alias("sample")
       )
       for p in sorted(glob.glob("results/*.stats.tsv"))
   ])
   mapped = cohort.filter(
       (pl.col("scope") == "run") & (pl.col("metric") == "mapped_fraction")
   ).select("sample", "value")

A sample whose ``mapped_fraction`` is an order of magnitude below its cohort is the one to look at
first: wrong organism, wrong regime, or a library that did not work.

Export the reference for a genome browser
-----------------------------------------

The markup arda transfers is itself exportable — every in-frame V·J scaffold with IgBLAST-quality
FR1–4 / CDR1–3 coordinates, in three kinds × four formats:

.. code-block:: bash

   arda export-ref --kind scaffolds --locus TRB --format gff3 -o trb.gff3
   arda export-ref --kind alleles --locus IGH --format tsv -o igh_alleles.tsv

Coordinates are **1-based closed**, as AIRR and GFF3 both are, so they load unshifted.
:doc:`reference_export` covers the kinds, the formats and the round-trip guarantee.
