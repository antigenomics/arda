# arda — Snakemake workflow

> **The sample sheet is nf-core/airrflow's.** This workflow, the Nextflow module
> (`integrations/nextflow/arda/`) and `arda cluster` (SLURM) all parse it through
> `arda.samples.read_sheet`, which speaks both airrflow's `sample_id` / `filename_R1` /
> `filename_R2` + AIRR metadata and nf-core's generic `sample` / `fastq_1` / `fastq_2`.
> One sheet, three runners, no translation step. A sheet carrying **both** id columns is refused —
> it does not say which column names the repertoire.
>
> `species` is read per sample (a cohort may mix organisms); `--config organism=` is the fallback
> for a sheet that does not carry it. `single_cell=TRUE` is refused: `arda cells` takes a
> per-molecule UMI consensus FASTQ, not a read pair.

Sample-sheet driven, **one job per read group**.

```sh
snakemake -s Snakefile --config samples=sheet.tsv outdir=results -c 32
snakemake -s Snakefile --config samples=sheet.tsv outdir=results --profile slurm
```

## Why the read group is the unit

Stage 1 (`arda map`) is per-read and shards perfectly, so every read group of every sample is an
independent job: a sheet of 5 samples × 4 lanes is **20 concurrent jobs, not 5**. Stages 2–3
(`assemble`, `correct`) are global **per sample** — a clone split across read groups would be
counted once per group, and contigs that tile across them would never be built — so they run once
per sample, in `arda cluster reduce`, over that sample's merged Stage-1 AIRR.

Read groups are consumed in **sheet order**, which is what makes each sample's result
byte-identical to its reads concatenated into one file. Do not sort them: `A_L010` sorts before
`A_L002`, and the clonotype fold is not permutation-invariant.

## The sheet

The same one `arda rnaseq --samples` and `arda cluster plan --samples` read — nf-core's spelling,
so an existing nf-core samplesheet works unmodified:

```tsv
sample	fastq_1	fastq_2
PT01	PT01_S1_L001_R1_001.fastq.gz	PT01_S1_L001_R2_001.fastq.gz
PT01	PT01_S1_L002_R1_001.fastq.gz	PT01_S1_L002_R2_001.fastq.gz
PT02	PT02_S2_L001_R1_001.fastq.gz	PT02_S2_L001_R2_001.fastq.gz
```

Two samples, three read groups. `PT01` is one repertoire and gets one `PT01.clones.tsv`.
`.csv` is read as CSV, anything else as TSV. Leave `fastq_2` blank for single-end. Relative paths
resolve against the sheet's own directory, so a sheet travels with its data.

## Config

| key | default | what |
|---|---|---|
| `samples` | *required* | the sheet |
| `outdir` | `results` | per-sample outputs land here |
| `workdir` | `<outdir>/work` | Stage-1 parts, one dir per sample |
| `regime` | `rnaseq` | `rnaseq` (bulk) or `amplicon` (targeted RepSeq) |
| `organism` | `human` | reference organism |
| `arda` | `arda` | the executable, for a pinned env |
| `map_threads` | `8` | cores per Stage-1 job |
| `reduce_threads` | `8` | cores per Stage-2/3 job |
| `extra` | `""` | extra flags appended to `arda map` |

**Raise `map_threads` before raising `-c`.** arda is CPU-bound on the MMseqs2 search and threads
internally, so a few big jobs beat many small ones.

**Name the regime, never the flags.** arda's two speed levers do not compose and each is a loss in
the other's regime; and any denoising preset but `fast` reads a Stage-1 column that only
`--junction-quality` writes. The workflow asks `arda.cluster.regime_flags` for both halves rather
than restating them, so the two cannot drift apart.

## Outputs

Per sample, under `outdir`:

| file | what |
|---|---|
| `<sample>.clones.tsv` | corrected clonotype table |
| `<sample>.airr.tsv` | Stage-1 mapped reads (AIRR Rearrangement) |
| `<sample>.assembled.airr.tsv` | Stage-3 long-CDR3 contigs |
| `<sample>.arda.json` | merged run report |
| `<sample>.stats.tsv` | run QC, long format |

## Resources

`arda map` is flat at **300–650 MB** whatever the depth. `arda cluster reduce` holds the clone
set: budget ~4 GB, more for a B-cell-rich sample (2,071.7 MB measured on 28,444 clonotypes from
105 M reads). A SLURM profile wants roughly:

```yaml
set-resources:
  map_read_group:
    mem_mb: 8000
    runtime: 240
  reduce_sample:
    mem_mb: 16000
    runtime: 480
```

## Not Snakemake?

`arda cluster plan --samples sheet.tsv --work-dir W` writes the same DAG as two TSV manifests —
one row per read group, one per sample — and prints the two commands filled in. Any scheduler can
drive arda from those and nothing else. `arda cluster submit-samples` renders them as SLURM arrays
directly.
