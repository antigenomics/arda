# ARDA — Nextflow module

A drop-in, nf-core-style local module that runs arda's RNA-seq mode on each sample and publishes
per-sample **AIRR clonotype tables** to `${params.outdir}/arda/`. It wraps a single
`arda rnaseq run` call (map + assemble + correct) and emits a `versions.yml`, so it composes with any DSL2
pipeline the same way STAR/Salmon/fastp do.

## What it produces (per sample `<id>`)

| file | contents |
|---|---|
| `<id>.clones.tsv` | corrected clonotype table: `junction, junction_aa, v_call, j_call, c_call, locus, duplicate_count, consensus_count, d_call, d2_call, d_support, d2_support` |
| `<id>.airr.tsv` | mapped reads, AIRR Rearrangement schema |
| `<id>.assembled.airr.tsv` | Stage-3 long-CDR3 reads rescued by contig assembly |
| `<id>.arda.json` | run report (reads mapped, per-locus counts, isotype/constant, timing, peak RSS) |

## Runtime & resources

arda is **CPU-bound** (the MMseqs2 search dominates) and **very low-memory** (< 400 MB peak RSS,
independent of read depth) — so give the process cores, not RAM. Measured on bulk tumor RNA-seq,
32 cores:

| reads | cores | wall time | throughput | peak RSS |
|---|---|---|---|---|
| 104.9 M (52.4 M pairs, 2×150) | 32 | 44 min | ~39,600 reads/s (~2.4 M/min) | 371 MB |

Throughput scales roughly linearly with cores. The module is labelled `process_high`; for full-depth
samples give it 16–32 cpus (`withName: 'ARDA' { cpus = 32 }`). `--threads` follows `task.cpus`.

## Requirements

arda is pip-installable and needs the `mmseqs2` binary — both are declared in `environment.yml`.

- **`-profile conda`** works out of the box (Nextflow builds the env from `environment.yml`).
- **`-profile docker`/`singularity`**: build the image from the `Dockerfile` here, push it to your
  registry, and point the module's `container` at it (see the Dockerfile header). A pinned image is
  the reproducible choice for a shared pipeline.

## Try it standalone (one process, no full pipeline)

```nextflow
// test.nf
include { ARDA } from './main.nf'

workflow {
    Channel
        .fromPath(params.input)                          // a CSV: sample,fastq_1,fastq_2
        .splitCsv(header: true)
        .map { row ->
            def single = !row.fastq_2
            [ [id: row.sample, single_end: single],
              single ? [file(row.fastq_1)] : [file(row.fastq_1), file(row.fastq_2)] ]
        }
        .set { reads }
    ARDA(reads)
}
```

```bash
nextflow run test.nf -profile conda \
    --input samplesheet.csv --outdir results
# -> results/arda/<sample>.clones.tsv, .airr.tsv, .arda.json
```

## Drop into an nf-core/rnaseq (v3.x) pipeline

The module consumes the same per-sample FASTQ channel the aligners do, so it needs no sample-sheet
changes. Five edits, all mirroring how an existing tool is wired:

1. **Copy** this folder to `modules/local/arda/` in your pipeline checkout.

2. **Include and call** it in `workflows/rnaseq/main.nf`. The cleanest input is the trimmed,
   aligner-independent FASTQ channel `ch_strand_inferred_filtered_fastq` (already `[meta, reads]`):
   ```nextflow
   include { ARDA } from '../../modules/local/arda'
   // ... after the trim/QC subworkflow ...
   if (params.run_arda) {
       ARDA(ch_strand_inferred_filtered_fastq)
       ch_versions = ch_versions.mix(ARDA.out.versions.first())
   }
   ```

3. **Aggregate config**: add
   `includeConfig "../../modules/local/arda/nextflow.config"`
   to the `includeConfig` block at the top of `workflows/rnaseq/nextflow.config` (where the other
   module configs are pulled in). This sets `ext.args` and the `${params.outdir}/arda` publishDir.

4. **Register the toggle** so strict schema validation accepts it: add `run_arda = false` to the
   `params { }` block in `nextflow.config`, and a matching boolean property in `nextflow_schema.json`
   (copy any existing `skip_*`/`run_*` boolean entry as a template).

5. **Container override** (only for `-profile docker/singularity/<your-profile>`): add
   `withName: 'ARDA' { container = '<your-registry>/arda-mapper:2.5.1' }` to your deployment config
   (e.g. `conf/<profile>.config`), exactly as the other tools' images are pinned there.

Run with `--run_arda`:
```bash
nextflow run . -profile <your-profile> --input samplesheet.csv --outdir results --run_arda
```

## Tuning

All arda flags pass through `ext.args` (set in `nextflow.config`). Common ones:

| goal | `ext.args` |
|---|---|
| mouse reference | `--organism mouse` |
| merge overlapping mates first | `--organism human --reconstruct` |
| keep every mapped read (recall-max) | `--organism human --min-score 0` |
| cap memory harder / looser | `--organism human --kmer 11` (or `--kmer 0` for the mmseqs default) |

`--threads` is wired to `task.cpus` automatically; do not set it in `ext.args`.
