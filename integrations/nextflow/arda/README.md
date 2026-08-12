# ARDA — Nextflow module

A drop-in, nf-core-style local module that runs arda's RNA-seq mode on each sample and publishes
per-sample **AIRR clonotype tables** to `${params.outdir}/arda/`. It wraps a single
`arda <mode>` call (map + assemble + correct) and emits a `versions.yml`, so it composes with
any DSL2 pipeline the same way STAR/Salmon/fastp do.

Pinned to **arda 2.20.0** (`environment.yml`, the `container` tag, and the `Dockerfile`).

> ⛔ **2.16.0 is a hard minimum, and it is a BREAKING one.** `arda rnaseq run` — the command every
> earlier version of this module invoked — was removed there. The regime is now the **command
> name** (`arda rnaseq` / `arda amplicon`), and each mode owns its own speed configuration, so this
> module and the CLI move together: an older arda fails with *"Got unexpected extra argument
> (run)"*, and an older module against 2.16.0 fails the same way.
>
> It is also the first release with `--shm` (SHM scoped to the framework rather than to the V/J
> segment) and with per-mode `--ec-mode` defaults.
>
> The benchmark tables below are labelled with the version they were **measured** on; the *pin* is
> the release that ships the commands.

## What it produces (per sample `<id>`)

| file | contents |
|---|---|
| `<id>.clones.tsv` | corrected clonotype table: `junction, junction_aa, v_call, j_call, c_call, locus, duplicate_count, consensus_count, d_call, d2_call, d_support, d2_support` |
| `<id>.airr.tsv` | mapped reads, AIRR Rearrangement schema |
| `<id>.assembled.airr.tsv` | Stage-3 long-CDR3 reads rescued by contig assembly |
| `<id>.arda.json` | run report (reads mapped, per-locus counts, isotype/constant, timing, peak RSS, `fast_fraction`) |

## ⛔ Pick the regime — it is one parameter and it is easy to get backwards

arda has two tuning paths and **they do not compose**. Choosing the wrong one is not an error; it
is a silent 2–4× slowdown. The module therefore selects the combination by name:

| `--regime` | arda command | speed configuration it implies | use for |
|---|---|---|---|
| `amplicon` | `arda amplicon` | `--two-pass --fast-segments --v-only-on-segment` | targeted RepSeq / 5′RACE libraries |
| `bulk` | `arda rnaseq` | `--prefilter` | whole-transcriptome RNA-seq (the default) |
| `default` | `arda rnaseq --exact` | *(none)* | the shipped one-pass path, for reproducing older runs |

⛔ **`--two-pass` on its own is a LOSS** — 0.762× on bulk and 0.87× on an IGH amplicon — and it is
no longer reachable by accident: arda owns the combination behind the mode name. Do not
hand-assemble these flags in `ext.args`.

⚠ `default` is **not** `arda rnaseq`. That mode turns `--prefilter` on, which costs ~0.15 % of
mapped reads (122 bulk datasets; up to 2.46 % on one library). `--exact` is what reproduces the
pre-2.16.0 default output.

A sheet may legitimately mix the two library types: put a `regime` key in the meta map and it wins
over `params.regime` for that sample.

⚠ `regime` is a short name in a shared params namespace. If your pipeline already uses it for
something else, rename the param in `nextflow.config` and in the one `params.getOrDefault('regime',
…)` call in `main.nf`; the per-sample `meta.regime` route is unaffected either way.

If you are unsure which regime a library wants, run one sample either way and read `fast_fraction`
from `<id>.arda.json`: it is the predictor, not the library's name. High (~.85) means the amplicon
regime pays; low (~.05) means it is overhead.

## Runtime & resources

arda is **CPU-bound** — the aligner dominates — so give the process cores. `--threads` follows
`task.cpus` automatically.

**TRA amplicon, 100,000 reads, 8 threads, same input file in the same job:**

| tool | config | wall (s) | CPU (s) | peak RSS (MB) |
|---|---|---|---|---|
| arda 2.11.1 | `--regime amplicon` | **5.35** | **12.73** | **631** |
| MiXCR 4.7.0 | `align --preset rna-seq --species hsa` | 5.90 | 45.24 | 3,027 |

1.10× faster on wall, **3.6× less CPU, 4.8× less RSS**.

**Bulk RNA-seq, 100,000 reads:**

| tool | wall (s) | CPU (s) | peak RSS (MB) |
|---|---|---|---|
| arda 2.11.1 | 2.51 | 5.4 | 234 |
| MiXCR 4.7.0 | 4.54 | 31.8 | 3,022 |
| TRUST4 | **1.91** | **4.36** | **192** |

⚠ TRUST4's stage here is **candidate read extraction**, not a per-read AIRR record with a junction;
arda's and MiXCR's are. The three numbers are not like-for-like work.

**IGH RepSeq amplicon, 100,000 pairs, 32 threads** — what the regime is worth on a real repertoire:

| dataset | config | wall (s) | peak RSS (MB) |
|---|---|---|---|
| IGH_repertoire | `--regime default` | 316.44 | 4,018 |
| IGH_repertoire | `--regime amplicon` | **76.25** | **1,479** |
| IGH_naive | `--regime default` | 305.32 | 3,736 |
| IGH_naive | `--regime amplicon` | **64.86** | **1,363** |

That is **4.15×** and **4.71×**, at roughly a third of the memory.

### Memory budget

**Stage 1 (`map`) is flat: 300–650 MB at any read depth** — it streams. What scales is **Stage 3
(`correct`), which holds the whole clone set in memory**: it peaked **2,071.7 MB** on a B-cell-rich
tumour with **28,444 clonotypes** from 105 M reads, while a *colder* **139 M-read** sample — more
reads, almost no repertoire — peaked **549 MB**.

**Budget ~4 GB per task**, and size it by expected repertoire richness, not by FASTQ size. If you
must cap tightly, `--no-assemble` (via `--arda_args`) keeps the run on the flat mapping-only
profile, at the cost of the long CDR3s no single read spans.

The module is labelled `process_high` because that is the nf-core label that hands out cores — but
that label also reserves 72 GB, roughly 18× what arda needs. arda wants **many cores and little
RAM**, which no standard label expresses, so set both explicitly:

```groovy
withName: 'ARDA' { cpus = 32; memory = 8.GB; time = 4.h }
```

## Accuracy

Against an IgBLAST truth on the same 100,000-read TRA amplicon:

| metric | arda 2.11.1 | MiXCR 4.7.0 |
|---|---|---|
| `v_gene` recall | .9867 | **.9973** |
| `v_gene` precision | **.9996** | .9978 |
| `v_allele` resolved | **.9868** | n/a |
| `j_gene` recall | .9892 | **.9904** |
| `j_gene` precision | .9953 | **.9995** |
| junction precision among emitted | .99919 | **.99991** |

arda is **more precise on the V call and declines rather than guessing**; MiXCR recalls slightly
more. MiXCR emits `*00` and so makes no allele call at all — across 25 cluster datasets arda's
median `v_allele` is **.9763** resolved (**.8328** by exact string, the difference being ambiguous
allele tie-lists, which are a scoring convention, not a call).

⛔ **A V/J boundary disagreement *inside* a junction is not an error.** V(D)J recombination is
probabilistic — exonuclease chew-back plus N/P-nucleotide addition mean the V-end / NDN / J-start
partition of a junction is frequently not identifiable from sequence alone, and the ground truth is
unknown. Overlapping V/J/NDN assignments are acceptable. What is checkable, and what the table
above scores, is the junction's **outer bounds** (Cys104 and [FW]118), the **gene/allele calls**,
and whether a tool **invents a junction it has no anchor for**.

## Requirements

arda is pip-installable and needs the `mmseqs2` binary — both are declared in `environment.yml`.

- **`-profile conda`** works out of the box (Nextflow builds the env from `environment.yml`) — once
  arda 2.20.0 is on PyPI; see the note at the top.
- **`-profile docker`/`singularity`**: build the image from the `Dockerfile` here, push it to your
  registry, and point the module's `container` at it (see the Dockerfile header). A pinned image is
  the reproducible choice for a shared pipeline.

⛔ **The aligner is pinned, deliberately.** An mmseqs index is only reusable by the release that
built it, and a cluster's cached index marker can differ from the shipped one. Left alone, arda
would reject the precompiled reference index and rebuild a private cache per task — or auto-fetch a
third build — with no error, and with results that are not comparable to anyone else's. The module
therefore exports **`ARDA_MMSEQS`**, pointing at whatever mmseqs the task's own conda env or
container provides. Set `--arda_mmseqs /abs/path/to/mmseqs` only if mmseqs lives outside that
environment. `environment.yml` pins `mmseqs2 =18.8cc5c` exactly for the same reason.

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
# amplicon library (targeted RepSeq / 5'RACE)
nextflow run test.nf -c nextflow.config -profile conda \
    --input amplicon.csv --outdir results_amplicon --regime amplicon

# bulk RNA-seq
nextflow run test.nf -c nextflow.config -profile conda \
    --input bulk.csv --outdir results_bulk --regime bulk
```

Both write `results_*/arda/<sample>.{clones.tsv,airr.tsv,assembled.airr.tsv,arda.json}`.

`amplicon.csv` is single-end here, `bulk.csv` paired:

```
sample,fastq_1,fastq_2
tra_amplicon,/data/tra.fastq,
```

```
sample,fastq_1,fastq_2
bulk_rnaseq,/data/SRR5233637_1.fq,/data/SRR5233637_2.fq
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
   module configs are pulled in). This sets the params defaults, `ext.args`, the genome gate and the
   `${params.outdir}/arda` publishDir.

4. **Register the toggles** so strict schema validation accepts them: add `run_arda = false` to the
   `params { }` block in `nextflow.config`, and matching properties in `nextflow_schema.json` for
   `run_arda` (boolean), `regime` (string enum `bulk|amplicon|default`), `arda_indel_rescue`
   (boolean), `arda_organism` (string), `arda_mmseqs` (string), `arda_args` (string).

5. **Container override** (only for `-profile docker/singularity/<your-profile>`): add
   `withName: 'ARDA' { container = '<your-registry>/arda-mapper:2.20.0' }` to your deployment
   config (e.g. `conf/<profile>.config`), exactly as the other tools' images are pinned there.

Run with `--run_arda`:
```bash
nextflow run . -profile <your-profile> --input samplesheet.csv --outdir results \
    --run_arda --regime bulk
```

## Organism follows the genome automatically

The shipped `nextflow.config` reads `params.genome` — the iGenomes assembly key nf-core sets from
`--genome`. arda's reference is IMGT-derived, so only the **species** matters; GRCh37 and GRCh38 are
the same reference to it.

| `--genome` | ARDA runs? | `--organism` |
|---|---|---|
| `GRCh38` | yes | `human` |
| `GRCh37` | yes | `human` |
| `GRCm39` | yes | `mouse` |
| `GRCm38` | yes | `mouse` |
| unset | yes | `human` |
| any other | **skipped** (pipeline still completes) | — |

Other assemblies are skipped because arda ships full references only for human and mouse. Setting
`--arda_organism` explicitly overrides both the gate and the mapping.

## Tuning

The regime is the main knob (above). Everything else goes through `--arda_args`, appended last:

| goal | flag |
|---|---|
| force an organism | `--arda_organism mouse` |
| merge overlapping mates first | `--arda_args '--reconstruct'` |
| keep every mapped read (recall-max) | `--arda_args '--min-score 0'` |
| cap memory harder / looser | `--arda_args '--kmer 11'` (or `--kmer 0` for the mmseqs default) |
| skip Stage 3 assembly (flat memory) | `--arda_args '--no-assemble'` |
| gapped rescue for hypermutated IG | `--arda_indel_rescue` (**requires `--regime amplicon`**) |

`--indel-rescue` is deliberately *not* part of the amplicon preset: its value tracks somatic
hypermutation load (+181 reads on a hypermutated repertoire, **−14 on a naive one**), so it is a
per-library call. The module raises an error rather than accepting it in a non-amplicon regime,
where arda would ignore it silently.

`--threads` is wired to `task.cpus` automatically, and the regime flags are owned by the module — do
not set either in `--arda_args`.
