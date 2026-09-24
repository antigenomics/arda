# `ARDA_ASSIGN` — arda as a V(D)J assignment step for nf-core/airrflow

A drop-in alternative to `CHANGEO_ASSIGNGENES` + `CHANGEO_MAKEDB`: the same `[meta, reads]`
channel in, a spec-valid **AIRR Rearrangement TSV** out, plus a clonotype table.

arda does the IgBLAST work **once, offline**, when its reference is built — every in-frame V·J
germline scaffold with FR1–4 / CDR1–3 markup — then maps reads onto it with MMseqs2 and transfers
the markup through the alignment in C++. So there is no per-run germline database to stage:
`--reference_igblast`, `--reference_fasta` and `--fetch_germlines` are not consulted.

Pinned to **arda 2.27.0** (`environment.yml`, the `container` tag, the `Dockerfile`).
`tests/unit/test_nextflow_integration.py` asserts all four against `arda.__version__`.

---

## ⛔ Breaking change from the pre-2.28 module

The old module invented its own vocabulary. This one speaks airrflow's.

| before | now | why |
|---|---|---|
| `process ARDA` | `process ARDA_ASSIGN` | names the step it replaces |
| `--regime bulk\|amplicon\|default` | `--library_generation_method` | airrflow's own field; one name for the protocol, not two |
| `meta.regime` per sample | the samplesheet's `pcr_target_locus` / protocol | the sheet already carries it |
| `params.arda_organism` (or `params.genome`) | the samplesheet's **`species`** column | `arda_organism` survives as an override only |
| `params.arda_shm`, `arda_call_level`, `arda_ec_mode`, `arda_clonotype_key` | one `params.arda_args` passthrough | each is a *decision* about what a clonotype is, not a tuning knob; spelling them out in a pipeline config invited setting them unread |
| `params.arda_indel_rescue` | `params.arda_args = '--indel-rescue'` | same reason; it is a per-library call whose value tracks SHM load |
| iGenomes `params.genome` gate | none | arda's reference is IMGT-derived, so only the species matters — GRCh37 and GRCh38 are the same reference to it |

**Migrating.** Delete `--regime`; set `--library_generation_method` to the protocol you actually
used. Move any `arda_*` tuning params into `arda_args`. Make sure your samplesheet has the
`species` column airrflow's schema already requires.

| old `--regime` | new `--library_generation_method` |
|---|---|
| `amplicon` | `specific_pcr`, `specific_pcr_umi`, `dt_5p_race`, `dt_5p_race_umi` |
| `bulk` | `trust4` |
| `default` | `trust4` plus `--arda_args '--exact'` |

---

## Protocol → mode

The two speed paths in arda do **not** compose, and choosing the wrong one is a silent 2–4×
slowdown rather than an error. So the choice is *derived* from the protocol, and arda itself owns
which flags that implies.

| `library_generation_method` | runs | why |
|---|---|---|
| `specific_pcr`, `specific_pcr_umi` | `arda amplicon` | the read already spans V into J |
| `dt_5p_race`, `dt_5p_race_umi` | `arda amplicon` | same |
| `trust4` | `arda rnaseq` | whole-transcriptome; 0.02–3 % receptor |
| `sc_10x_genomics` | **refused** | see below |

### Single cell is refused, on purpose

`arda cells` exists, but its input is **one per-molecule UMI consensus FASTQ with the cell barcode
in the record name** — what `migec assemble` or Cell Ranger writes — not a raw 10x read pair. arda
does no barcode demultiplexing and no UMI collapse; both belong upstream. Mapping
`sc_10x_genomics` here would hand `arda cells` reads it cannot interpret, and the failure would
look like a bad repertoire rather than a wiring error. Run that path separately; see
`docs/singlecell.rst`.

---

## Samplesheet

airrflow's schema, unchanged. The columns this module reads:

| column | used for |
|---|---|
| `sample_id` | `meta.id`, the output prefix. Repeated rows **merge in row order** |
| `species` | `arda --organism`. `human` / `mouse` (arda also ships rat, rabbit, rhesus) |
| `single_cell` | `TRUE` is refused here — see above |
| `filename_R1`, `filename_R2` | the reads. Blank `filename_R2` is single-end |
| `subject_id`, `tissue` | carried as `arda qc batch`'s grouping labels |

⚠ **arda reads this sheet natively** — `arda.samples.read_sheet` speaks both the airrflow dialect
and nf-core's generic `sample` / `fastq_1` / `fastq_2`. That is what lets **one** samplesheet drive
this module, the Snakemake workflow in `integrations/snakemake/arda/`, and `arda cluster` (SLURM)
without a translation step anywhere.

**A sample may arrive in several files.** One FASTQ per lane is still one repertoire and must give
one clonotype table. arda maps each read group and concatenates before the global stages, which is
byte-identical to the same reads in one file — so pass them all rather than `cat`-ing.
Never sorted by filename: `A_L010` sorts before `A_L002`, and the clonotype fold is not
permutation-invariant.

---

## Wiring it in

```groovy
include { ARDA_ASSIGN } from '../modules/local/arda/main'

ARDA_ASSIGN( ch_reads )            // [ meta, reads ] — the same channel CHANGEO_ASSIGNGENES takes
ch_versions = ch_versions.mix( ARDA_ASSIGN.out.versions.first() )
```

```bash
nextflow run . \
    --input samplesheet.tsv \
    --library_generation_method specific_pcr_umi \
    --outdir results \
    -profile conda
```

`includeConfig` the `nextflow.config` beside this file for `publishDir` and the `arda_args`
passthrough.

### Execution profiles

- **conda** — works from `environment.yml` as shipped.
- **docker / singularity** — build the image beside this module, push it to your own registry, and
  point the module at the tag:

  ```bash
  docker build -t arda-mapper:2.27.0 integrations/nextflow/arda
  docker tag  arda-mapper:2.27.0 <your-registry>/arda-mapper:2.27.0
  docker push <your-registry>/arda-mapper:2.27.0
  ```

  then `withName: 'ARDA_ASSIGN' { container = '<your-registry>/arda-mapper:2.27.0' }`.

---

## Resources

arda is CPU-bound — the MMseqs2 search dominates. ~40–50k reads/s on 32 cores; a full-depth
~100 M-read sample takes ~45 min. The module is `label 'process_high'`.

**Memory scales with repertoire richness, not with FASTQ size.** Stage 1 (`map`) is flat at
**300–650 MB** at any depth because it streams. Stage 3 (`correct`) holds the whole clone set:
it peaked **2,071.7 MB** on a B-cell-rich tumour (28,444 clonotypes, 105 M reads), while a colder
**139 M**-read sample — more reads, almost no repertoire — peaked **549 MB**. Budget ~4 GB.

---

## What is deliberately not wired

| airrflow param | why |
|---|---|
| `productive_only` | arda has **no** productive filter and does not pretend to. It emits the AIRR `productive` column and leaves the decision to the consumer — the same stance as its QC surface, which flags and never filters. airrflow's own filter runs downstream, unchanged. |
| `reference_igblast`, `reference_fasta`, `fetch_germlines` | arda's reference is built offline, once. Set the species instead. |
| `clonal_threshold` | arda's Stage 2 is an abundance + quality error model over exact junctions, not a distance threshold over a clone. Different question; airrflow's clonal step runs on arda's AIRR output as usual. |

Each is **named** rather than silently ignored: a parameter that is accepted and does nothing is
the failure mode this project keeps hitting, and `tests/unit/test_nextflow_integration.py` asserts
they stay named and stay unread.
