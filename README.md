<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/arda_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/arda_light.svg">
    <!-- Absolute PNG fallback: PyPI strips <picture>/<source> and cannot render a relative or
         raw-served SVG, so the logo must be an absolute-URL raster here. GitHub uses the SVG sources. -->
    <img alt="arda" src="https://raw.githubusercontent.com/antigenomics/arda/master/assets/arda_dark.png" width="340">
  </picture>
</p>

<h1 align="center">arda — Antigen Receptor Domain Annotation</h1>

<p align="center">
  <a href="https://pypi.org/project/arda-mapper/"><img alt="PyPI" src="https://img.shields.io/pypi/v/arda-mapper"></a>
  <a href="https://github.com/antigenomics/arda/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/arda/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/arda/"><img alt="docs" src="https://github.com/antigenomics/arda/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <img alt="license" src="https://img.shields.io/badge/license-GPLv3-green">
</p>

**Versatile, fast, exact** FR/CDR annotation of **TCR** and **BCR** sequences —
mRNA and protein in FASTA, and reads in FASTQ from both **amplicon** and **bulk
RNA-seq** — for nucleotide *and* amino-acid input, across all loci at once.

`arda` does the expensive IgBLAST work **once, offline** — building a pre-aligned
reference database of every in-frame V·J germline scaffold with FR1–4 / CDR1–3
markup — then at runtime maps your sequences to that database with **MMseqs2** and
transfers the markup through the alignment in a small **C++** hot path. The result
is a **spec-valid [AIRR](https://docs.airr-community.org/) Rearrangement** annotation
that matches IgBLAST (98–99.7% region concordance on real GenBank mRNA), from a plain
CLI + Python library — no Docker, no workflow engine.

It also annotates records that have **no read behind them** — a CDR3 amino acid, a V
call and a J call, as in a VDJdb row — marking up which residues each germline
templates, repairing the junction, and inferring the D gene from the junction's length.

## Install

```bash
pip install arda-mapper   # from PyPI (imports as `arda`); binary wheels ship the C++ extension
```

`mmseqs2` (the search backend) is fetched/managed by arda at runtime. For development — and to
get the committed germline references on disk — use `setup.sh`:

```bash
bash setup.sh            # uv .venv, builds the C++ extensions, fetches IgBLAST + static mmseqs
source .venv/bin/activate
```

Needs [uv](https://docs.astral.sh/uv/). Flags: `--build-db` (rebuild references after install),
`--tests` (run the unit + synthetic suites — and **fail** if they fail). It installs
`.[test,dev]`, wipes any stale `build/` first, and then verifies three things that a bare
`import arda` does **not**: that `_markup`, `_segmap` and `_denoise` actually compiled (arda falls
back to a pure-Python markup path otherwise, so a failed build looks like a successful install and
surfaces later as a silent slowdown), that `mmseqs` resolves, and that every mode and stage
command resolves on the CLI. The committed `database/vdj/<organism>/`
references mean **most users never need to build anything**. A `pip install arda-mapper`
with no source checkout **auto-fetches** the curated references into `~/.cache/arda` on
first use and builds the MMseqs2 index there — no `$ARDA_HOME`, no build step
(set `ARDA_NO_AUTO_FETCH` for air-gapped runs with a pre-populated cache).

Supported organisms: **human, mouse** (full IG + TR), **rat, rabbit, rhesus_monkey**
(IG only — IgBLAST ships no TR internal annotation for these).

## The three modes

The speed configurations are **regime-specific and do not compose**, and picking the wrong one is
slower than picking none. So the regime is the **command name**, and arda owns what it implies:

| command | the library it is for | configuration it carries | denoising default |
|---|---|---|---|
| `arda rnaseq` | whole-transcriptome bulk RNA-seq (0.02–3 % receptor) | `--prefilter` | `--ec-mode rnaseq` |
| `arda amplicon` | targeted RepSeq / 5′RACE (reads span V into J) | `--two-pass --fast-segments --v-only-on-segment` | `--ec-mode amplicon` |
| `arda cells` | single cell — **UMI consensus per molecule**, barcode in the record name | reference-free per-cell assembly, then chain pairing | — (clonotypes are per cell) |

```bash
arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # bulk / whole-transcriptome
arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # targeted RepSeq / 5'RACE
arda cells    asm/PBMC.consensus.fq.gz -p out/PBMC            # single cell (see below)
arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/ --exact   # no speedups at all
```

`rnaseq` and `amplicon` run `map` → `assemble` → `correct` over raw paired FASTQ and write
`<prefix>.airr.tsv` (one AIRR row per mapped read), `<prefix>.clones.tsv` (the clonotype table),
`<prefix>.assembled.airr.tsv`, `<prefix>.arda.json` (run report) and `<prefix>.stats.tsv` /
`.stats.json` (run QC). Progress goes to **stderr**; the output paths, one per line, to
**stdout**.

`arda cells` is the single-cell entry point and starts one step later: its input is **one UMI
consensus per molecule** with the cell barcode in the record name — what `migec assemble` writes.
arda does no demultiplexing, no barcode correction, no UMI collapse and no cell calling; the
upstream tool does all four. From there it assembles each cell's contigs **with no germline
reference**, annotates them, pairs the chains and flags multiplets. Against Cell Ranger on
`sc5p_v2_hs_PBMC_1k` VDJ-T, **98.2 %** of its 943 CDR3s appear verbatim inside one of arda's
contigs. Guide: [single cell](https://docs.isalgo.dev/arda/singlecell.html).

Never: Until 2.16.0 the only entry point was `arda rnaseq run`, used for amplicon as well, with the
regime spelled out as four loose flags. **`--two-pass` alone is a LOSS** — 0.762× on bulk, 0.87× on
an IGH amplicon — and it was the one flag that entry point exposed for four releases. Naming the
mode makes the dominated combination unreachable by accident. `arda rnaseq run` no longer exists.

The predictor is not the library's name but whether a read hits **both** a V and a J segment —
`fast_fraction` in the run report. Primer-anchored amplicon reads do; bulk reads land anywhere in a
transcript and mostly do not, which is why the segment path is overhead there and the 16-mer
prefilter (a scan-term optimisation) is the bulk lever instead.

⚠ `arda rnaseq` enables `--prefilter`, which costs **~0.15 % of mapped reads** (122 bulk datasets;
up to 2.46 % on one library), concentrated in J→C and hypermutated IGH. Use `--exact` where that
matters. `--indel-rescue` (amplicon only) reroutes indel-bearing reads to the gapped path; its
value tracks SHM load, so it stays a per-library call and never rides the preset — and arda now
**refuses** it outside the amplicon mode rather than ignoring it. Per-flag measurements:
`arda amplicon --help`, `arda map --help` and the
[usage guide](https://docs.isalgo.dev/arda/usage.html).

### Flags that add or scope a column

`map` → `assemble` → `correct` are separate commands as well as stages inside a mode
([detail below](#the-three-stages-in-detail)); these are the flags that change what they emit:

| stage / flag | what it does |
|---|---|
| `--assemble` *(default on)* | contig assembly for long CDR3s no single 100–150 bp read spans; recovers ~95 % of the abundant long clones a filter-only pass misses |
| `--shm framework` *(default)* | SHM (`v_identity`, `v_mutations`, `j_mutations`) scoped **outside the junction** — IGH/IGK/IGL is where it is real |
| `--shm both` | also emit the pre-2.16.0 junction-inclusive values as `*_full` columns |
| `--isotype` *(default on)* | IGH isotype: `c_call`/`c_class`, voted per fragment then per clonotype, reported as **class** never subclass |
| `--call-level gene` | drop the allele suffix before the clonotype key, collapsing allele-level call splits |
| `--map-d` *(default on)* | D and tandem D-D alignment into the junction |
| `map --junction-quality` | Phred+33 string over exactly the bases of `junction`; needed by `correct --min-junction-q` and by `stats`' junction-quality metrics |
| `map --mutation-quality` | Phred behind each `v_mutations` / `j_mutations` entry, comma-joined and one-for-one; what `stats` scores `allele_candidate` on |
| `map --cell-from <dialect>` | lift the cell barcode out of `sequence_id` into `cell_id`; arda does **no** demultiplexing and no barcode correction — the upstream tool did both |
| `correct --cell-from migec` | AIRR `umi_count`: **distinct UMIs**, not distinct records, so a saturated barcode migec split into two molecules counts once. Omitted — never 0 or 1 — when the dialect names no UMI or the ids do not parse |

`arda shm -i in.airr.tsv -o out.airr.tsv` does the SHM recount standalone, needing **no reference
and no re-map** — the germline anchors are already in the file.

`arda scenarios -i clones.tsv -o d_prior.tsv` estimates the **generative model of V(D)J
recombination** from nucleotide junctions — trimming and insertion distributions, P(D|J) — by EM
over recombination scenarios. `arda.dpost` currently marginalises a model borrowed from OLGA;
this is how arda estimates its own, and the output is the same long table, drop-in.

`arda shm-model -i mapped.airr.tsv -o shm.tsv --locus IGH` fits the other half of that picture:
**where a mutation is expected**, as `P(substitution | 5-mer germline context, region)`. It reads
what every mode already writes under the default `--shm framework` and takes under a second on a
bulk library.

⛔ **Not a per-allele-per-position table**, and the measurement is why. Across two donors a
per-allele-per-position profile transfers at Pearson **r .26–.56** against **r .74–.79** for 5-mer
context; within one donor the positional table repeats at r ≈ .99, because it is a portrait of
that donor's expanded clones rather than a property of the allele. AID's own WRCY/RGYW motifs
carry **4.66–5.00×** the rate of every other covered position in all four libraries measured.
⚠ The fitted **scale** is the sample's, not the model's — two donors differ **1.6×** in overall
rate with the same shape — so it is written as provenance and a consumer rescales.

`arda scenarios --shm-model shm.tsv` is what reads it. The templated V length was bounded by an
**exact** match, so on IGH one substitution in the V tail forced the rest of that tail to be
explained as insertion — and `delV` and `insVD` are the distributions being fitted. Measured on
26,619 real IGH junctions with the model fitted on a **different donor**: `delV` mean
**4.148 → 2.708** (its mode 1 → 0), `insVD` **12.065 → 10.972**, and `delJ` **12.247 → 12.253**,
which does not move because the model reaches the V side only. Leave the flag off for TR and
unmutated IG, where it costs ~5 % wall and changes nothing.

Never: a scenario is **not identifiable from sequence** — 4,346 tuples reproduce one real human
TRB junction exactly — so the counts are expected counts summed over scenarios, never one MAP
reading. And an insertion costs its own sequence (`0.25^len`), not just its length: without that,
EM explains the whole junction as insertion. [Detail](https://docs.isalgo.dev/arda/scenarios.html).

### Accuracy regimes: which knob for which question

Separate from the speed flags, and set from **what you are going to do with the answer**. All are
off by default; the shipped output does not move unless you ask.

| question | flags | what it buys |
|---|---|---|
| Repertoire-level D usage | *(default, `E ≤ 0.2`)* | highest D recall |
| A D call you will act on, or a tandem D-D you will report | `--d-max-evalue 0.01` | gene agreement vs IgBLAST **.9765 → .9985** (TRB amplicon), at ~⅓ the call rate |
| Low-frequency variant recovery (spike-in, MRD, monoclonal control) | `map --junction-quality` + `correct --error-rate 1e-5 --ec-mode accurate` | keeps both published MIGEC spike-ins **and** monoclonal purity **.96034 → .99530** |
| SHM / lineage trees | *(default)* | `v_mutations`/`j_mutations` in germline coordinates, `+36 ms` per 100 k reads |
| Monoclonal QC / cell-line purity | `--ec-mode amplicon --clonotype-key junction` | Jurkat **90 → 10** clonotypes, TRB purity **.98963 → .99990**, **reads unchanged at 14,531**, and **98.50 %** of them on the two published clones |
| A targeted library that is deep | `--ec-mode amplicon` | quality-directed rescue at 12 subs / 50× — reaches the class the abundance model structurally cannot |
| Bulk RNA-seq, where singletons are mostly real | `--ec-mode rnaseq` | the same rescue kept narrow (6 subs / 200×) |

`--d-max-evalue` is a **recall/precision dial**, not a bug fix: the shipped 0.2 is deliberately the
loosest band, because dropping two thirds of the D calls is the wrong default for a repertoire
tool. `--ec-mode accurate` is `--min-junction-q 20` — it judges the one base that discriminates a
clonotype from its parent on its **Phred score** rather than on abundance, which is a measurement
the abundance model does not have.

**Never: Every denoising mode MOVES reads onto a parent and never discards them** — the sum of
`duplicate_count` is invariant across modes, and a clonotype with no qualifying parent keeps its
reads and is reported as an orphan. That is not caution: on a polyclonal hypermutated repertoire a
plain quality *filter* at the same threshold would strand **3.70 %** of all junction-bearing reads
with no parent to inherit them. If the read total moves when you change `--ec-mode`, that is a
defect — please report it.

⚠ The modes are **off by default** (`fast` = arda's historical behaviour), because whether the far
class they collapse is badly-sequenced SHM or error is not settled: on a hypermutated IGH
repertoire `amplicon` removes 178 clonotypes carrying 179 reads, 177 of them singletons, and on the
matched naive library it removes **zero**. Depth: [D segments](https://docs.isalgo.dev/arda/d_segments.html),
[SHM](https://docs.isalgo.dev/arda/shm.html),
[error correction](https://docs.isalgo.dev/arda/error_correction.html).

## CLI reference

The **stages** are separate commands, so any one can be run, inspected or replaced:

```bash
arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv   # Stage 1: receptor reads out of RNA-seq
arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --junction-quality  # + Phred over the junction
arda assemble -i mapped.airr.tsv -o assembled.airr.tsv         # Stage 3: CDR3s no single read spans
arda correct  -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv   # Stage 2: clonotypes
arda correct  -i mapped.airr.tsv -o clones.tsv --ec-mode accurate    # quality-gated correction
arda correct  -i mapped.airr.tsv -o clones.tsv --call-level gene     # collapse allele-level call splits
arda shm      -i mapped.airr.tsv -o rescoped.airr.tsv          # recount SHM outside the junction
arda stats    -i mapped.airr.tsv -c clones.tsv -r SAMPLE.arda.json -o SAMPLE.stats.tsv   # run QC
arda qc batch  -d results/ -o results/batch --samples sheet.tsv   # every sample's QC, one table
arda qc report -i results/batch.qc.json -o results/batch.qc.html  # ...as one self-contained page
```

Everything else:

```bash
arda info                                   # resolved paths + tool availability
arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa
arda annotate -i reads.fastq -o out.airr.tsv --strand forward   # plus-strand only
arda annotate -i reads.fastq -o out.airr.tsv --d-max-evalue 0.01  # the strict D band
arda markup -i junctions.tsv -o marked.tsv --report -           # mark up + repair bare (CDR3aa, V, J) records
arda resolve-ties -i mapped.airr.tsv -o widened.airr.tsv        # every germline the read cannot rule out
arda resolve-ties -i mapped.airr.tsv -o widened.airr.tsv --loci IGK,IGL   # ...on the loci where it helps
arda genotype -i mapped.airr.tsv -o donor.genotype.tsv --loci TRB          # which V alleles this donor carries
arda resolve-ties -i mapped.airr.tsv -o narrowed.airr.tsv --genotype donor.genotype.tsv  # ...and apply one
arda scenarios -i clones.tsv -o d_prior.tsv                     # EM over the recombination scenario set
arda shm-model -i mapped.airr.tsv -o shm.tsv --locus IGH        # where hypermutation is EXPECTED (5-mer context)
arda cluster submit --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE --shards 20 --partition cpu
arda igblast -i reads.fastq -o truth.airr.tsv                   # gold-standard IgBLAST (all loci)
arda export-ref --kind segments --locus TRB --format fasta      # the reference, out of the CLI
arda build-db   --organism all              # rebuild references (needs IgBLAST)
arda build-index --organism all             # (re)build the precompiled mmseqs DBs
```

### Quality control

Every run — `rnaseq`, `amplicon` and `cells` alike — writes **`<prefix>.stats.tsv`** and
**`<prefix>.stats.json`**: the numbers that decide whether a sample is usable, without re-reading
the FASTQ, and in the same scopes for every mode so one cohort table holds all of them.
`arda stats` produces the same table from any subset of a run's artifacts, so it also works on a
bare `arda annotate` output:

```bash
arda stats -i SAMPLE.airr.tsv -c SAMPLE.clones.tsv -r SAMPLE.arda.json \
           --r1 R1.fq.gz --r2 R2.fq.gz -o SAMPLE.stats.tsv
```

Four columns — `scope`, `key`, `metric`, `value` — one value per cell, so a metric can be grepped,
`join`ed across samples or plotted without reshaping:

| scope | key | what |
|---|---|---|
| `run` | `map` / `correct` / `assemble` | the run report verbatim: total and mapped reads, FASTQ bytes, read length, paired, threads, wall time, peak RSS |
| `sample` | — | library-wide totals, junction lengths and quality, SHM rate, V/J gene coverage |
| `chain` | `TRB`, `IGH`, … | per locus, **reads and clonotypes**: functional / non-functional, stop codons, truncated junctions, min/max/mean junction length, junction quality, SHM rate, chimeras |
| `v_gene` / `j_gene` | `TRBV19` | reads and clonotypes per germline gene |
| `allele_candidate` | `TRBV19*01:G45A` | a recurrent, high-quality V mutation, with its frequency and mean Phred |
| `junction_aa_len` · `read_len` · `clone_size` · `isotype` · `chain_support` | `IGH:17` | the distributions, keyed `locus:bucket` — junction length in residues, aligned read length in 10-nt buckets, clonotype size in powers of two, constant-region class, molecules per chain |

```
$ awk -F'\t' '$1=="chain" && $2=="IGH"' SAMPLE.stats.tsv
chain  IGH  reads                       104
chain  IGH  reads_truncated_junction    1
chain  IGH  junction_nt_mean            48.75
chain  IGH  junction_quality_mean       35.6718
chain  IGH  shm_rate                    0.0413
chain  IGH  clonotypes_chimeric         2
```

A metric with no input is **omitted, never reported as 0** — a run without `--junction-quality`
does not read as "mean quality 0", and a locus no read spanned has no minimum junction length
rather than one of zero. Two columns feed the quality metrics and both are opt-in on `map`:
`--junction-quality` (Phred over the junction bases) and `--mutation-quality` (the Phred behind
each `v_mutations` / `j_mutations` entry, one-for-one).

The distribution scopes are there because min/max/mean cannot show a **shape**, and shape is most
of the diagnosis: a bimodal junction length is two primer sets in one tube, a read-length cliff is
an adapter left on, and a clone-size distribution with no singletons is a library amplified before
it was sequenced. Each reads as an unremarkable mean.

#### A cohort, and a page

```bash
arda qc batch  -d results/ -o results/batch --samples sheet.tsv --report results/batch.qc.html
```

`qc batch` reads **only** the per-sample `stats.tsv` files — never an AIRR or a clonotype table —
so a thousand-sample cohort is one concatenation and runs on results copied off a cluster with the
bulk data left behind. It writes `batch.qc.tsv` (long, with each metric's group median, MAD and
robust z), `batch.qc.wide.tsv` (one row per sample) and `batch.qc.json`. The sample sheet's
optional `project` and `batch` columns name the group a sample is compared within; they are labels
and change nothing about a run.

⚠ **No threshold is shipped.** There is no pass/fail column and no constant saying what a good
`mapped_fraction` is — that depends on the library, the organism and the depth. What a batch can
say is whether a sample looks like the batch it came in with, so each metric carries its group's
median, MAD and a robust z (`|z| ≥ 3.5`, Iglewicz–Hoaglin); a group under 5 samples gets a median
and no z, because MAD over four points is not a scale. Repertoire biology — diversity, clonality,
rarefaction, overlap, cross-sample clonotype matching — is deliberately absent; that is
[vdjtools](https://github.com/antigenomics/vdjtools)', which takes arda as a dependency.

`arda qc report` renders any of it as **one HTML file with the data inlined and no external
reference of any kind** — no CDN, no bundle, no fetch — so it opens on an air-gapped login node or
as an email attachment long after the results directory is gone. Sortable sample table shaded by
z, any metric against its batch median, the distributions overlaid per sample, and per-sample
provenance. arda gains no plotting dependency for it.

⚠ **`allele_candidate` is a shortlist, not a genotype call.** A novel allele, somatic
hypermutation and a base miscall are the same string in the mutation list; what separates them is
how often the mutation recurs across an allele's reads and how good the base is. Both are reported
per variant and the thresholds are yours (`--allele-min-frac`, `--allele-min-reads`). arda does not
genotype. Likewise the chimera, non-functional and stop-codon counts are **flags, never filters** —
nothing in `stats` removes a row from any output.

#### Verbosity and logging

Logging is a stdlib logger configured by three **global** options, before the subcommand:

```bash
arda -v --log-file run.log amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/
arda -q rnaseq --r1 R1.fq.gz -p SAMPLE -d out/          # warnings and errors only
```

Default prints the stage lines and a throttled progress line (reads seen, reads mapped, reads/s,
RSS). `-v` adds DEBUG with the level and module name. `--log-file` is **always DEBUG whatever the
console level is**, and stamps every line with a timestamp and the process peak RSS — so a quiet
cluster job still leaves a full record.

`arda cluster` holds every sharded/SLURM helper. `arda cluster submit` shards Stage 1 across an
array and runs Stages 2–3 **once** over the merged output — byte-identical to a single-node run.
`arda cluster split-fasta` / `submit-fasta` are the **amplicon / single-end** FASTA path: they drop
quality and separate mates, so paired RNA-seq must not use them.

**`arda export-ref`** dumps arda's most valuable offline artifact — every in-frame V·J
scaffold with IgBLAST-quality FR1–4 / CDR1–3 coordinates — in 3 kinds × 4 formats:

```bash
arda export-ref --kind scaffolds --locus TRB --format gff3 -o trb.gff3   # V×J (and J+C) reference
arda export-ref --kind segments  --locus TRB --format fasta              # collapsed per-allele V/J/C
arda export-ref --kind anchors   --locus TRB                             # per-allele CDR3 anchors (tsv)
```

Coordinates are **1-based closed** (AIRR), so they pass through GFF3 unchanged; `--format airr`
shapes a scaffold as an AIRR Rearrangement row, feedable straight into anything reading arda's
own output. Details: [reference export](https://docs.isalgo.dev/arda/reference_export.html).

[`examples/`](examples/) is a runnable tour, every artifact derived from real data committed to
this repo and regenerated by `python examples/regenerate.py`: one real mRNA per locus; the two
human reads (of 7,341, across five organisms) that carry a **tandem D-D**; seven VDJdb records
covering every junction-repair outcome; and a 1,035-read FASTQ that runs the whole bulk RNA-seq
pipeline in ~6 s. See [`CHANGELOG.md`](CHANGELOG.md) for what changed per release.

Input may be FASTA or FASTQ, plain or gzipped. Nucleotide input is searched on **both strands**
by default (reverse-complement reads are re-oriented and flagged `rev_comp=T`); a single search
annotates a mixed bulk RNA-seq file across all loci.

## Widening the V call: `arda resolve-ties`, and the loci where it helps

A read aligned over `[germline_start, germline_end]` is explained exactly as well by any germline
carrying that same stretch, so naming one gene is a claim the data does not support.
`arda resolve-ties` widens `v_call`/`j_call` to every germline the alignment cannot rule out. It
adds no alignment — the tie is a string comparison against the reference over the span already
aligned — and it is **off by default**, because it changes `v_call` on every library.

⚠ **Whether it helps is a property of the locus.** Measured against an `arda igblast` truth
(`v_score >= 70`) on eleven arms across five geometries, scored on **exact `v_gene`-set
agreement** rather than intersection:

| locus | arms | exact before | exact after | change |
|---|---|---|---|---|
| **IGK** | 2 bulk | .5707 / .5558 | **.9373 / .9308** | **+36.7 / +37.5 points** |
| **IGL** | 2 bulk | .8586 / .8577 | **.9137 / .9115** | **+5.5 / +5.4 points** |
| **TRB** | 1 amplicon | .9299 | **.9694** | **+4.0 points** |
| IGH | 6 (bulk, 5'RACE, multiplex) | .8255–.9755 | .6860–.9608 | **−0.6 to −14.9 points** |

So `--loci IGK,IGL` widens only those loci and copies every other row through untouched: one
command on a mixed library gets the win where it exists and leaves IGH byte-identical.

The reason is mechanical rather than mysterious. The deciding quantity is the **ambiguity
deficit** — IgBLAST's genes per read minus arda's, before widening. **arda's IGK call is 0.42
genes/read too narrow** (1.18 against IgBLAST's 1.60), and closing that gap *is* the 37-point win.
**arda's IGH call is already as wide as IgBLAST's** (deficit 0.00–0.07 on every arm, 1.589 against
1.64 on a 5'RACE library), so there is nothing to recover and each arm only pays the overshoot.
**You can check this without a truth file**: compare mean genes per `v_call` before and after, and
if it moves far more than a few hundredths, the rule is overshooting on that locus.

⚠ The cost is that fewer reads name a single gene: on IGK the share falls **.8195 → .3883**. That
is the honest statement of what an IGK read supports, and it is why this is an opt-in command
rather than a default. Clonotype cost is +2 of 362 on IGK, +3 of 23,559 on TRB, and **zero** on
IGL. ⚠ Intersection recall rises on *all eleven* arms including IGH — only the exact-set number
shows the IGH regression.

## Personalized germline: `arda genotype`, `resolve-ties --genotype`

A reference catalogues every allele anyone carries; no donor carries all of them, and none carries
more than two per gene. Restricting V calls to the carried set removes ambiguity that was never
real — TIgGER measured 11.2 % → 1.5 % ambiguous assignments doing this on full-length BCR.

```bash
# apply an allele set you already trust -- one column of allele names is a valid genotype file
arda resolve-ties -i sample.airr.tsv -o sample.genotyped.tsv --genotype genotype.tsv

# or infer one from reads arda has already mapped
arda genotype -i sample.airr.tsv -o sample.genotype.tsv --loci TRB
```

Restriction adds `v_call_genotyped` and leaves `v_call` byte-identical. It **never re-aligns and
never rebuilds a reference**: given the span a read already aligned over, it is a set intersection.
An OGRDB set, a MiXCR-inferred library or a hand-written list are all valid input, and every named
allele is validated against the reference — an unknown one raises rather than silently restricting
nothing.

The inference groups reads by junction — **junctional diversity makes a junction de facto a UMI**,
so one junction is one rearrangement and one independent draw from the donor's two chromosomes, and
disagreement *within* a junction is error rather than allele (which is where the error rate is
measured from). The call is a **likelihood ratio between diploid genotypes**, not a coverage rule:
homozygous and heterozygous are the same multinomial formula at genotype size 1 and 2, and a gene
is called only if the winner beats the runner-up by `--min-log10-bf`.

⚠ **Allele-level genotyping is a read-length feature, and most libraries do not have it.**
Separating a gene's alleles needs a median of 150 nt of TRBV (175 TRAV, 230 IGHV) from the 3' end;
only 16 of 44 multi-allele human TRBV genes separate within 100 nt. On a real TRB amplicon
(`SRR5233641`, 151 nt paired, 33,440 clonotypes) reads cover a median of **72 nt** of V — 56 after
anchor-clipping, **none** reaching 150 — so 33.1 % of clonotypes can be assigned an allele and
**17 of 53 genes are called**: 11 trivially (one catalogued allele) and **6 genuinely inferred**,
at log₁₀ Bayes factors up to 251. The other 36 are refused, including one with 1,037 clonotypes
that still cannot separate its alleles. The refusals are the point. Details in
[`docs/genotype.rst`](docs/genotype.rst).

⚠ **The same axis sets what a `v_call` is worth on IG, one level coarser.** Naming the *gene* is
much cheaper than separating its alleles, but it is not free: measured against an IgBLAST truth
(`v_score >= 70`) on a human IGH 5'RACE library of 98,639 truth reads, `v_gene` recall is
**.1170 on reads covering under 60 nt of V germline and .9896 on reads covering 200 nt or more** —
and **.9896 is the TRA amplicon's .9867**, so there is no IG-specific accuracy deficit at that
coverage. The 5.9 % of reads under 60 nt produce about 56 % of every V miss.

**Position beats length**: IGHV genes diverge in FR1/CDR1/CDR2 and are conserved through FR3 near
Cys104, so a short **5'** alignment identifies the gene and a short **3'** one does not — the same
`< 60 nt` bin scores **.1170** on 5'RACE (whose reads run C → J → V) against **.8472** on bulk
RNA-seq of the same locus. Somatic hypermutation does *not* order this: stratified by IgBLAST's
own `v_identity` the deficit is non-monotonic, with the nearly-unmutated `97–99 %` bin the worst
of six. arda ships **no threshold** on it — a span gate was measured at 7.6 : 1 in favour on IGH
5'RACE and a clear loss on bulk IGH/IGK/IGL — so `v_germline_start` / `v_germline_end` are there
for you to filter on. Tables in
[`docs/usage.rst`](docs/usage.rst#what-an-ig-v-call-means-when-the-read-is-short).

## Productivity: `productive`, `stop_codon`, `vj_in_frame`

These three AIRR columns are the most misread ones arda writes, because each is scoped
differently and none of them means "the junction looks like germline".

| column | what arda puts in it |
|---|---|
| `vj_in_frame` | `T` when the V's frame, carried through the junction, arrives at FR4 on a codon boundary — so FR4 **and the constant exon spliced to it** are translated in their own frame. `phase = (fwr4_start - v_coding_start) % 3`, `T` iff `phase == 0`. |
| `stop_codon` | `T` when any annotated region — FR1–FR3, CDR1–CDR2, the junction, **or FR4** — carries a `*`. Scoped to the annotated span, not the whole read. |
| `productive` | The conjunction: in frame **and** no annotated region carries a stop. |
| all three | **Empty** when the read never reached both the junction and FR4 — *unevaluable*, not non-productive. On real bulk RNA-seq that is ~72 % of mapped reads, and writing `F` there would look like a repertoire of broken rearrangements. |

They are **flags, never filters**: no stage drops a non-productive read. One allele of a T cell's
two is usually out of frame, so a sample with none has been filtered upstream.

`vj_in_frame` says **where the frame lands at FR4**, which is what decides whether C translates.
It does not say the junction's 3′ residues match the called J's germline residues. Those are
different tests, and on `JX277388.1` (mouse TRB, `TRBJ2-5*01`) they disagree — the read carries
two nucleotides `TRBJ2-5*01` does not have, so the J's 5′ germline residues and its own FR4
cannot both be in the germline frame:

```
arda:  CTCSAGGPRHPVLF GPGTR          FR4 reproduces germline GPGTR exactly
       + TRBC2*01  ->  ...FGPGTRLLVL EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARG...   0 stop codons
```

The receptor translates, so `vj_in_frame=T` is right. Matching the junction's tail against
germline instead is not a frame test at all: it fires on **51.8 % of IGL, 51.4 % of IGK and
39.8 % of IGH** rows (mean `v_identity` .955–.968) against **3.4 % of TRB and 7.1 % of TRA**
(.992) — it is detecting hypermutation reaching the J, which is biology.

Worked examples, including a J truncated too far for any tool to call:
[productivity](https://docs.isalgo.dev/arda/productivity.html).

## Why

IgBLAST is the gold standard but is slow to invoke per-batch and awkward to embed.
`arda` keeps IgBLAST-quality region calls while being:

- **Fast & scalable** — MMseqs2 search + a C++ projection step; on a TRA amplicon it
  matches MiXCR's wall clock at **3.6× less CPU and 4.8× less RSS**
  ([below](#performance)); multiprocessing and SLURM-friendly from small FASTA to large FASTQ.
- **Embeddable** — `import arda; arda.annotate_sequences(...)`.
- **Honest** — a D call is gated on an E-value that ships with it (`d_support`), a
  germline allele with no derivable anchor is flagged rather than guessed, a junction
  whose Cys104 anchor is not actually in the read is **not emitted**, and a `j_call`
  requires J evidence rather than being inherited from the scaffold.
- **Easy to install** — `pip install arda-mapper` (binary wheels ship the C++
  extension); the `mmseqs` binary is fetched as a static build at runtime — no
  conda. IgBLAST is fetched on first use the same way, and is only needed to
  (re)build the reference DB or to run `arda igblast`, never for annotation.

## Performance

### Head-to-head: the full pipeline, to a clonotype table

⛔ **A stage comparison is not a benchmark.** Every leg below runs **end to end and emits
clonotypes**; the only latitude a tool gets is picking its best-fitting preset for the library
type. One job, six legs alternating, three reps, medians reported. M3 Mac, 8 threads, arda 2.27.0
· MiXCR 4.7.0 · TRUST4. TRUST4 rows count only its **complete** CDR3s (`C…[FW]`, no `_`/`?`), so
the clonotype and read columns mean the same thing in all three rows.

**TRA amplicon, 100,000 reads** — `arda amplicon` · `mixcr analyze generic-amplicon --rna` ·
TRUST4 defaults:

| pipeline | wall (s) | CPU (s) | peak RSS (MB) | clonotypes | reads in clonotypes |
|---|---:|---:|---:|---:|---:|
| MiXCR `generic-amplicon` | **7.82** | 42.91 | 3,052 | 19,697 | 42,712 |
| **arda `amplicon`** | 14.40 | **30.49** | 965 | **19,841** | **43,503** |
| TRUST4 | 73.77 | 117.96 | **490** | 18,559 | 37,688 |

**MiXCR is 1.84× faster on wall** on a primer-anchored amplicon — that is its regime and the
number is not disputed here. arda spends **1.41× less CPU** and **3.16× less RSS** to get there,
and returns the most clonotypes (+0.7 % over MiXCR, +6.9 % over TRUST4) over the most reads
(+1.9 %, +15.4 %). TRUST4 is **5.1× slower than arda** on this arm.

**Bulk RNA-seq, 660,000 pairs (SRR5233639)** — `arda rnaseq` · `mixcr analyze rna-seq` · TRUST4
defaults:

| pipeline | wall (s) | CPU (s) | peak RSS (MB) | clonotypes | reads in clonotypes |
|---|---:|---:|---:|---:|---:|
| TRUST4 | **13.68** | **54.53** | **460** | 1,941 | 6,247 |
| **arda `rnaseq`** | 21.95 | 219.68 | 1,033 | **2,213** | **8,484** |
| MiXCR `rna-seq` | 30.23 | 233.11 | 2,849 | 1,732 | 4,288 |

**arda finds the most on the regime it exists for**: +27.8 % clonotypes and **+97.9 % reads
assigned** against MiXCR, +14.0 % and +35.8 % against TRUST4, at **1.38× MiXCR's wall** and
**2.76× less RSS** for comparable CPU. ⚠ TRUST4 is genuinely **1.61× faster on wall at 4.0× less
CPU** on this arm, and it is the cheapest of the three on memory in both — it reaches 87.7 % of
arda's clonotypes and 73.6 % of its assigned reads to do it. All three walls were stable across
three reps (arda 21.88–24.04, MiXCR 30.09–30.83, TRUST4 13.65–14.22).

**No regression across eight releases.** arda **2.18.0** (the version these arms were last
measured on) against this one, same job, legs alternating, 3 reps, the same committed reference
and the same mmseqs binary — `database/` has not changed since `v2.18.0`, which is what makes the
comparison valid:

| arm | 2.18.0 wall (s) | 2.27.0 wall (s) | 2.18.0 RSS | 2.27.0 RSS | clonotypes / reads |
|---|---:|---:|---:|---:|---|
| amplicon 100 k | 13.94 | **13.66** | 939 MB | 939 MB | 19,841 / 43,503 |
| bulk 660 k pairs | **21.45** | 21.74 | 1,033 MB | 1,033 MB | 2,213 / 8,484 |

Medians of three. Bulk is **1.4 % slower on wall and 0.8 % more CPU**, inside the rep spread
(2.18.0 21.06–22.26, 2.27.0 21.52–22.07) and accounted for by the per-run QC stage that 2.20.0
added; amplicon is 2.0 % faster. ⛔ **A count-equal leg can still be a changed leg**, so the check
is a call digest over `(locus, v_call, j_call, junction, duplicate_count)`, not a row count: both
arms are **byte-identical** between the two versions.

**IGH RepSeq amplicon, 100,000 pairs, 32 threads** (aldan3). What the amplicon configuration is
worth against the shipped one-pass default on hypermutated IGH:

| dataset | config | wall (s) | peak RSS (MB) |
|---|---|---:|---:|
| IGH_repertoire | one-pass default | 316.44 | 4,018 |
| IGH_repertoire | `--two-pass --fast-segments --v-only-on-segment` | **76.25** | **1,479** |
| IGH_naive | one-pass default | 305.32 | 3,736 |
| IGH_naive | `--two-pass --fast-segments --v-only-on-segment` | **64.86** | **1,363** |

4.15× and 4.71×, at ~2.7× less memory.

**vs TRUST4 on IGH amplicon at 32 threads** (round 20, aldan3; every leg of a tier ran on the same
staged input, so no ratio here is cross-job):

| dataset | reads | arda `amplicon` wall (s) | TRUST4 wall (s) |
|---|---:|---:|---:|
| IGH_repertoire | 100,000 | **201.05** | 615.04 |
| IGH_naive | 100,000 | **136.51** | 359.66 |
| migec_exp1_TCR | 500,000 | **316.25** | 423.96 |
| migec_exp1_IGH | 500,000 | 223.77 | **225.10** |

⚠ Wall clock only, and the full-depth (hours-scale) head-to-head is still cluster work.


### Synthetic benchmarks vs IgBLAST

⚠ The two tables below are **synthetic** — generated human IGH sequences, not a real library —
from [`scripts/bench_vs_igblast.py`](scripts/bench_vs_igblast.py) and
[`scripts/bench_prefilter.py`](scripts/bench_prefilter.py). They measure scaling shape, not
head-to-head standing; use the tables above for that. 16 threads.

| sequences | arda wall | arda rate | speedup vs IgBLAST | region concordance |
|----------:|----------:|----------:|-------------------:|-------------------:|
| 10,000    | 5.5 s     | ~1.8k/s   | 4.4×               | 98.9%              |
| 50,000    | 16 s      | ~3.0k/s   | 7.3×               |                    |
| 100,000   | 30 s      | ~3.3k/s   | 7.9×               |                    |

Bulk RNA-seq is faster per read than amplicon, because mmseqs prefilters by k-mer matching —
reads with no receptor k-mer are rejected before alignment. 150 nt reads, 16 threads:

| receptor content | throughput |
|-----------------:|-----------:|
| 100% (amplicon)  | ~5.7k reads/s |
| 10%              | ~19k reads/s |
| 1% (blood RNA-seq) | ~25k reads/s |

### Memory

arda is **CPU-bound**; large FASTQ is streamed in bounded chunks (a background reader prefetches
the next chunk while the current one is annotated), so **mapping is flat at ~300–650 MB at any
read depth**. Peak RSS tracks **repertoire richness**, not reads: Stage 3 holds the clone set, so
on a B-cell-rich tumour (28,444 clonotypes from 105 M reads) `correct` peaked at **2,071.7 MB**,
while a colder sample with *more* reads (139 M) peaked at **549 MB**. **Budget ~4 GB**, and size a
SLURM `--mem` from Stage 3, never Stage 1.

## Accuracy

### Recall on bulk RNA-seq — arda vs TRUST4 vs MiXCR

The three-way comparison, on the **same 16 datasets where every tool ran**, 5,273 real fragments,
Wilson 95 % CIs. Best bold; ties bold together.

| tool | config | recall | precision (lower bound) | false positives | FP / 1M reads | peak RSS |
|---|---|---|---|---:|---:|---:|
| **arda** | shipped defaults | **0.986** [.982–.989] | **0.889** [.881–.897] | **70** | **21.9** | **254 MB** |
| **TRUST4** | defaults | **0.987** [.984–.990] | 0.160 [.156–.164] | 22,185 | 6,932.8 | 306 MB |
| MiXCR | `-OallowNoCDR3PartAlignments=true -OminSumScore=40` | 0.964 [.958–.969] | 0.051 [.050–.052] | 91,229 | 28,509.1 | 1,213 MB |
| MiXCR | `align --preset rna-seq` as shipped | 0.193 [.183–.204] | 0.565 [.542–.588] | 91 | 28.4 | 1,213 MB |
| arda | `--reconstruct` (6 paired datasets) | 0.995 | 0.905 | 14 | 11.7 | 249 MB |

**Recall is a statistical tie between arda and TRUST4** — a 6-fragment gap with overlapping CIs, so
both are marked winners and neither should be quoted as "higher recall" than the other.
**Precision is not a tie**: arda's lower bound (0.881) sits above every competitor's *upper* bound.

⚠ **Benchmark MiXCR at its best config, not its default.** `align --preset rna-seq` as shipped
gives 0.193 recall; the two free options above take it to 0.964. Quoting the default alone is the
mistake this project made and had to retract.

⚠ Precision here is a **lower bound** — the grey band `30 ≤ v_score < 70` is scored under neither
metric. Never: And TRUST4's recall and FP are scored on its **candidate extraction**, a different (and
much less filtered) stage than arda's post-`--min-score` output; the two are not like-for-like on
FP, which is why the FP column carries a per-1M normalisation rather than a bare ratio.

What explains the table is the **J→C class — 22 % of real fragments** on bulk:

| tool | V-covered reads (n = 4,117) | J→C agreement (n = 1,156) |
|---|---:|---:|
| arda | 0.9864 | **0.9844** |
| TRUST4 | 0.9944 | 0.9611 |
| MiXCR (recall config) | **0.9990** | 0.8382 |
| MiXCR (default) | 0.1482 | 0.3538 |

On V-covered reads every tool is ≥ 0.986. arda's overall recall rests on the 345 `J + C` scaffolds
in its reference, which took that class from 0.0606 to 0.9844 and overall recall from 0.78 to 0.986.
⚠ Reported as *agreement*, not recall: arda now ships C scaffolds, so the adjudicator is no longer
independent of it there.

### Gene calls on a targeted amplicon

Against an IgBLAST truth on the TRA amplicon, 100,000 reads, arda 2.27.0 and MiXCR 4.7.0 at its
best amplicon preset, both scored per read from the same truth file in the same job.

⛔ **Print the coverage before the rates.** A per-tool inner join hands each tool its own
denominator — a truth read the tool emitted no row for simply vanishes instead of counting as a
miss. Over the 48,033 truth reads at `v_score ≥ 70`, **arda emits a row for 48,030 (99.99 %) and
MiXCR for 46,503 (96.81 %)**. So both denominators are shown:

| metric | arda 2.27.0 | MiXCR 4.7.0 |
|---|---:|---:|
| **all 48,033 truth reads** | | |
| `v_gene` recall | **.9867** | .9660 |
| `v_gene` precision | **.9996** | .9977 |
| `j_gene` recall | **.9892** | **.9892** |
| `j_gene` precision | .9953 | **.9996** |
| `junction` recall (nt, exact) | .9473 | **.9708** |
| **common subset, 46,502 reads** | | |
| `v_gene` recall | .9869 | **.9977** |
| `v_gene` precision | **.9997** | .9977 |
| `j_gene` recall | .9959 | **.9996** |
| `junction` recall (nt, exact) | .9533 | **.9778** |

The two views say different and equally true things. **On the reads it emits, MiXCR is the more
accurate caller.** **Over the whole library, arda recalls more V genes** — .9867 against .9660 —
because MiXCR emits nothing at all for 1,530 truth reads and arda for 3. arda's V calls are also
the more **precise** of the two under either denominator: it declines rather than guessing.

Junction, stated the same way: of the 46,787 truth junctions, **arda emits one for 94.81 % at
.99919 precision among emitted; MiXCR for 97.09 % at .99989.** The 5.19 % arda declines are reads
that *have* an anchor pair and lost the projection — a known and specific gap, not a calling
error. MiXCR suffixes every allele `*00`, i.e. makes no allele call at all, so there is no
`v_allele` comparator; arda's is **.9868 resolved** (.9461 by exact string, the 4-point difference
being ambiguous-allele tie lists, which are a scoring convention and not a call).

⚠ These five arda figures are **unchanged from 2.11.1**, fifteen releases back: `v_gene` recall
.9867, precision .9996, `j_gene` recall .9892, precision .9953, junction precision among emitted
.99919 all reproduce to every digit published.

**Score alleles as tie lists, not exact strings.** IgBLAST and arda both return an ambiguous
allele as a comma-joined set; scoring that as a miss is a scoring artifact, not an error. Across
25 datasets the median is `v_allele_exact` **.8328** against `v_allele_resolved` **.9763**.

**A V/J boundary disagreement *inside* the junction is not an error.** V(D)J recombination is
probabilistic: exonuclease chew-back and N/P-nucleotide addition mean the V-end / N-D-N / J-start
partition is often not identifiable from sequence alone, so overlapping V/J/NDN assignments are
acceptable and the ground truth is unknown. What *is* checkable, and what the table scores: the
junction's **outer bounds** (Cys104 and [FW]118), the **gene/allele calls**, and whether a tool
invents a junction it has no anchor for.

On ~7.3k real GenBank mRNA records spanning **all five organisms and their loci** (committed,
gzipped test fixtures), region concordance with IgBLAST on productive records is **98–99.7%** per
organism, and `junction_aa`/`cdr3_aa` match IgBLAST ~99% while satisfying the AIRR invariants
exactly. (GenBank also contains genomic/partial/non-productive entries that confuse both tools;
those are excluded.)


## The three stages, in detail

The modes run these for you; run them separately when you want to inspect, replace, or shard one
of them. `arda rnaseq` is a recall-first pipeline for extracting the repertoire from bulk RNA-seq,
where 1–5% of reads are receptor-derived:

```bash
arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --report run.json
arda assemble -i mapped.airr.tsv -o assembled.airr.tsv
arda correct  -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv
```

`--extra-airr` is what folds Stage 3 back in; **without it the assembled reads are silently
discarded.** `arda rnaseq` / `arda amplicon` do all three in one call and wire it for you.

- **`map`** streams paired FASTQ, keeps only reads mapping to a receptor scaffold, and writes
  them as AIRR. The reference includes **`J + C` constant-region scaffolds**, so a read spanning
  the J→C splice — which ends in the constant region and has no V to anchor — still maps and
  carries a **`c_call`** (the CH1 exon) plus a **`c_class`** isotype (`IGHG`/`IGHM`/`IGHA` … —
  the class, never the noise-prone subclass). In paired mode the isotype of a CDR3-bearing read
  is recovered from its constant-region mate. `--reconstruct` merges each overlapping mate pair
  into one fragment, resolving overlap mismatches by the higher-Phred base.
- **`assemble`** reconstructs clonotypes whose CDR3 is too long for any single 100–150 bp read to
  span (V(DD)J ultralong, ~20–40 aa) by greedy overlap-extension anchored on Stage-1's per-read
  `cdr3_start`.
- **`correct`** aggregates reads into clonotypes and collapses sequencing-error CDR3 variants.
  Abundance is the AIRR **`duplicate_count`** — every read that *encompasses* the junction,
  the true expression estimate — with **`consensus_count`** for distinct fragments.

`arda igblast -i reads.fastq -o truth.airr.tsv` runs IgBLAST across all loci as a gold-standard
reference for benchmarking (see the `arda-benchmark` project). `--receptor ig` or `tr` skips a
whole pass over every read when the library's receptor type is known; the default `both` runs
both and keeps whichever scores higher.

### Error correction: `--error-rate` is a per-library calibration

`correct`'s error model is per-base with a **length-scaled threshold** (a mismatch over a longer
junction is likelier an error) and is SHM-indel-tolerant. Its one knob is `--error-rate`, and the
default (1e-3, ~Phred 30) is **not** universally safe:

On the MIGEC spike-in library (PRJNA239303), at the default `--error-rate 1e-3` arda erases both
published spike-in variants. That is a signal-to-noise limit, not a defect: on the paper's own
metric computed over raw reads, **V1/Err1 = 1.35 and V2/Err2 = 0.28** — the second variant is
*less* abundant than the worst 2-substitution PCR error in the same library, so **no
abundance-based method, at any threshold, can separate them.** That is precisely why UMI
consensus exists; it moves V1/Err1 to 26–76 before any correction runs.

What to do about it:

```bash
arda correct -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv --error-rate 1e-5
```

`1e-5` recovers both spike-in variants exactly. On an independent error cloud, `1e-4` kept both
while still removing 72% of the real PCR errors. `--error-rate` has a physical reading — **~1e-3
for raw reads** (the sequencer's Phred-30 substitution rate) and **~1e-5 for UMI-consensus input**
— so calibrate per library rather than trusting one default.

**And it no longer has to cost precision.** `1e-5` used to buy the variants at 3.5 points of
purity on a monoclonal control (.99540 → .96034), because it also stops collapsing real PCR error.
`correct --min-junction-q` gates on the Phred score of the one base that discriminates a clonotype
from its putative parent — a measurement abundance does not have, and one that separates cleanly
(the published spike-ins read median Q 34–35 there, the error cloud around them **median Q 24**):

```bash
arda map     --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv --junction-quality
arda correct -i mapped.airr.tsv -o clones.tsv --error-rate 1e-5 --ec-mode accurate
```

Both variants kept, monoclonal purity back to **.99530**, spurious junctions 297 → 62, distinct
error clonotypes 1,630 → 124. What it cannot do is rescue a variant below the **template-error
floor**: RT and early-PCR errors happen *before* the UMI is attached, so consensus cannot remove
them and they are high-Q by construction. Full write-up, including why `binom`/`betabinom` are not
`accurate`: [error correction](https://docs.isalgo.dev/arda/error_correction.html).

## Library

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],  # strings or (id, seq) pairs
    seqtype="nt", organism="human",
)
# -> list of AIRR record dicts: v_call, d_call/d2_call, j_call, c_call/c_class,
#    fwr1..fwr4, cdr1..cdr3, *_start/*_end (1-based closed), *_aa, junction(_aa),
#    np1/np2/np3, d_support/d2_support, {v,j,c,d}_cigar, v_mutations/j_mutations,
#    sequence_alignment, germline_alignment, productive, ...
# The TSV is a spec-valid AIRR Rearrangement file (passes airr.schema validation).
```

### Analysing the output

Every arda TSV is written with quoting **off**, so read it with quoting off — polars' default
turns a blank field into the two-character value `""`:

```python
import polars as pl

clones = pl.read_csv("results/PT01.clones.tsv", separator="\t", quote_char=None)
top = (clones.filter(pl.col("locus") == "TRB")
       .with_columns((pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias("frequency"))
       .sort("duplicate_count", descending=True)
       .select("junction_aa", "v_call", "j_call", "duplicate_count", "frequency"))
```

Clonal fractions per locus, V-gene usage, repertoire overlap, SHM per read, gating a run on
`<prefix>.arda.json` and cohort QC from the `stats.tsv` files, all runnable:
[recipes](https://docs.isalgo.dev/arda/examples.html).

### Records with no read behind them

A VDJdb row is a CDR3 amino acid, a V call and a J call. There is nothing to align, but the
germlines still template a known run of residues into each end of the junction:

```python
from arda.cdr3fix import markup_cdr3
from arda.annotate.dmap import map_d_junction
from arda.dpost import posterior_d

mk = markup_cdr3("CAIRDDKII", "TRAV12-3*01", "TRAJ30*01", "HomoSapiens")
mk.cdr3_repaired            # 'CAIRDDKIIF'  -- the Phe118 anchor restored
[str(e) for e in mk.errors] # ["J del@8 missing 'F' d=0"]
mk.good                     # True: both sides repaired, both anchors present

junction_nt = "TGTGCTCTTGGGCCCCGGCCTTCCTACAGCGAGGAGTTGGGGGATACCCATCGGGCCGATAAACTCATCTTT"
map_d_junction(junction_nt, "TRDV1*01", "TRDJ1*01", "human").d2_call      # 'TRDD3*01' (tandem D-D)
posterior_d("CASSPLGQAYEQYF", "TRBV5-1*01", "TRBJ2-7*01", "human").d_call  # 'TRBD1'
```

**Coordinates here are *junction* space** — Cys104 through Phe/Trp118, **both anchors
included**. That is what VDJdb's `cdr3` column holds, and it is *not* arda's `cdr3` field,
which excludes both. Conflating them silently corrupts every coordinate.

Repair is conservative: only edits adjacent to a conserved anchor are applied, everything
deeper is *reported* and left alone, and a repair is refused outright unless the result opens
with Cys104 and closes with Phe/Trp118.

### Annotating bare germline segments

There is no coverage filter, so a **V-only** or **J-only** query maps to its scaffold and only
the regions inside the query's coverage are returned — a bare V yields `fwr1..fwr3`, a bare J
yields `fwr4`:

```python
from arda.annotate.mapper import annotate_records

recs = annotate_records(
    [("TRBV9*01", v_germline_nt), ("TRBJ2-7*01", j_germline_nt)],
    organism="human", seqtype="nt", strand="forward", map_d=False,
)
```

(mirpy uses exactly this to bake per-allele FR/CDR subsequences into its gene library; see
`tests/synthetic/test_germline_segments.py`. `arda export-ref` is the CLI equivalent.)

## How it works

1. **Reference build** (`arda.refbuild`, offline): download IMGT/V-QUEST germlines
   → enumerate deduplicated in-frame **V×J** scaffolds (D only affects CDR3
   interior, so it isn't enumerated) plus **`J + C` constant-region scaffolds** (the
   CH1 exon spliced onto each J, so J→C reads have somewhere to land) → annotate with
   `igblastn -outfmt 19` → translate → write `database/vdj/<organism>/{alleles.fasta,
   alleles.aa.fasta, markup.tsv, markup.aa.tsv, combinations.tsv, d_germlines.fasta,
   cdr3_anchors.tsv, d_prior.tsv, build.log}`.
2. **Runtime** (`arda.annotate`): MMseqs2 search query→scaffolds → best hit →
   C++ `transfer_regions` projects scaffold region coordinates onto the query
   (handling indels, truncation, mid-codon alignment starts, reverse strand) → for
   VDJ loci a gapless C++ local alignment of the CDR3 interior against the D
   germlines adds `d_call`/`d2_call` + `np*`; a hit on a `J + C` scaffold adds
   `c_call`/`c_class` → AIRR TSV. Ambiguous D and C calls are comma-joined allele
   lists, as V/J already are. Out-of-frame junctions are reported with an N-bridge
   (`_`) so FR4 still reads.

   The V..J interior is bounded by the **per-allele junction anchors** in
   `cdr3_anchors.tsv`, not by the scaffold projection — a scaffold has a 9 nt N-pad
   where a read has a 20–40 nt N-D-N region, so the projection collapses the very window
   the D lives in. A junction is emitted only when the read actually reaches its Cys104
   anchor. The D call is accepted on a Karlin–Altschul E-value (`d_support`, re-thresholdable
   with `--d-max-evalue`) rather than a per-locus score floor, and is constrained by germline
   geometry **before** the statistics: `/OR` orphons sit outside the locus and cannot rearrange
   at all; TRBD2 lies 3′ of the entire TRBJ1 cluster, so a TRBJ1 rearrangement can never be
   assigned TRBD2; and a tandem D-D must run in genomic 5′→3′ order, since deletional joining
   cannot produce any other. D mapping also runs on `--seqtype aa`, against each D germline's
   three translated frames.
3. **Bare records** (`arda.cdr3fix`, `arda.dpost`): a VDJdb-style row — CDR3 amino acid,
   V, J, species, and no read — is marked up against the same anchors (`arda markup`),
   its errors located and conservatively repaired, and optionally given a D gene inferred
   from the junction *length* (`--d-posterior`).

Fast sequence primitives (`translate`, `detect_coding_frame`, `reverse_complement`,
`back_translate`) live in the C++ extension and are re-exported from
`arda.refbuild.translate` — mirpy-API-compatible, so mirpy can `import arda` and reuse them.

The reference ships with **precompiled MMseqs2 indexes** (`database/vdj/<organism>/mmseqs/`),
used automatically when the local MMseqs2 version matches; otherwise arda transparently
rebuilds a private cache on first run (`arda build-index` regenerates the shipped DBs).
`segments.fasta` — the collapsed per-allele reference the two-pass path uses — is *generated*,
not shipped, and is built on demand when missing.

## Pipeline integration

`arda rnaseq` / `arda amplicon` write `<prefix>.clones.tsv` (AIRR clonotypes), `<prefix>.airr.tsv` (mapped
reads), `<prefix>.assembled.airr.tsv` and `<prefix>.arda.json` (run report). Because it is a
plain CLI over named files, it drops into any workflow engine with no glue code.

A ready-to-use **Nextflow module** lives in [`integrations/nextflow/arda/`](integrations/nextflow/arda/):
copy it to `modules/local/arda/`, feed it the trimmed per-sample FASTQ channel the aligners
already use, and it publishes per-sample clonotype tables to `${params.outdir}/arda/`. It ships a
conda `environment.yml` (works with `-profile conda`) and a `Dockerfile`, and emits a
`versions.yml`. See its [README](integrations/nextflow/arda/README.md) and the
[pipeline-integration guide](https://docs.isalgo.dev/arda/pipeline_integration.html).

The module exposes the regime and the denoising framework as params, so nothing needs
`task.ext.args` surgery:

```groovy
params {
    regime             = 'amplicon'   // or 'bulk' (default); per-sample via meta.regime
    arda_ec_mode       = 'amplicon'   // fast (default) | accurate | amplicon | rnaseq
    arda_clonotype_key = 'junction'   // full (default) | junction
    arda_mmseqs        = '/opt/conda/bin/mmseqs'
}
```

It validates both against their allowed sets and **warns** when `arda_ec_mode` and `regime`
disagree (e.g. the amplicon preset on a bulk library), because those presets are tuned to opposite
clonotype-size distributions and picking the wrong one is a real cost rather than a no-op.

**On a cluster**, `arda slurm` writes a `submit.sh` chaining split → `sbatch --array` →
merge with an `afterok` dependency, and a sharded run is **byte-identical** to a single-node one.
Never: Shard Stage 1 only: error correction compares a clonotype against its neighbours by abundance,
so running it per shard asks the question against a fraction of the evidence. See the
[cluster guide](https://docs.isalgo.dev/arda/cluster.html).

### A sample split across several FASTQs

Lanes or chunks are **read groups** of one sample, not separate samples. `--r1`, `--r2` and `--id`
are repeatable and matched by position (repeat an id to merge), or pass an nf-core-shaped
`--samples sheet.tsv`. The result is byte-identical to the same reads concatenated, and the unit
of parallel work becomes the read group — `arda cluster submit-samples` and `arda cluster plan`
schedule one task each.

```bash
arda rnaseq --samples sheet.tsv -d out/                          # sample,fastq_1,fastq_2
arda cluster submit-samples --samples sheet.tsv -d out/ --submit # SLURM, one task per read group
```

Full guide: [samples split across files](https://docs.isalgo.dev/arda/samples.html).
## Roadmap / TODO

See [`ROADMAP.md`](ROADMAP.md) for the full list. Done: the V·J reference build across
5 organisms, MMseqs2 mapping with C++ markup transfer, all-loci querying, streaming I/O,
**D-segment mapping incl. D-D fusions** (nt *and* aa input, E-value gated, genomic-order
constrained), **constant-region `J + C` scaffolds** (`c_call`/`c_class` isotype), **bulk
RNA-seq mode** with long-CDR3 contig assembly and coverage-based expression, **junction
markup + repair on bare records**, **multi-node (SLURM) sharding**, **the segment-based
fast paths** (`--two-pass`, `--fast-segments`, `--v-only-on-segment`, `--prefilter`,
`--indel-rescue`), **per-segment SHM in germline coordinates** (`v_mutations`/`j_mutations`),
**quality-gated error correction** (`--junction-quality`, `--min-junction-q`, `--ec-mode`)
and **`arda export-ref`**.

Next: full-depth clonotype benchmarking; a segment-native per-read path that removes MMseqs2
from the amplicon hot loop; and an `arda.hmm` semi-Markov model of V→N→D→N→J that would
subsume the E-value gate, the genomic-order constraint and the D posterior into one
forward–backward pass.

## Development

```bash
pip install -e .                                      # rebuilds the C++ ext on import
python -m pytest tests/unit tests/synthetic -q        # fast suite (needs mmseqs)
python -m pytest tests/realworld -q                   # vs IgBLAST, on committed fixtures — offline
env RUN_BENCHMARK=1 python -m pytest tests/benchmark -s   # timing/memory/scaling
```

Optional extras gate optional suites: `.[groundtruth]` (`olga`) for the generative ground-truth
tests that keep `arda.cdr3fix` honest, `.[test]` for `airr` schema validation. Without them those
tests **skip**, so `pip install -e '.[test]'` before reading a green suite as full coverage.

Layout: `src/arda/{refbuild,annotate}`, C++ in `src/_markup/markup.cpp` and
`src/_segmap/segmap.cpp`, references in `database/`, downloads in gitignored `bin/` + `data/`.
