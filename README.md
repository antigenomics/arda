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

<p align="center">
  <a href="https://docs.isalgo.dev/arda/quickstart.html"><b>Quickstart</b></a> ·
  <a href="https://docs.isalgo.dev/arda/"><b>Documentation</b></a> ·
  <a href="https://docs.isalgo.dev/arda/cli.html"><b>Commands</b></a> ·
  <a href="https://docs.isalgo.dev/arda/outputs.html"><b>Output reference</b></a> ·
  <a href="https://docs.isalgo.dev/arda/benchmarks.html"><b>Benchmarks</b></a>
</p>

---

**arda annotates T and B cell receptor sequences.** Give it sequencing reads — bulk RNA-seq,
targeted amplicon or single cell — and it returns [AIRR
Rearrangement](https://docs.airr-community.org/) records with V/D/J/C gene calls, FR1–FR4 and
CDR1–CDR3 coordinates and the junction, plus a clonotype table. One command, no workflow engine,
no container.

It works on nucleotide and amino-acid input, on FASTA and FASTQ, across all loci in one pass, for
**human** and **mouse** (full IG + TR) and for **rat, rabbit and rhesus macaque** (IG only — IMGT
ships no TR reference for those three). It also annotates records with no read behind them — a
CDR3 amino acid plus a V and J call, as in a VDJdb row — marking up which residues each germline
templates, repairing the junction, and inferring the D gene from junction length.

## Install

```bash
pip install arda-mapper
```

That is all. The distribution is `arda-mapper` and imports as `arda`; binary wheels carry its C++
extensions for CPython 3.10–3.13 on Linux, macOS (arm64) and Windows. The MMseqs2 search binary and
the curated germline reference are fetched automatically on first use into `~/.cache/arda` — no
conda, no manual reference build. Run `arda info` to see what got resolved.

For a development checkout, `bash setup.sh` builds a uv venv, compiles the extensions and fetches
IgBLAST and MMseqs2. See the [installation guide](https://docs.isalgo.dev/arda/installation.html).

## Quick start

Pick the mode that matches your library. The mode name carries that library's speed configuration,
so there are no flags to get right:

```bash
arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # bulk / whole-transcriptome
arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # targeted RepSeq / 5'RACE
arda cells    asm/PBMC.consensus.fq.gz -p out/PBMC            # single cell, from UMI consensus
```

Each mode runs mapping → assembly → error correction and writes:

| file | contents |
|---|---|
| `SAMPLE.clones.tsv` | one row per clonotype — the table most analyses start from |
| `SAMPLE.airr.tsv` | one AIRR Rearrangement row per mapped read |
| `SAMPLE.arda.json` | run report: reads in, reads mapped, why the rest were not, wall time, peak RSS |
| `SAMPLE.stats.tsv` | run QC in long format, for per-sample and cohort checks |

Progress goes to stderr; the output paths, one per line, to stdout.

Already have assembled sequences? Skip the pipeline:

```bash
arda annotate -i sequences.fasta -o out.airr.tsv --organism human
```

Full walkthrough with real output: **[Quickstart](https://docs.isalgo.dev/arda/quickstart.html)**.

## Why arda

**IgBLAST-quality regions, without IgBLAST at runtime.** arda does the expensive work once,
offline, building a reference of every in-frame V·J germline scaffold already annotated by IgBLAST.
At runtime it maps a query to that reference with MMseqs2 and projects the markup through the
alignment in C++. Region concordance with IgBLAST is **98–99.7 %** on ~7,300 real GenBank mRNA
records across five organisms. IgBLAST is needed only to rebuild the reference.

**Precise on bulk RNA-seq.** Over 16 datasets where every compared tool ran (5,273 real fragments,
Wilson 95 % CIs): recall **0.986** [.982–.989], precision lower bound **0.889** [.881–.897]. Recall
ties TRUST4; precision does not — arda's lower bound sits above every competitor's upper bound.

**Cheap, not only fast.** Mapping holds a flat **300–650 MB** at any read depth, which is what lets
many samples share a node. On a 100,000-read TRA amplicon arda uses **3.2× less peak memory** and
**1.4× less CPU** than MiXCR, and returns the most clonotypes over the most reads.

**It declines rather than guessing.** A D call carries the E-value it was accepted on. A junction
whose Cys104 anchor is not in the read is not emitted. A germline with no derivable anchor is
flagged, not invented. A `j_call` requires J evidence rather than inheriting it from the scaffold.

**Embeddable and scriptable.** `import arda; arda.annotate_sequences(...)` returns AIRR record
dicts — no subprocess, no temp files. The CLI is plain files in, plain files out, so it drops into
any workflow engine; the same sample sheet drives the CLI, a SLURM array, the Nextflow module and
the Snakemake workflow, and a sharded run is byte-identical to a single-node one.

Every figure above is measured. Methods, repetition spread, and what each comparison is and is not
entitled to claim: **[Benchmarks](https://docs.isalgo.dev/arda/benchmarks.html)**.

## Documentation

| | |
|---|---|
| **[Quickstart](https://docs.isalgo.dev/arda/quickstart.html)** | Install, run a library, read the table that comes out |
| **[Choosing a mode](https://docs.isalgo.dev/arda/usage.html)** | Bulk, amplicon or single cell, and what each writes |
| **[Output reference](https://docs.isalgo.dev/arda/outputs.html)** | Every file a run produces, and every column in it |
| **[Command reference](https://docs.isalgo.dev/arda/cli.html)** | All 22 commands, grouped by task |
| **[Recipes](https://docs.isalgo.dev/arda/examples.html)** | Copy-paste shell and polars snippets for post-run analysis |
| **[Use cases](https://docs.isalgo.dev/arda/use_cases.html)** | Monoclonal QC, low-frequency variants, SHM, tool comparisons |
| **[Single cell](https://docs.isalgo.dev/arda/singlecell.html)** | Reference-free per-cell contigs, chain pairing, doublets |
| **[Quality control](https://docs.isalgo.dev/arda/qc.html)** | Per-sample metrics, cohort roll-up, the HTML dashboard |
| **[At scale](https://docs.isalgo.dev/arda/cluster.html)** | SLURM sharding, [sample sheets](https://docs.isalgo.dev/arda/samples.html), [Nextflow and Snakemake](https://docs.isalgo.dev/arda/pipeline_integration.html) |
| **[How it works](https://docs.isalgo.dev/arda/how_it_works.html)** | The offline reference, the runtime path, the C++ extensions |
| **[Background](https://docs.isalgo.dev/arda/error_correction.html)** | [Error correction](https://docs.isalgo.dev/arda/error_correction.html), [SHM](https://docs.isalgo.dev/arda/shm.html), [D segments](https://docs.isalgo.dev/arda/d_segments.html), [productivity](https://docs.isalgo.dev/arda/productivity.html), [genotypes](https://docs.isalgo.dev/arda/genotype.html), [recombination models](https://docs.isalgo.dev/arda/scenarios.html) |
| **[Glossary](https://docs.isalgo.dev/arda/glossary.html)** | Junction vs CDR3, scaffold, read group, regime, and the rest |
| **[API](https://docs.isalgo.dev/arda/api.html)** | The Python library |

## From Python

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],   # strings or (id, seq) pairs
    seqtype="nt", organism="human",
)
```

Each record is a dict of AIRR fields — `v_call`, `d_call`/`d2_call`, `j_call`, `c_call`/`c_class`,
`fwr1`–`fwr4`, `cdr1`–`cdr3` with 1-based closed `*_start`/`*_end` and `*_aa`, `junction(_aa)`,
`np1`–`np3`, `d_support`, the per-segment CIGARs, `v_mutations`/`j_mutations`, `productive`, and
more. The TSV form passes `airr.schema` validation.

Reading arda's output — note that quoting is off, both when writing and when reading:

```python
import polars as pl

clones = pl.read_csv("out/SAMPLE.clones.tsv", separator="\t", quote_char=None)
top = (clones.filter(pl.col("locus") == "TRB")
       .with_columns((pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias("frequency"))
       .sort("duplicate_count", descending=True))
```

## Examples

[`examples/`](examples/) is a runnable tour, every artifact derived from real data committed to this
repository and regenerated by `python examples/regenerate.py`: one real mRNA per locus; the two
human reads (of 7,341, across five organisms) carrying a tandem D-D; seven VDJdb records covering
every junction-repair outcome; and a 1,035-read FASTQ that runs the whole bulk RNA-seq pipeline in
about six seconds.

## Development

```bash
bash setup.sh                                         # uv venv, C++ extensions, IgBLAST, mmseqs
source .venv/bin/activate

python -m pytest tests/unit tests/synthetic tests/realworld -q   # the CI gate
ruff check src/
make -C docs html                                                # zero warnings required
```

`tests/realworld` compares against IgBLAST on committed fixtures and runs offline.
`tests/benchmark` is opt-in with `RUN_BENCHMARK=1`. Install `-e '.[test]'` before reading a green
suite as full coverage — optional extras gate optional suites.

Layout: the package is `src/arda/` (`cli.py`, `annotate/`, `rnaseq/`, `refbuild/`), the four
nanobind C++ extensions are `src/_markup/`, `src/_prefilter/`, `src/_segmap/`, `src/_denoise/`, the
committed references are in `database/`, and `integrations/` holds the Nextflow module and Snakemake
workflow. Contributor guidance is in [`CLAUDE.md`](CLAUDE.md); design records are in
[`project/`](project/); planned work is in [`ROADMAP.md`](ROADMAP.md); data provenance is in
[`SOURCES.md`](SOURCES.md); per-release changes are in [`CHANGELOG.md`](CHANGELOG.md).

## Related tools

arda is part of the [antigenomics](https://github.com/antigenomics) stack and stays deliberately
inside annotation. Repertoire biology — diversity, clonality, rarefaction, overlap, cross-sample
clonotype matching — belongs to [vdjtools](https://github.com/antigenomics/vdjtools), which takes
`arda-mapper` as a dependency.

## License

GPLv3. See [`LICENSE`](LICENSE).
