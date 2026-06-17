<p align="center">
  <picture>
    <source media="(prefers-color-scheme: light)" srcset="assets/arda_light.svg">
    <img alt="arda" src="assets/arda_dark.svg" width="340">
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
is an [AIRR](https://docs.airr-community.org/)-formatted annotation that matches
IgBLAST (≈97% region concordance on real GenBank mRNA), from a plain CLI + Python
library — no Docker, no workflow engine.

## Why

IgBLAST is the gold standard but is slow to invoke per-batch and awkward to embed.
`arda` keeps IgBLAST-quality region calls while being:

- **Fast & scalable** — MMseqs2 search + a C++ projection step; multiprocessing
  and SLURM-friendly from small FASTA to large FASTQ.
- **Embeddable** — `import arda; arda.annotate_sequences(...)`.
- **Easy to install** — conda for the `mmseqs` binary, `pip install -e .` for the
  package + C++ extension; IgBLAST is fetched into a gitignored `bin/` and is only
  needed to (re)build the reference DB, not at runtime.

## Install

```bash
pip install arda-mapper   # from PyPI (imports as `arda`); binary wheels ship the C++ extension
```

`mmseqs2` (the search backend) is fetched/managed by arda at runtime. For development — and to
get the committed germline references on disk — use `setup.sh`:

```bash
bash setup.sh            # creates conda env `arda`, fetches IgBLAST, pip install -e .
conda activate arda
```

Flags: `--no-conda` (use the active env), `--build-db` (rebuild references after
install), `--tests` (run the fast suites). The committed `database/vdj/<organism>/`
references mean **most users never need to build anything**; a `pip`-only install builds them
on first use (or set `$ARDA_HOME` to a checkout).

Supported organisms: **human, mouse** (full IG + TR), **rat, rabbit, rhesus_monkey**
(IG only — IgBLAST ships no TR internal annotation for these).

## CLI

```bash
arda info                                   # resolved paths + tool availability
arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa
arda annotate -i reads.fastq -o out.airr.tsv --strand forward   # plus-strand only
arda build-db   --organism all              # rebuild references (needs IgBLAST)
arda build-index --organism all             # (re)build the precompiled mmseqs DBs
arda slurm -i big.fastq -o big.airr.tsv --shards 50 --partition cpu   # cluster scale
```

See [`examples/`](examples/) for a runnable per-locus demo and
[`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) for measured speed/accuracy.

The reference database ships with **precompiled MMseqs2 indexes**
(`database/vdj/<organism>/mmseqs/`), so annotation runs out of the box with no
build step. They are used automatically when the local MMseqs2 version matches the
shipped one; otherwise arda transparently rebuilds a private cache on first run
(`arda build-index` regenerates the shipped DBs for your version).

Input may be FASTA or FASTQ, plain or gzipped. Nucleotide input is searched on
**both strands** by default (reverse-complement reads are re-oriented and flagged
`rev_comp=T`); a single search annotates a mixed bulk RNA-seq file across all loci.

## Library

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],  # strings or (id, seq) pairs
    seqtype="nt", organism="human",
)
# -> list of AIRR record dicts: v_call, d_call/d2_call, j_call, fwr1..fwr4,
#    cdr1..cdr3, *_start/*_end (1-based closed), *_aa, junction(_aa), np1/np2/np3,
#    v_sequence_end, j_sequence_start, productive, rev_comp, ...
```

### Annotating bare germline segments

There is no coverage filter, so a **V-only** or **J-only** query maps to its
scaffold and only the regions inside the query's coverage are returned. This lets
you annotate isolated germline V or J alleles without synthesising a
rearrangement — a bare V yields `fwr1..fwr3`, a bare J yields `fwr4`:

```python
from arda.annotate.mapper import annotate_records

recs = annotate_records(
    [("TRBV9*01", v_germline_nt), ("TRBJ2-7*01", j_germline_nt)],
    organism="human", seqtype="nt", strand="forward", map_d=False,
)
# V record -> fwr1/cdr1/fwr2/cdr2/fwr3 (+ v_sequence_end = CDR3 start)
# J record -> fwr4 (+ j_sequence_start = CDR3 end / FR4 start)
```

(mirpy uses exactly this to bake per-allele FR/CDR subsequences into its gene
library; see `tests/synthetic/test_germline_segments.py`.)

## How it works

1. **Reference build** (`arda.refbuild`, offline): download IMGT/V-QUEST germlines
   → enumerate deduplicated in-frame **V×J** scaffolds (D only affects CDR3
   interior, so it isn't enumerated) → annotate with `igblastn -outfmt 19` →
   translate → write `database/vdj/<organism>/{alleles.fasta, alleles.aa.fasta,
   markup.tsv, markup.aa.tsv, combinations.tsv, build.log}`.
2. **Runtime** (`arda.annotate`): MMseqs2 search query→scaffolds → best hit →
   C++ `transfer_regions` projects scaffold region coordinates onto the query
   (handling indels, truncation, mid-codon alignment starts, reverse strand) → for
   VDJ loci a gapless C++ local alignment of the CDR3 interior against the D
   germlines adds `d_call`/`d2_call` + `np*` → AIRR TSV. Out-of-frame junctions are
   reported with an N-bridge (`_`) so FR4 still reads.

See [`memory/`](memory/) for design rationale and gotchas. Fast sequence
primitives (`translate`, `detect_coding_frame`, `reverse_complement`,
`back_translate`) live in the C++ extension and are re-exported from
`arda.refbuild.translate` — mirpy-API-compatible, so mirpy can `import arda` and
reuse them.

## Performance

Exact annotation that matches IgBLAST while being several times faster, scaling to
large FASTQ. Synthetic human IGH, 16 threads (`scripts/bench_vs_igblast.py`):

| sequences | arda | arda rate | speedup vs IgBLAST | region concordance |
|----------:|-----:|----------:|-------------------:|-------------------:|
| 10,000    | 5.5s | ~1.8k/s   | 4.4×               | 98.9%              |
| 50,000    | 16s  | ~3.0k/s   | 7.3×               |                    |
| 100,000   | 30s  | ~3.3k/s   | 7.9×               |                    |

On ~7.3k real GenBank mRNA records spanning **all five organisms and their loci**
(committed, gzipped test fixtures), region concordance with IgBLAST on productive
records is **98–99.7%** per organism; `junction_aa`/`cdr3_aa` match IgBLAST ~99%
and satisfy the AIRR invariants exactly. V-gene assignment agrees ~100%. (GenBank
also contains genomic/partial/non-productive entries that confuse both tools; those
are excluded from the comparison.)

**Bulk RNA-seq is much faster than amplicon**, because mmseqs prefilters by k-mer
matching — reads with no receptor k-mer are rejected before alignment. At 150 nt
reads, 16 threads (`scripts/bench_prefilter.py`):

| receptor content | throughput |
|-----------------:|-----------:|
| 100% (amplicon)  | ~5.7k reads/s |
| 10%              | ~19k reads/s |
| 1% (blood RNA-seq) | ~25k reads/s |

Extrapolated to a **32-core node**, a 30M-read bulk RNA-seq library (~1% receptor)
annotates in roughly **10–20 min** — the same order of magnitude as a STAR genome
alignment pass on the same data (STAR is faster per read, but arda maps only to a
tiny germline DB and the non-receptor majority costs just prefilter rejection).
Large FASTQ is **streamed in bounded chunks** (a background reader prefetches the
next chunk while the current one is annotated), so memory stays flat regardless of
input size — `--chunk-size` tunes it.

## Roadmap / TODO

See [`ROADMAP.md`](ROADMAP.md). Done: V·J reference build (5 organisms), MMseqs2
mapping, C++ markup transfer, reverse-complement, all-loci querying, streaming I/O,
out-of-frame junctions, **D-segment mapping incl. D-D fusions**, precompiled
indexes, **multi-node (SLURM) sharding**. Next: full AIRR productivity.

## Development

```bash
pip install -e .                                  # rebuilds the C++ ext on import
python -m pytest tests/unit tests/synthetic -q    # fast suite
env ARDA_REALWORLD=1 python -m pytest tests/realworld -s   # vs IgBLAST (network)
env RUN_BENCHMARK=1   python -m pytest tests/benchmark -s  # timing/memory/scaling
```

Layout: `src/arda/{refbuild,annotate}`, C++ in `src/_markup/markup.cpp`,
references in `database/`, downloads in gitignored `bin/` + `data/`.