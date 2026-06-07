# arda — Antigen Receptor Domain Annotation

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
bash setup.sh            # creates conda env `arda`, fetches IgBLAST, pip install -e .
conda activate arda
```

Flags: `--no-conda` (use the active env), `--build-db` (rebuild references after
install), `--tests` (run the fast suites). The committed `database/vdj/<organism>/`
references mean **most users never need to build anything**.

Supported organisms: **human, mouse** (full IG + TR), **rat, rabbit, rhesus_monkey**
(IG only — IgBLAST ships no TR internal annotation for these).

## CLI

```bash
arda info                                   # resolved paths + tool availability
arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa
arda annotate -i reads.fastq -o out.airr.tsv --strand forward   # plus-strand only
arda build-db --organism all                # rebuild references (needs IgBLAST)
```

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
# -> list of AIRR record dicts: v_call, j_call, fwr1..fwr4, cdr1..cdr3,
#    *_start/*_end (1-based closed), *_aa, junction, productive, ...
```

## How it works

1. **Reference build** (`arda.refbuild`, offline): download IMGT/V-QUEST germlines
   → enumerate deduplicated in-frame **V×J** scaffolds (D only affects CDR3
   interior, so it isn't enumerated) → annotate with `igblastn -outfmt 19` →
   translate → write `database/vdj/<organism>/{alleles.fasta, alleles.aa.fasta,
   markup.tsv, markup.aa.tsv, combinations.tsv, build.log}`.
2. **Runtime** (`arda.annotate`): MMseqs2 search query→scaffolds → best hit →
   C++ `transfer_regions` projects scaffold region coordinates onto the query
   (handling indels, truncation, mid-codon alignment starts) → AIRR TSV.

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

- **D-segment mapping** — scaffolds are V·J only, so `v_call`/`j_call` and all
  FR/CDR coordinates are assigned but `d_call` is not yet. Planned: align the CDR3
  interior against a D germline DB after V/J transfer, including **double D-D
  junctions** (D-D fusions → `d_call` + `d2_call`).
- Multi-node sharding (single-node streaming + threading is implemented).

## Development

```bash
pip install -e .                                  # rebuilds the C++ ext on import
python -m pytest tests/unit tests/synthetic -q    # fast suite
env ARDA_REALWORLD=1 python -m pytest tests/realworld -s   # vs IgBLAST (network)
env RUN_BENCHMARK=1   python -m pytest tests/benchmark -s  # timing/memory/scaling
```

Layout: `src/arda/{refbuild,annotate}`, C++ in `src/_markup/markup.cpp`,
references in `database/`, downloads in gitignored `bin/` + `data/`.