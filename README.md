# arda — Antigen Receptor Domain Annotation

Fast FR/CDR region annotation for **TCR** and **BCR** nucleotide *and* amino-acid
sequences (MHC groove/helix annotation is scaffolded for a future release).

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
arda build-db --organism all                # rebuild references (needs IgBLAST)
```

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

See [`memory/`](memory/) for design rationale and gotchas.

## Development

```bash
pip install -e .                                  # rebuilds the C++ ext on import
python -m pytest tests/unit tests/synthetic -q    # fast suite
env ARDA_REALWORLD=1 python -m pytest tests/realworld -s   # vs IgBLAST (network)
env RUN_BENCHMARK=1   python -m pytest tests/benchmark -s  # timing/memory/scaling
```

Layout: `src/arda/{refbuild,annotate}`, C++ in `src/_markup/markup.cpp`,
references in `database/`, downloads in gitignored `bin/` + `data/`.