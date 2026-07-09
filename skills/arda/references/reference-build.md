# Reference Build (offline)

How the curated `database/vdj/<organism>/` references are produced. **Most users
never run this** — references ship with the package. Build only to add or refresh
an organism. Requires IgBLAST (fetched into `bin/` by `setup.sh`).

**Contents:** pipeline · CLI · outputs · manifest · when to rebuild.

## Pipeline (`arda.refbuild`)

1. **Download** IMGT/V-QUEST germlines for the organism.
2. **Enumerate** deduplicated, in-frame **V×J scaffolds**. D segments are *not*
   enumerated (they only affect the CDR3 interior and are mapped at runtime).
   Scaffolds are `V + N*pad + J`, with padding chosen to preserve the J frame.
3. **Annotate** each scaffold once with `igblastn -outfmt 19` (the expensive,
   offline step).
4. **Translate** and **write** `database/vdj/<organism>/`:
   - `alleles.fasta`, `alleles.aa.fasta` — scaffold sequences (nt / aa)
   - `markup.tsv`, `markup.aa.tsv` — per-scaffold region coordinates + sequences
   - `combinations.tsv` — scaffold → (V, J) allele pairs, padding
   - `d_germlines.fasta` — D germlines for VDJ loci
   - `loci_manifest.tsv` — per-locus reference coverage, one row per defined locus
   - `build.log`
   - `mmseqs/<nt|aa>/` — precompiled MMseqs2 indexes (+ `VERSION`)

At runtime, arda projects a scaffold's markup onto each query via the C++
`transfer_regions` hot path — no IgBLAST involved.

## Locus coverage manifest

`loci_manifest.tsv` records reference coverage for **every** defined locus (TRA,
TRG, IGK, IGL, IGH, TRB, TRD), one row each — so a locus with no reference is
*visible* instead of silently absent. Columns: `locus`, `group`,
`n_vj_scaffolds`, `n_jc_scaffolds`, `n_d_germlines`, `unreachable_d_germlines`,
`status` (`ok` / `EMPTY`). A locus is `EMPTY` when the build produced no
scaffolds — reads from it fall through unannotated.

`build_species` warns at build end on:

- **Empty loci** — any locus with no scaffolds (`N/7 loci have NO reference
  scaffolds`). Only human and mouse have a TCR reference; IMGT has no TR loci for
  rat / rabbit / rhesus, so their `TRA`/`TRG`/`TRB`/`TRD` come out `EMPTY`.
- **Unreachable D germlines** — a VDJ locus that ships D germlines but no
  scaffolds. Runtime D lookup is keyed by a hit scaffold's locus, so with no
  scaffold those germlines can never be reached (dead weight).

## CLI

```bash
arda build-db --organism all        # full offline rebuild (needs IgBLAST in bin/)
arda build-db --organism human
arda build-index --organism all     # only rebuild the mmseqs indexes for your mmseqs version
```

`build-index` is the lightweight one: use it when the shipped indexes don't match
your local MMseqs2 version (arda otherwise rebuilds a private cache in `data/`).

## IgBLAST

IgBLAST binaries are downloaded into `bin/` (gitignored) by `setup.sh` /
`scripts/fetch_igblast.py` — the conda-forge build lags the NCBI release, so arda
fetches the NCBI release directly. IgBLAST is a **build-time-only** dependency.

## When to rebuild

- Adding a new organism or refreshing IMGT germlines → `build-db`.
- Local mmseqs version differs from the shipped indexes → `build-index` (or just
  let arda build a private cache on first run).
- A normal install / annotation run needs neither — the committed references are
  authoritative.
