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

## Why

IgBLAST is the gold standard but is slow to invoke per-batch and awkward to embed.
`arda` keeps IgBLAST-quality region calls while being:

- **Fast & scalable** — MMseqs2 search + a C++ projection step; multiprocessing
  and SLURM-friendly from small FASTA to large FASTQ.
- **Embeddable** — `import arda; arda.annotate_sequences(...)`.
- **Honest** — a D call is gated on an E-value that ships with it (`d_support`), a
  germline allele with no derivable anchor is flagged rather than guessed, and a
  repair that would not produce a canonical junction is refused.
- **Easy to install** — `pip install arda-mapper` (binary wheels ship the C++
  extension); the `mmseqs` binary is fetched as a static build at runtime — no
  conda. IgBLAST is fetched into a gitignored `bin/` and is only needed to
  (re)build the reference DB, not at runtime.

## Install

```bash
pip install arda-mapper   # from PyPI (imports as `arda`); binary wheels ship the C++ extension
```

`mmseqs2` (the search backend) is fetched/managed by arda at runtime. For development — and to
get the committed germline references on disk — use `setup.sh`:

```bash
bash setup.sh            # uv .venv, fetches IgBLAST + static mmseqs, editable install
source .venv/bin/activate
```

Needs [uv](https://docs.astral.sh/uv/). Flags: `--build-db` (rebuild references after
install), `--tests` (run the fast suites). The committed `database/vdj/<organism>/`
references mean **most users never need to build anything**. A `pip install arda-mapper`
with no source checkout **auto-fetches** the curated references into `~/.cache/arda` on
first use (the `arda-reference-vdj.tar.gz` release asset) and builds the MMseqs2 index
there — **no `$ARDA_HOME` and no build step required** (set `ARDA_NO_AUTO_FETCH` for
air-gapped runs with a pre-populated cache).

Supported organisms: **human, mouse** (full IG + TR), **rat, rabbit, rhesus_monkey**
(IG only — IgBLAST ships no TR internal annotation for these).

## CLI

```bash
arda info                                   # resolved paths + tool availability
arda annotate -i reads.fastq -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta  -o out.airr.tsv --organism human --seqtype aa
arda annotate -i reads.fastq -o out.airr.tsv --strand forward   # plus-strand only
arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report -       # mark up + repair bare (CDR3aa, V, J) records
arda rnaseq map --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv  # filter receptor reads from bulk RNA-seq
arda rnaseq assemble -i mapped.airr.tsv -o assembled.airr.tsv   # rescue CDR3s no single read spans
arda rnaseq correct -i mapped.airr.tsv -o clones.tsv            # collapse CDR3 errors into clonotypes
arda rnaseq run --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # one-shot map+assemble+correct for pipelines
arda igblast -i reads.fastq -o truth.airr.tsv                   # gold-standard IgBLAST (all loci)
arda build-db   --organism all              # rebuild references (needs IgBLAST)
arda build-index --organism all             # (re)build the precompiled mmseqs DBs
arda slurm -i big.fastq -o big.airr.tsv --shards 50 --partition cpu   # cluster scale
```

[`examples/`](examples/) is a runnable tour, every artifact derived from real data committed to
this repo and regenerated by `python examples/regenerate.py`: one real mRNA per locus; the two
human reads (of 7,341, across five organisms) that carry a **tandem D-D**; seven VDJdb records
covering every junction-repair outcome, including one arda reports and refuses to rewrite and one
it refuses outright; and a 1,035-read FASTQ that runs the whole bulk RNA-seq pipeline in ~6 s.
Two tests re-run that script and fail if a committed artifact stops reproducing.

See [`CHANGELOG.md`](CHANGELOG.md) for what changed per release and
[`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) for measured speed and accuracy.

The reference database ships with **precompiled MMseqs2 indexes**
(`database/vdj/<organism>/mmseqs/`), so annotation runs out of the box with no
build step. They are used automatically when the local MMseqs2 version matches the
shipped one; otherwise arda transparently rebuilds a private cache on first run
(`arda build-index` regenerates the shipped DBs for your version).

Input may be FASTA or FASTQ, plain or gzipped. Nucleotide input is searched on
**both strands** by default (reverse-complement reads are re-oriented and flagged
`rev_comp=T`); a single search annotates a mixed bulk RNA-seq file across all loci.

## Pipeline integration

`arda rnaseq run` is a one-shot `map`+`assemble`+`correct` for bulk RNA-seq: given paired (or single)
gzipped FASTQ it writes `<prefix>.clones.tsv` (AIRR clonotypes), `<prefix>.airr.tsv` (mapped reads),
`<prefix>.assembled.airr.tsv` (assembled long-CDR3 reads) and `<prefix>.arda.json` (run report).
Because it is a plain CLI over named files, it drops into any workflow engine with no glue code.

A ready-to-use **Nextflow module** lives in [`integrations/nextflow/arda/`](integrations/nextflow/arda/):
copy it to `modules/local/arda/` in an nf-core/rnaseq (or similar) checkout, feed it the trimmed
per-sample FASTQ channel the aligners already use, and it publishes per-sample clonotype tables to
`${params.outdir}/arda/`. It ships a conda `environment.yml` (works with `-profile conda` out of the
box) and a `Dockerfile`, and emits a `versions.yml`. See its
[README](integrations/nextflow/arda/README.md) and the
[pipeline-integration guide](https://docs.isalgo.dev/arda/pipeline_integration.html) for the five-line drop-in.

arda is **CPU-bound**: ~40–50k reads/s on 32 cores, so a full-depth bulk RNA-seq sample (~50 M read
pairs) finishes in ~45 min. Throughput scales with cores. Peak memory tracks **repertoire richness**,
not read depth — mapping is flat (~300–400 MB at any depth), while Stage 3 holds the clone set, so a
B-cell-rich tumour peaked at 2.7 GB (28k clonotypes) and a colder sample with *more* reads used 314 MB.
Budget ~4 GB.

## Library

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],  # strings or (id, seq) pairs
    seqtype="nt", organism="human",
)
# -> list of AIRR record dicts: v_call, d_call/d2_call, j_call, c_call/c_class,
#    fwr1..fwr4, cdr1..cdr3, *_start/*_end (1-based closed), *_aa, junction(_aa),
#    np1/np2/np3, d_support/d2_support, {v,j,c,d}_cigar, sequence_alignment,
#    germline_alignment, productive, ...
# The TSV is a spec-valid AIRR Rearrangement file (passes airr.schema validation).
```

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

map_d_junction(junction_nt, "TRDV1*01", "TRDJ1*01", "human").d2_call   # 'TRDD3*01'
posterior_d("CASSPLGQAYEQYF", "TRBV5-1*01", "TRBJ2-7*01", "human").d_call  # 'TRBD1'
```

**Coordinates here are *junction* space** — Cys104 through Phe/Trp118, **both anchors
included**. That is what VDJdb's `cdr3` column holds, and it is *not* arda's `cdr3` field,
which excludes both. Conflating them silently corrupts every coordinate.

Repair is conservative: only edits adjacent to a conserved anchor are applied, everything
deeper is *reported* and left alone, and a repair is refused outright unless the result opens
with Cys104 and closes with Phe/Trp118.

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

## Bulk RNA-seq mode

`arda rnaseq` is a recall-first pipeline for extracting the receptor repertoire from
bulk RNA-seq, where 1–5% of reads are receptor-derived:

- **`map`** streams paired FASTQ, keeps only reads that map to a receptor scaffold, and
  writes them as AIRR. The reference includes **`J + C` constant-region scaffolds**, so a
  read spanning the J→C splice — which ends in the constant region and has no V to anchor —
  still maps, and carries a **`c_call`** (the CH1 exon) plus a **`c_class`** isotype
  (`IGHG`/`IGHM`/`IGHA` … — the class, never the noise-prone subclass). In paired mode the
  isotype of a CDR3-bearing read is recovered from its constant-region mate. `--reconstruct`
  merges each overlapping mate pair into one fragment, resolving overlap mismatches by the
  higher-Phred base.
- **`assemble`** (Stage 3) reconstructs clonotypes whose CDR3 is too long for any single
  100–150 bp read to span (V(DD)J ultralong, ~20–40 aa) by greedy overlap-extension anchored on
  Stage-1's per-read `cdr3_start`, and folds the recovered reads back into `correct`.
- **`correct`** aggregates reads into clonotypes and collapses sequencing-error CDR3 variants.
  Abundance is the AIRR **`duplicate_count`** — every read that *encompasses* the junction
  (spanning or partial), the true expression estimate — with **`consensus_count`** for distinct
  fragments. The error model is per-base with a **length-scaled threshold** (a mismatch over a
  longer junction is likelier an error) and is SHM-indel-tolerant, keeping only complete junctions.

```bash
arda rnaseq run --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # one-shot map + assemble + correct
```

`arda igblast -i reads.fastq -o truth.airr.tsv` runs IgBLAST across all loci as a
gold-standard reference for benchmarking (see the `arda-benchmark` project).

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
   the D lives in. The D call is then accepted on a Karlin–Altschul E-value (`d_support`)
   rather than a per-locus score floor, and is constrained by germline geometry: TRBD2
   lies 3′ of the entire TRBJ1 cluster, so a TRBJ1 rearrangement can never be assigned
   TRBD2. D mapping also runs on `--seqtype aa`, against each D germline's three
   translated frames.
3. **Bare records** (`arda.cdr3fix`, `arda.dpost`): a VDJdb-style row — CDR3 amino acid,
   V, J, species, and no read — is marked up against the same anchors (`arda markup`),
   its errors located and conservatively repaired, and optionally given a D gene inferred
   from the junction *length* (`--d-posterior`).

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
out-of-frame junctions, **D-segment mapping incl. D-D fusions** (IGH/TRB/TRD, on
nucleotide *and* protein input, E-value gated, genomic-order constrained),
**constant-region `J + C` scaffolds** (`c_call`/`c_class` isotype), **bulk RNA-seq
mode** (`rnaseq map`/`assemble`/`correct`/`run`), **long-CDR3 contig assembly**,
**coverage-based expression** (`duplicate_count`/`consensus_count`), **junction markup
+ repair on bare records** (`arda markup`, `arda.cdr3fix`), **D on a bare junction**
(`arda.annotate.dmap`) and **the aa D posterior** (`arda.dpost`), precompiled indexes,
**multi-node (SLURM) sharding**. Next: full-depth clonotype benchmarking, and an
`arda.hmm` semi-Markov model of V→N→D→N→J that would subsume the E-value gate, the
genomic-order constraint and the D posterior into one forward–backward pass.

## Development

```bash
pip install -e .                                      # rebuilds the C++ ext on import
python -m pytest tests/unit tests/synthetic -q        # fast suite (needs mmseqs)
python -m pytest tests/realworld -q                   # vs IgBLAST, on committed fixtures — offline
env RUN_BENCHMARK=1 python -m pytest tests/benchmark -s   # timing/memory/scaling
```

Optional extras gate optional suites:
`.[groundtruth]` (`olga`) for the generative ground-truth tests that keep `arda.cdr3fix`
honest, `.[test]` for `airr` schema validation. Without them those tests **skip**, so
`pip install -e '.[test]'` before reading a green suite as full coverage.

Layout: `src/arda/{refbuild,annotate}`, C++ in `src/_markup/markup.cpp`,
references in `database/`, downloads in gitignored `bin/` + `data/`. Design rationale and
the gotchas that cost us the most live in [`memory/`](memory/) — read
[`memory/d-mapping.md`](memory/d-mapping.md) and
[`memory/junction-markup.md`](memory/junction-markup.md) before touching either.