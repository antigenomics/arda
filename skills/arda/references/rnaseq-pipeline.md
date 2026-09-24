# The three-stage pipeline, run reports and run QC

Loaded on demand from `SKILL.md`. `map` -> `assemble` -> `correct`, what each stage owns, what the
run report and the QC table carry, and which numbers to size a job from.

## Bulk RNA-seq mode (`arda rnaseq`)

For libraries where only 1–5% of reads are receptor-derived. Three stages, run separately or in
one shot with `arda rnaseq` / `arda amplicon` (which do all three by default). `seqtree` and
`dnaio` are core deps; the `rnaseq` extra is an empty alias kept so old pins resolve.

**`map`** — streams paired FASTQ (`--r1`/`--r2`), keeps only reads mapping to a receptor
scaffold, writes them as AIRR. Recall-first, with `--min-score`/`--kmer`/`--max-seqs` around
one default preset.

- The reference includes `J + C` constant-region scaffolds, so a read spanning the J→C splice
  (no V, hence no junction) still maps and carries `c_call`/`c_class`. In paired mode a
  CDR3-bearing read gets its isotype from its constant-region mate.
- `--reconstruct` merges overlapping mates into one fragment, giving a short read the mate's
  V/J context; overlap mismatches resolve to the higher-Phred base. FASTQ quality is read only
  on this path, so the default stays fast.

**`assemble`** (Stage 3) — recovers clonotypes whose CDR3 no single 100–150 bp read spans
(V(DD)J ultralong, ~20–40 aa), by anchored greedy overlap-extension over Stage-1's per-read
`cdr3_start`. It carries the contig's D call onto every member read: an ultralong CDR3 is
where a tandem D-D is both most likely and least visible to one read.

> `annotate.contig` gives an assembled contig its AIRR cigars two ways, producing the same
> record: `reannotate_contigs` (re-align it — what `assemble` uses) and `merge_contig` (stitch
> the reads' existing alignments via C++ `_markup.merge_alignment`). Merge is ~9× faster at
> ~10⁵ contigs/sample (scRNA-seq) and is the intended default once the assembler emits read
> layouts.

**`correct`** — collapses sequencing-error CDR3 variants into clonotypes keyed by
`(locus, v_call, j_call, junction)`.

- Abundance is the AIRR **`duplicate_count`** (every read encompassing the junction), with
  **`consensus_count`** for distinct fragment consensuses. There is no `count` column.
- A neighbour is an error *child* when `count[parent] * p_sub**n_subs * p_ind**n_indel >=
  count[child]`. Knobs: `--max-subs`, `--max-indel`, `--error-rate`, `--indel-rate` (per-BASE,
  length-scaled), `--require-vj`, `--error-method` (`simple|binom|betabinom`), `--complete-only`
  (on by default).
- Row order is deterministic — abundance ties break on `(junction, v_call, j_call)`.
- Each clonotype's D is mapped once into its *corrected* junction (`d_call`/`d2_call`/
  `d_support`), not voted over reads: D is a function of the junction, and a read's copy of it
  carries sequencing error.

## Run reports: `peak_rss_mb` is monotone, by design

Every stage reports `wall_seconds`, `peak_rss_mb` and `rss_gain_mb`. `peak_rss_mb` is the
**whole-process** high-water mark *as of that stage's end* (getrusage offers no per-stage
reset), which is exactly what a SLURM `--mem` or Nextflow `memory` directive must cover;
`rss_gain_mb` is that stage's contribution. **Budget for Stage 3, not Stage 1**: mapping is flat
at ~300–650 MB at any depth, but the clone set scales with repertoire richness — Stage-3
`correct` peaked at **2,071.7 MB** on a B-cell-rich tumour (28,444 clonotypes from 105 M reads),
versus **549 MB** for a colder sample with *more* reads (139 M). Budget ~4 GB.

## Run QC: `arda stats`

Every mode run also writes `<prefix>.stats.tsv`; `arda stats` builds the same table from any
subset of a run's artifacts (`--airr`, `--clones`, `--report`, `--r1`/`--r2`), reading only what
already exists — no re-read of the FASTQ, no alignment.

Four columns, `scope` / `key` / `metric` / `value`, one value per cell. Scopes: `run` (the report
verbatim — reads, FASTQ bytes, read length, paired, threads, wall time, peak RSS), `sample`,
`chain` (per locus, **reads AND clonotypes**: productive/non-functional, stop codons, out-of-frame,
truncated junctions, junction length min/max/mean, junction quality, SHM rate, chimeras),
`v_gene`/`j_gene`, and `allele_candidate`.

Never: Long, not wide — the metric set differs per scope, so a wide table is mostly empty cells.
Never: A metric with **no input is omitted, never emitted as 0**: a run without `--junction-quality`
has no `junction_quality_mean` row rather than a zero that reads like a terrible library.
Never: Truncation, a stop codon and an out-of-frame junction are counted **separately**; `_COMPLETE`
folds them together and a QC table must not.
⚠ `allele_candidate` is a **shortlist, not a genotype call** — arda does not genotype. A novel
allele, SHM and a miscall are the same string in the mutation list; recurrence within the allele
(`--allele-min-frac`) and Phred are what separate them, and both are reported per variant.
Never: Chimera / non-functional / stop-codon counts are **flags, never filters**.

The two quality columns feeding it are opt-in on `map` and use **different encodings**:
`--junction-quality` writes raw Phred+33 *characters* over `junction`; `--mutation-quality` writes
comma-joined *integers*, one per `v_mutations` / `j_mutations` entry, in the same order.
