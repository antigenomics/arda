# The three-stage pipeline, run reports and quality control

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
`v_gene`/`j_gene`, `allele_candidate`, and five **distribution** scopes keyed `locus:bucket` —
`junction_aa_len`, `read_len` (10-nt buckets of the aligned length), `clone_size` (powers of two,
weighted by clonotypes and by reads), `isotype`, and `chain_support` (`arda cells`). They exist
because min/max/mean cannot show a shape: a bimodal junction length is two primer sets in one
tube, a read-length cliff is an adapter left on, and a clone-size distribution with no singletons
is an over-amplified library. Each reads as an ordinary mean. Only occupied buckets get a row.

`arda cells` writes the same table in the same scopes, so a cohort can mix single-cell and bulk;
it adds `pairing_rate`, `doublet_rate` and `molecules_placed_fraction`. Every run also writes
`<prefix>.stats.json` — the same rows, typed and nested — built from the same list, so the two
cannot diverge.

Never: Long, not wide — the metric set differs per scope, so a wide table is mostly empty cells.
Never: A metric with **no input is omitted, never emitted as 0**: a run without `--junction-quality`
has no `junction_quality_mean` row rather than a zero that reads like a terrible library.
Never: Truncation, a stop codon and an out-of-frame junction are counted **separately**; `_COMPLETE`
folds them together and a QC table must not.
⚠ `allele_candidate` is a **shortlist, not a genotype call** — arda does not genotype. A novel
allele, SHM and a miscall are the same string in the mutation list; recurrence within the allele
(`--allele-min-frac`) and Phred are what separate them, and both are reported per variant.
Never: Chimera / non-functional / stop-codon counts are **flags, never filters**.

## A cohort: `arda qc batch` and `arda qc report`

`arda qc batch -d results/ -o batch [--samples sheet.tsv]` joins every sample's QC table. It reads
**only** the `*.stats.tsv` files — never an AIRR or a clonotype table — so a 1,000-sample cohort
is one `concat` on a laptop against results copied off a cluster. Writes `<prefix>.qc.tsv` (long,
with median/MAD/z), `.qc.wide.tsv` (one row per sample) and `.qc.json`. The sample sheet's optional
`project` and `batch` columns are labels that name the group a sample is compared within; nothing
in the pipeline reads them.

Never: **no threshold is shipped.** There is no pass/fail column and no constant saying what a good
`mapped_fraction` is — it depends on the library, the organism and the depth. Each metric carries
its group's median, MAD and robust z (`0.6745·(x − median)/mad`, flagged at `|z| ≥ 3.5`); a group
under 5 samples gets a median and **no z**, because MAD over four points is not a scale estimate.
Flags, never filters.

Never: **repertoire biology is not here.** Diversity, clonality, rarefaction, overlap and
cross-sample clonotype matching belong to `vdjtools`, which takes arda as a base dependency. arda's
QC answers "did this run work, and is this sample like its batch".

`arda qc report -i batch.qc.json -o batch.qc.html` renders it (or one sample's `.stats.json`) as a
single HTML file with the data inlined and **no external reference of any kind** — no CDN, no
bundle, no new dependency — so it opens air-gapped and survives the results directory.

The two quality columns feeding it are opt-in on `map` and use **different encodings**:
`--junction-quality` writes raw Phred+33 *characters* over `junction`; `--mutation-quality` writes
comma-joined *integers*, one per `v_mutations` / `j_mutations` entry, in the same order.
