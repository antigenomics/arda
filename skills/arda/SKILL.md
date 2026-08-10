---
name: arda
description: >
  Fast TCR/BCR FR/CDR region annotation (Antigen Receptor Domain Annotation).
  Use whenever the user wants to annotate immune-receptor sequences with
  framework/CDR regions, V/D/J gene calls, constant-region isotype (c_call/c_class),
  junction/CDR3 boundaries, or AIRR output — for TCR (TRA/TRB/TRG/TRD) or BCR
  (IGH/IGK/IGL), nucleotide or amino acid, single sequences or large FASTA/FASTQ.
  Also use for: marking up, validating or REPAIRING a bare (CDR3 amino acid, V, J)
  record that has no read behind it — a VDJdb-style row — including locating where it
  disagrees with germline and restoring a missing Cys104 / Phe118 anchor; calling the D
  gene (and tandem D-D) on a bare junction, or inferring it from junction length when the
  protein shows none of it; extracting the receptor repertoire (clonotypes, isotype usage)
  from bulk RNA-seq; getting germline FR1-FR3/CDR1-CDR2 (V) or FR4 (J) subsequences for
  individual alleles; building or rebuilding the reference database from IMGT
  germlines via IgBLAST; or diagnosing mmseqs2 setup. arda runs MMseqs2 + a C++ coordinate
  projection (IgBLAST is offline/build-time only), so it is much faster than IgBLAST at
  annotation time. Load references/ files for the detailed API, region/junction
  semantics, reference-build pipeline, or mmseqs install/troubleshooting.
license: GPL-3.0
compatibility: >
  Python 3.10+; `pip install arda-mapper` (>=2.5.0 for junction markup / repair, the aa D
  posterior, and D on protein input) — no source checkout and **no `ARDA_HOME`**: the curated
  `vdj/` reference auto-fetches into `~/.cache/arda` on first use
  (set `ARDA_NO_AUTO_FETCH` for air-gapped runs), and the `mmseqs` binary auto-fetches a
  static build into the cache if missing — so a bare `pip install` annotates out of the box.
  Bulk RNA-seq needs nothing extra: `seqtree` is a core dependency since 2.5.5. (Before that it
  lived in an optional `[rnaseq]` extra, and a plain install would map and assemble a whole sample
  and only then die, before writing any clonotype table.)
  A source checkout / `$ARDA_HOME` still uses the committed `database/`. Shell is fish — use
  fish syntax in terminal commands.
metadata:
  repo: https://github.com/antigenomics/arda
---

# arda Skills Guide

arda annotates the framework (FR1–FR4) and complementarity-determining (CDR1–CDR3)
regions of TCR/BCR sequences. The expensive IgBLAST markup is done **once,
offline**, when the reference database is built; at annotation time arda only runs
an MMseqs2 search + a C++ routine that projects the reference region coordinates
onto each query. That makes it embeddable and ~4–8× faster than IgBLAST, with 98–99.7%
region concordance on real GenBank mRNA across all five organisms.

It also handles records with **no read behind them** — a CDR3 amino acid plus a V and J
call, as in VDJdb — marking up which residues each germline templates, repairing the
junction, and inferring the D gene from the junction's length.

## Core API

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],  # raw strings or (id, seq) pairs
    seqtype="nt",          # "nt" or "aa"
    organism="human",      # human | mouse | rat | rabbit | rhesus_monkey
    map_d=True,            # map D segments for VDJ loci — works on aa input too
)
# -> list of AIRR record dicts (one per query)
```

For explicit control of strand / sensitivity / in-memory vs file streaming, use
the mapper directly:

```python
from arda.annotate.mapper import annotate_records, annotate_file

recs = annotate_records(queries, organism="human", seqtype="nt",
                        strand="forward", map_d=False, sensitivity=7.0)
annotate_file("reads.fastq.gz", "out.airr.tsv", organism="human")  # streamed, memory-flat
```

Each record dict carries (1-based closed coords, query space): `locus`,
`v_call`/`d_call`/`d2_call`/`j_call`, the constant-region `c_call`/`c_class`
(isotype), `productive`/`stop_codon`/`vj_in_frame`, `rev_comp`, `v_identity`,
`sequence_alignment`/`germline_alignment`, `{v,j,c,d}_cigar`, `*_germline_start/end`,
`v_sequence_end`, `j_sequence_start`, `np1/np2/np3`, `d_support`/`d2_support`,
`junction(_aa)`, and per region in `(fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4)`:
`{r}_start`, `{r}_end`, `{r}`, `{r}_aa`. Ambiguous D and C calls are comma-joined allele
lists (as V/J are). The TSV is a **spec-valid AIRR Rearrangement** file (passes
`airr.schema` validation).

`d_support` is the Karlin–Altschul **E-value the D call was gated on** (accepted at
`<= 0.2` for nt, `<= 0.05` for aa). It ships so a consumer can re-threshold: keeping rows
with `d_support <= x` for `x < 0.2` reproduces exactly what a stricter arda would have
called. A missing `d_call` on a VDJ locus usually means the best hit did not clear the
gate, not that mapping was skipped.

Read [references/annotation.md](references/annotation.md) for the full field list,
parameter semantics (strand/sensitivity/threads/chunking), AIRR column order, the D
E-value gate + genomic-order constraint, and performance notes.

## Batch annotation — never loop (use mmseqs2's own parallelism)

**Always gather every sequence first, make ONE `annotate_sequences` call, then do
downstream analysis on the batch output.** Each `annotate_*` call pays a fixed ~825ms
mmseqs2 process+index-load cost; a batch of 300 sequences costs the *same* ~930ms total
because mmseqs2 parallelises internally across threads. So:

```python
# RIGHT — one batched call, mmseqs2 threads internally
recs = arda.annotate_sequences([(cid, seq) for cid, seq in all_chains], organism="human")
by_id = {r["sequence_id"]: r for r in recs}     # then map back per-item, downstream
```

Do **not** wrap per-item `annotate_*` in a Python `ProcessPoolExecutor`/`ThreadPoolExecutor`
or a loop: a process pool that forks after mmseqs2/BLAS have spawned threads **deadlocks**,
a thread pool just serialises on the same overhead, and either way you pay the fixed cost N
times instead of once. mmseqs2 is the parallel layer — Python orchestration is single-call.

## Region & junction semantics

- Region coordinates are projected through the MMseqs2 alignment, so they are
  correct even for truncated, mutated, or reverse-strand queries.
- There is **no coverage filter**: a partial read (or a bare germline V or J)
  maps to its scaffold and returns only the regions inside its coverage. A bare
  V → `fwr1..fwr3`; a bare J → `fwr4`. This is how callers get per-allele
  germline FR/CDR subsequences without synthesising a rearrangement.
- `junction` spans Cys104 through the [FW]118 that opens FR4; `cdr3` is
  J-anchored. Out-of-frame junctions are reported with an N-bridge (`_`).

Read [references/region-segments.md](references/region-segments.md) for the
bare-germline recipe, junction/CDR3 details, and coordinate round-trip rules.

## Bare records — a CDR3 amino acid, a V call, a J call, no read

VDJdb-style rows have nothing to align. The V and J germlines still template a known run of
residues into each end of the junction, and arda ships those per allele.

```python
from arda.cdr3fix import markup_cdr3, markup_records   # markup_records: a whole polars frame
from arda.annotate.dmap import map_d_junction          # D (+ tandem D-D) on a bare nt junction
from arda.dpost import posterior_d                     # D gene + position from junction LENGTH

mk = markup_cdr3("CAIRDDKII", "TRAV12-3*01", "TRAJ30*01", "HomoSapiens")
mk.cdr3_repaired             # 'CAIRDDKIIF'  -- the Phe118 anchor restored
mk.v_end, mk.j_start         # residues templated by V / index of the first J residue
[str(e) for e in mk.errors]  # ["J del@8 missing 'F' d=0"]
mk.good                      # both sides repaired AND both anchors present
mk.to_cdr3fix()              # VDJdb's `cdr3fix` JSON object, key for key
```

CLI: `arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report - [--d-posterior]`.

> **The single biggest correctness trap.** These coordinates are **junction space**: Cys104
> through Phe/Trp118, **both anchors included**. That is what VDJdb's `cdr3` column holds. It
> is **not** arda's `cdr3` field, which excludes both — `junction_aa` is two residues longer
> than `cdr3_aa`. Conflating them silently corrupts every coordinate, and downstream corrupts
> Pgen, clustering and matching.

Repair is deliberately conservative and its two decisions are separate:

- Every germline disagreement is **reported** (side, kind, position, extent, distance from the
  anchor). Only edits *adjacent* to a conserved anchor are **applied**; deeper ones are left
  alone, because there a mismatch is as likely to be the real V/N boundary as a typo.
  `Cdr3Error.applied` is true only when the edit reached `cdr3_repaired`.
- **A repair always lands on a canonical junction.** If the result would not open with Cys104
  and close with Phe/Trp118, it is refused and the submission returned untouched. So `good`
  implies canonical. An allele with no derivable anchor gives `FailedBadSegment` — flagged,
  never guessed.

`posterior_d` infers the D gene *and where it sits* from the junction's nucleotide length,
which pins `insVD + |D surviving| + insDJ`. Shipped for human IGH/TRB/TRD and mouse TRB only
(the pairs with a published generative model); **every other pair returns `None` rather than
guessing** — do not substitute a human proxy.

## Organisms & loci

| Organism | Loci with full markup |
|----------|-----------------------|
| human, mouse | TRA, TRB, TRG, TRD, IGH, IGK, IGL |
| rat, rabbit, rhesus_monkey | IGH, IGK, IGL (IG only) |

VDJ loci (D segments mapped): IGH, TRB, TRD. D-D fusions sought in all three. D mapping
runs on **protein input too**, against each D germline's three translated frames — useful
for IGH (a call on ~36% of real records, agreeing with the nucleotide call on 98% of them),
mostly silent for the TR loci, whose D is too short to survive trimming into protein. On aa
input `d_germline_*` and `d_cigar` stay empty: those offsets index a reading frame, not the
germline.

**Genomic order constrains the call.** TRBD2 lies 3′ of the entire TRBJ1 cluster, and V(D)J
joining deletes the intervening DNA, so a TRBJ1 rearrangement is never assigned TRBD2 — in
any species with that architecture. IGH and TRD place every D 5′ of every J, so nothing is
excluded there.

Constant-region `J + C` scaffolds (isotype `c_call`/`c_class`) are built for every locus
with a CH1 exon in the bundle.

## CLI

```bash
arda info                                   # versions + available references
arda annotate -i reads.fastq.gz -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta -o out.tsv --seqtype aa --no-map-d
arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report -   # bare (CDR3aa, V, J) records
arda map --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv   # bulk RNA-seq: filter receptor reads
arda assemble -i mapped.airr.tsv -o assembled.airr.tsv    # rescue CDR3s no read spans
arda correct -i mapped.airr.tsv -o clones.tsv            # collapse CDR3 errors into clonotypes
arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # BULK mode: map+assemble+correct
arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # AMPLICON mode: same, other preset
arda shm -i mapped.airr.tsv -o rescoped.airr.tsv                # recount SHM outside the junction
arda igblast -i reads.fastq -o truth.airr.tsv                   # gold-standard IgBLAST (all loci)
arda build-db --organism all                # offline reference build (needs IgBLAST)
arda build-index --organism all             # rebuild mmseqs indexes for local mmseqs version
arda slurm -i big.fastq -o big.airr.tsv --shards 50   # multi-node AMPLICON: split → array → merge
arda cluster submit --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE --shards 8   # multi-node RNA-SEQ
```

## Cluster: two adapters, and picking the wrong one fails silently

| input | command | shard unit |
|---|---|---|
| amplicon / single-end FASTA | `arda slurm` (`arda split` + `arda merge`) | one record |
| bulk paired RNA-seq | `arda cluster submit` (`arda cluster split` + `cluster reduce`) | one read **pair** |

**Never point `arda split` / `arda slurm` at paired RNA-seq.** They write FASTA — dropping the
quality strings `--reconstruct` needs — and round-robin *records*, which puts a fragment's two
mates in different shards. There is no error; the numbers just come out wrong.

**Never shard Stage 2 or Stage 3.** `correct` counts distinct fragments and collapses error
variants globally; `assemble` grows contigs across reads. Per shard, a clone split across N
shards is counted N times and the long-CDR3 contigs Stage 3 exists for are never built, because
the reads that tile them never meet. `arda cluster submit` distributes only `map` and runs the
rest once, through the same function the mode commands use — so a sharded run is
**byte-identical** to a single-node one (verified on real data, all three artifacts).

## mmseqs: nothing to install, and nothing to pin

`pip install arda-mapper` auto-fetches a static binary on first use — there is no extra to
add, and the `[mmseqs]` extra installs nothing (see installation.rst). Candidates are
**version-matched** against the shipped indexes — an index is only reusable by the release that
built it, so an unrelated `mmseqs` on `PATH` would silently discard `database/`'s precompiled
DBs and rebuild a private cache. `$ARDA_MMSEQS` overrides everything and is not checked.

## Run reports: `peak_rss_mb` is monotone, by design

Every stage reports `wall_seconds`, `peak_rss_mb` and `rss_gain_mb`. `peak_rss_mb` is the
**whole-process** high-water mark *as of that stage's end* (getrusage offers no per-stage
reset), which is exactly what a SLURM `--mem` or Nextflow `memory` directive must cover;
`rss_gain_mb` is that stage's contribution. **Budget for Stage 3, not Stage 1**: mapping is flat
at ~300–650 MB at any depth, but the clone set scales with repertoire richness — Stage-3
`correct` peaked at **2,071.7 MB** on a B-cell-rich tumour (28,444 clonotypes from 105 M reads),
versus **549 MB** for a colder sample with *more* reads (139 M). Budget ~4 GB.

## Bulk RNA-seq mode (`arda rnaseq`)

For libraries where only 1–5% of reads are receptor-derived. Three stages, run separately or
in one shot with `rnaseq run` (which does all three by default). Needs the `rnaseq` extra:
`pip install arda-mapper`.

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

## mmseqs2 (auto-installed)

Annotation needs the `mmseqs` binary. Resolution order: `$ARDA_MMSEQS` →
`<project>/bin/mmseqs` → `mmseqs` on PATH → **auto-fetch** a static binary into
`bin/mmseqs`. So neither conda nor pip users must install it manually. The conda
env (`environment.yml`) also ships `mmseqs2` from bioconda.

Read [references/install-mmseqs.md](references/install-mmseqs.md) for env vars
(`ARDA_MMSEQS`, `ARDA_MMSEQS_ASSET`, `ARDA_NO_AUTO_FETCH`), the shipped/precompiled
indexes, and version-mismatch handling.

## Rebuilding the reference

Most users never build anything — `database/vdj/<organism>/` ships with
precompiled markup and MMseqs2 indexes. Rebuild only when adding/refreshing an
organism (needs IgBLAST).

IgBLAST needs no setup step. `setup.sh` puts a release in `bin/`; every other install —
including a plain `pip install` — fetches one into `<cache>/igblast` on first use. So
`arda igblast`, the gold standard every benchmark is scored against, works out of the box.
`$ARDA_IGBLAST` reuses an existing install, and `arda.igblast.igblast_version()` reports which
NCBI release is in play, which belongs in any results record.

Every build writes `loci_manifest.tsv` — one row per defined locus (V-J / J+C
scaffold counts, D germlines, unreachable-D count, `ok`/`EMPTY` status) — and warns at
build end on any `EMPTY` locus or unreachable D germlines. This is what makes an absent
reference visible rather than silent: rat/rabbit/rhesus have no TR loci in IMGT, so
their TCR loci build `EMPTY` (the IG-only limitation in the table above).

Read [references/reference-build.md](references/reference-build.md) for the
`arda.refbuild` pipeline (IMGT germlines → V×J scaffolds → IgBLAST → markup TSVs)
and `build-db` / `build-index`.

### `--prefilter`: the bulk lever

`arda map --prefilter` drops reads that share no exact 16-mer with the reference **before**
they reach MMseqs2. On bulk RNA-seq the receptor fraction is 0.02–3 %, so `mmseqs search` spends
essentially all of its time proving reads are *not* receptor reads; the fitted cost model says so
directly — `wall ≈ reads/46,353 + hits/350`, dominated by the **read count**, not the answer.

| sample | receptor % | before | after | speedup | vs TRUST4 |
|---|---|---|---|---|---|
| SRR10611239 | 0.024 | 82.58 s | 7.79 s | **10.6×** | **0.332** |
| SRR6926533 | 0.123 | 39.11 s | 5.26 s | **7.4×** | 0.603 |
| SRR8363894 | 0.772 | 46.86 s | 13.19 s | **3.6×** | 1.149 |

**Off by default**, because it costs ~0.5 % of real reads. The shape of that loss is what makes it
usable: every lost read scored **75–79 bits** against the `--min-score 75` cutoff (nothing
confident is ever dropped), `junction_aa` moved on **zero** shared reads, and the loss is
**entirely IG with zero across all four TR loci** — one substitution destroys k consecutive exact
windows, and SHM supplies substitutions. Check `prefilter_stats` (`seen`/`passed`) in the run
report to see whether it earned its keep on a given library.

**Measured against an independent IgBLAST truth on 344,554 real IGH-amplicon mates**, the loss is
now on the axis that causes it rather than a locus label: recall **.99209** vs the default's .99582,
and **1,283 of the 1,285 lost reads are in the `<90 %` V-identity bin** — 0 above 95 %, 2 in
90–95 %. On the cell lines the two reads it drops from Raji sit at 89.66 % and 77.78 % identity
(mean identity lost 83.72 vs 93.01 kept), and the 77.78 % read is the most mutated in the sample.
So it is safe for germline-proximal work and progressively costly with SHM load.

⚠ **Bulk only, and that is now quantitative.** Amplicon runs 46–90 % receptor, so almost everything
passes and the scan is pure overhead: on a 90 %-receptor IGH amplicon it is worth **exactly 1.00×**
(319.15 s vs 319.74 s) while still losing those 1,285 reads. Its value is confined to cold bulk
(**3.17×** on a 0.15 %-receptor library), where IG content is negligible anyway. ⚠ **The index is built from `Reference.target_fasta`, the same FASTA MMseqs2 searches.**
Never build it from a hand-listed segment set: against a `V+pad+J`-only reference the loss is
16.29 %, **69.27 % of it J→C reads**. Deriving it from the search target makes that unreachable.

⚠ **Do not reach for an MMseqs2 flag instead.** A 15-setting sweep found **no lossless candidate
above 1.05×**. `--min-ungapped-score 30` is free and gives *zero* speedup — which is the proof
that the cost is the k-mer stage, not the ungapped extension. And MMseqs2 can only prefilter reads
already in a DB, so the FASTA write and `createdb` (19.6 % of a 4 M-read run) are paid regardless.

⚠ **`--adaptive` and `--two-pass` are both WORSE on top of `--prefilter`** (measured: 11.13/16.01/
8.19 s and 14.69/22.84/9.17 s against 10.20/14.66/6.42). Once the scan term is gone, `--adaptive`'s
re-search is pure overhead, and prefiltered survivors being 84 % hitting does *not* put them in
`--two-pass`'s amplicon regime.

### The segment reference and the rescue guarantee

`build-index` also writes `segments.fasta` + `segments.markup.tsv` (`arda.refbuild.segments`):
every V allele, every J allele and every **C** allele as its own target — **924 human targets
(775 V + 124 J + 25 C) against 15,414 V×J scaffolds**, 16.7× fewer. It is *derived* from
`alleles.fasta` + `markup.tsv` in under a second, so it is neither committed nor shipped in the
release tarball, exactly like the mmseqs indexes.

The two-pass search uses it: hit the segment reference, take each read's best V and best J, look
the pair up in `combinations.tsv` — that names exactly **one** V×J scaffold, so the second
alignment is one target per read instead of ~277. A read running J into C names its J+C scaffold
the same way, through `Reference.jc_combinations()` — `(j_call, c_call) → scaffold id`.

⚠ **The constant region is its own target, and that was not always so.** Through 2.7.2 the 345 J+C
scaffolds were copied through verbatim, and they are a J×C product (IGH 14 J × 11 C, IGL 9 × 7,
TRB 16 × 2) in which every scaffold of a locus ends in the *same* constant sequence. A read
reaching C was therefore aligned against all of them to learn one `c_call`: measured on a TRA
amplicon, **27.7 % of the targets drew 76.4 % of the alignments**, 4,977 per target against 603 for
a V. Collapsing it (345 → 25) is **1.89× on the segment search and 1.33× on the whole `map`**, and
it *improves* the calls — a J call decided by a whole-scaffold bit score whose constant half is
arbitrary is a worse J call (`v_call` disagreements vs the one-pass 147 → 85, `j_call` 401 → 296).
A C target is kept for **every** locus, not only the informative ones: only IGH's 11 alleles
separate anything reportable (7 classes = the isotype), but a C target is also the sole segment
target a read lying wholly inside the constant region can hit, and dropping the other 14 makes 14
of 453 fixture reads vanish.

⚠ **`--two-pass` pays off when reads SPAN V INTO J — not by library type.** It needs a read to hit
both a V and a J segment, and the predictor is **`fast_fraction` in the report**, nothing else:

| library | receptor | fast path | vs one-pass |
|---|---|---|---|
| TCR amplicon (TRA) | 48 % | 85 % | **3.51×** |
| 5′-RACE human TRB | 100 % | **95.6 %** | **2.96×** |
| 5′-RACE mouse TRA | 100 % | 89 % | **2.64×** |
| 5′-RACE human **IGH** | 100 % | **16.3 %** | **1.03× slower** |
| IGH RepSeq amplicon | 90 % | **5.2 %** | **0.89× slower** |
| bulk RNA-seq | 2.74 % | 5 % | **0.762×, 31 % slower** |

⛔ This was documented as "an amplicon optimisation, do not reach for it on bulk" and that framing
is **wrong in both directions**. The two 100 %-receptor rows are the *same dataset, same tool*: IGH
reads there cover V and stop short of the short IGHJ target, so they hit a V and no J and the fast
path collapses, while TRB reads span the junction. Read as amplicon-only, the flag gets left off
exactly where it is worth 3×. Bulk is separately a *scan-term* problem, which is `--prefilter`'s
job; this lever only touches the align term.

⛔ **Every `fast path` number above is MMseqs2-specific.** `fast_fraction` is a property of
**(reads × segment mapper)**, not of the reads. On the same 100,000 IGH RepSeq pairs it is **0.052**
with the mmseqs segment search and **0.5018** with `--fast-segments` — 9.7×, from changing nothing
but the mapper, with `v_only` rescues falling 169,004 → 85,933. MMseqs2 misses the short IGHJ on
95 % of these reads and `_segmap`'s ungapped extension finds it on half. Re-read this table with
`--fast-segments` on before concluding a library is a bad fit for the fast path.

Off by default because there is no library type that predicts it — run a sample and read
`fast_fraction`.

### `--fast-segments`: the segment pass without a homology search

`arda map --two-pass --fast-segments` replaces the segment `mmseqs search` with
`arda._segmap`, a C++ seed → vote-by-diagonal → ungapped-extension mapper over the same
`segments.fasta`. It only **nominates**: every candidate is still aligned against the full
`V+pad+J` scaffold and scored by MMseqs2, so this is not a second aligner.

Why it can be exact: the segment pass asks only for the best V allele, the best J allele and their
coordinates (`_SEGMENT_FORMAT`) over a fixed 236 kb germline reference — **no cigar, no backtrace,
no gaps**. Germline V/J carry no indels relative to a read except sequencing error and IG SHM, so
ungapped extension is sufficient, and matches are near-identical to germline rather than remote
homologs.

| | segment step, 100 k amplicon reads, 8 threads | agreement with `mmseqs search` |
|---|---|---|
| `mmseqs search` | 2,770 ms | — |
| `_segmap` | **74 ms** | V allele .9997, J allele .9998, C allele 1.0000 |

End to end it is **1.49×** on `--two-pass` (50 k pairs, 9.10 s → 6.12 s) with `locus` identical on
every read, and on the 13-dataset benchmark panel **1.22× the default at identical recall, with
FP 62 → 59 and zero control FP** — better on every accuracy axis measured there.

**Where it is worth the most: IGH amplicon, 1.87× at 41 % less memory** (100 k pairs, 90 % receptor,
319.74 s → 170.77 s, RSS 4,016 → 2,382 MB), because it raises `fast_fraction` 0.052 → 0.5018 there.
It also **reverses** `--two-pass`'s bulk penalty: 0.73× → **1.24×** on a 0.15 %-receptor library,
where `no_segment_hit` is 99.87 % and those reads skip the full search entirely. `--two-pass` has
always been a structure-preserving prefilter that keeps the V/J call; it just was not worth paying
an `mmseqs search` for on bulk.

⚠ **Cost, measured against an independent IgBLAST truth on 344,554 real mates** (two IGH RepSeq
amplicons): recall **.99525** vs the default's .99582 — **197 reads**, of which **200 are in the
`<90 %` V-identity bin and 0 are above 90 %**. Those losses are *not* a seeding failure: `_segmap`
seeds more reads than mmseqs, not fewer (`no_segment_hit` 395,601 vs 395,682). It promotes ~100 more
reads onto the fast path, so a hypermutated V-only read gets seated on an implied V×J scaffold via a
weak J hit and the alignment then falls under `--min-score 75`. Same class as the round-5 narrowing
refutation — **do not "fix" it by lowering k.** Still off by default for that reason.

Requires the extension — check `arda.segmap.available()`. Without it the flag is a silent no-op, so
assert it in any job that claims to measure it.

⛔ **`_segmap` CANNOT rank `V×J` scaffolds — only segments.** Pointed at the 15,414-scaffold
reference it indexes and maps fine and calls garbage (`v_gene` agreement **.3430**): a
junction-spanning read sits on a scaffold at **two** diagonals, V at one offset and J shifted by the
N-pad plus the non-templated junction, so one ungapped extension scores `max(V, J)` and never their
sum. The true scaffold therefore does not outscore a wrong one sharing its better-covered half.
That is the same fact that makes the segment reference work, read backwards — and it is why the
rescue stays MMseqs2's job. Do **not** retry it with chaining or a banded aligner "to fix the
diagonal": chaining two anchors across a non-templated insert is precisely the general-aligner cost
the whole approach exists to avoid.

Two constants make it equivalent, and both were wrong first:

- **`K = 12`**, the `-k` arda already passes MMseqs2 — *not* 16, which is `prefilter`'s value and is
  calibrated for **rejection**. Seed length sets sensitivity to mismatches: at k=16 the mapper
  seeded 53,048 reads against mmseqs' 53,121, and a read with no segment hit is assumed
  non-receptor and **never rescued**, so those 73 were lost outright.
- **`MIN_SCORE = 40`.** mmseqs applies `-e 1e-3`; this scheme has no e-value, so with no floor
  43,010 reads pick up a constant-region hit against mmseqs' 473 — half of them scoring exactly 38,
  a bare seed plus a couple of flanking matches.

### `--indel-rescue`: what an ungapped extension structurally cannot score

`arda map --two-pass --fast-segments --indel-rescue`. An ungapped extension follows ONE
diagonal, so a read carrying an indel relative to germline scores only up to the indel — measured
in the unit test: a 120 nt read scores **240** clean and **120**, exactly half, with one base
deleted. Its scaffold is then chosen on truncated evidence.

Measured on **341,294 real IGH mates** (IgBLAST gapped-alignment strings; a `-` in either row *is*
an indel, so nothing is inferred):

| V identity | `>=98` | `95-98` | `90-95` | `<90` |
|---|---|---|---|---|
| reads carrying a V indel | 0.74 % | 1.63 % | 3.56 % | **8.00 %** |

3.18 % pooled, tracking SHM load because AID makes indels and not only substitutions.

The signature is in the seed votes before any extension runs — two well-supported diagonals on the
**same** target, offset by the indel length — and the votes are already sorted by
`(target, diagonal)`, so one pass reads it.

**Rerouted, never dropped.** A flagged read is demoted from `implied` to `rescue` and realigned
gapped, so a false positive costs a little speed and *cannot* cost a read. That asymmetry is why
`MIN_DIAG_SEEDS` is deliberately low. Counted as `indel_rescued` in the report.

⚠ **Recall cannot measure this flag.** Recall is identical with and without it (174,066 of 174,226
amplicon fragments either way) because these reads are found regardless. What moves is the CALL,
which is load-bearing: the clonotype key is `(locus, v_call, j_call, junction)` at allele level.
Adjudicated against IgBLAST at gene level on exactly the reads whose call moved:

| sample | v_call moved | correct without | correct with |
|---|---|---|---|
| IGH_repertoire (hypermutated) | 586 | 53.24 % | **84.13 %** |
| IGH_naive (low SHM) | 100 | **77.00 %** | 63.00 % |
| **pooled** | **686** | 56.71 % | **81.05 %** |

**+167 reads, +24.3 points** — and the sign flip is the mechanism confirming itself. Where real
indels are common the flag fixes truncated calls; where they are rare its false positives (repeats
read as two diagonals) dominate. So **turn it on for hypermutated IG and leave it off elsewhere**:
on 13 bulk RNA-seq datasets it demotes **zero** reads and the output is byte-identical, which is
correct — bulk TR carries ~0 indels.

**Reads are never dropped by the fast path.** `arda.annotate.shortlist.shortlist()` partitions
every read into `implied` (took the fast path) or `rescue` (goes back to the full reference), and
asserts the partition is total. Anything that does not resolve — V only, J only, a V×J pair the
reference does not contain, a failed second alignment — is realigned, not discarded.

Traps, each of which produced *correct output* while quietly breaking something:

- **`JC|` targets are named by scaffold id, not by allele.** Feed the raw target name to
  `shortlist()` and every J→C read returns `no_such_combination`. Resolve through
  `Reference.segment_j_call()`. Measured cost of getting this wrong: the fast path collapsed
  from 85.3 % to 0.1 %. (Only reachable on a pre-2.8.0 reference; `C|` targets are named by allele.)
- **Never `top_hit` the segment pass.** One best hit per read destroys the V+J pairing the whole
  scheme depends on — `implied` goes to 0.
- **`_SEGMENT_SIDE` is ONE mapping shared by `_segment_rows` and `_segment_best_hits`.** Spell the
  target-kind rule out twice — once in polars, once in Python — and they drift: when `C|` was added
  to the loop alone, the polars reduction discarded every C row, `best_c` stayed empty, no
  constant-only read was rescued, and 15 J→C reads vanished **without `no_segment_hit` moving**,
  because the rows were gone before any counter saw them.
- **Nominate the J+C contest from the J, never from a C hit.** A J→C read spans the J/C boundary,
  so its constant overlap is often below the search threshold on its own even though the
  concatenated scaffold cleared it. Gating on C evidence re-invented `TRBV12-3*02`, destroyed a
  `c_call`, and fabricated a `junction_aa` on a read the one-pass calls V-less.
- **`segments.fasta` is generated, and a deploy does not regenerate it.** A new mapper reads stale
  `JC|` targets through its back-compat path and *succeeds* while measuring the old reference.
  Assert the target composition, not the file's existence.

**TRD is TRAV/DV + TRDJ; the J (and C) decides the locus, not the V.** arda's reference already
encodes it — of 1,050 scaffolds built from a TRAV/DV segment, 1,005 are locus TRA (with a TRAJ)
and 45 are locus TRD (with a TRDJ). There is deliberately no "α/δ is ambiguous" rule:
`combinations.tsv` is the arbiter, a pair it contains is real biology and one it does not is a
genuine chimera.

### `annotate.project`: the junction from coordinates, no alignment

```python
from arda.annotate.project import project_junction, UNVALIDATED_LOCI

proj, refusal = project_junction(strand_seq, qlen, v_row=v_row, j_row=j_row,
                                 v_call=..., j_call=..., anchors=ref.anchors,
                                 split_checked=...)   # REQUIRED keyword, deliberately no default
if proj:
    proj.junction   # AIRR: Cys104 .. [FW]118 INCLUSIVE
    proj.cdr3       # IMGT: anchors excluded, two residues shorter
```

AIRR `junction` runs Cys104 -> [FW]118. The segment pass already returns `(target, tstart, qstart)`
per read per side and `cdr3_anchors.tsv` already records `anchor_nt`, so the position is arithmetic:

```
offset = (anchor_nt + 1) - tstart
pos    = qstart + offset                    # forward
pos    = (qlen - qstart + 1) + offset       # reverse complement
```

⛔ **Three coordinate systems disagree about their origin here.** `tstart` is **1-based** on the
forward target (`segmap.cpp`), `anchor_nt` is **0-based** in the germline (`cdr3fix.Anchor`), and a
minus-strand `qstart` is in **forward** coordinates while the sequence it indexes is the reverse
complement. Each is an off-by-one that still yields a plausible-looking junction — right length,
starts with a codon, ends with a codon. `strand_seq` must be the strand the hits were MEASURED on;
pass `revcomp(read)` for a minus-strand hit rather than reflecting coordinates afterwards.

Measured vs IgBLAST at `v_score >= 70` on 254,867 reads: IGH .99977 (naive) / .99634 (91.77 % median
V identity), TRA .99947, TRB .99949, IGK and IGL exact. >= .993 in every V-identity stratum.

**It refuses rather than degrades** — `no_anchor`, `unvalidated_locus`, `indel_split`,
`strand_mismatch`, `order`, `off_read`, `bad_codon`. A well-formed junction that is wrong is worse
than no junction: the reference-geometry bug shipped junctions that started `C`, ended `[FW]`,
passed `--complete-only`, and were short by exactly the allele's truncation.

⛔ **`UNVALIDATED_LOCI = {"TRD"}`.** Not because TRD is known bad — because it has **zero** coverage.
Across two TR amplicons the segment pass never handed a TRD read both anchors, so all 767 TRD
truth junctions went to the aligner and TRD never appeared. Absent reads like fine in every
aggregate. The locus is taken from the **J** anchor, never the V: TRAV/DV rearranges to either TRAJ
or TRDJ and the J decides.

⛔ **A `[FW]GXG` motif check is NOT equivalent to reading `anchor_nt`.** `TRAJ35*01`'s anchor codon
decodes **Cys (TGC)** — it is `status = ok` and a functional IMGT `F` gene, with a real Cys six
codons past the FGXG. A motif check deletes the gene silently (33/33 amplicon reads lost).

⚠ Yield, not accuracy, sets the scope: 87.0 % / 77.3 % of hit amplicon reads carry both anchors
against **7.2 % of bulk** reads, which mostly do not span a junction at all.

## Sequence primitives

`arda.refbuild.translate` exposes fast C++-backed helpers, mirpy-API-compatible:
`translate(nt, frame=0)`, `detect_coding_frame(nt)`, `reverse_complement(nt)`,
`back_translate(aa)`, `aa_coords_from_nt(nt_start, nt_end, coding_start)`.

## Gotchas

- **`junction` is not `cdr3`.** `junction`/`junction_aa` include both conserved anchors;
  `cdr3`/`cdr3_aa` exclude both, so `cdr3_aa == junction_aa[1:-1]` always. Everything in
  `arda.cdr3fix` / `dmap` / `dpost` works in *junction* space, matching VDJdb's `cdr3`
  column. Mixing the two conventions is the most expensive mistake available here.
- **An empty `d_call` is a decision, not a gap.** The call is gated on `d_support` (E-value
  ≤ 0.2 nt, ≤ 0.05 aa). Human TRB gets a D on only ~47% of junctions because an ordinary
  TRB interior is 11–21 nt and heavily-trimmed TRBD1 scores below the gate. VJ loci (TRA,
  TRG, IGK, IGL) have no D gene at all.
- aa input returns region `*_aa` directly with no frame bridging, so `stop_codon` and
  `vj_in_frame` stay empty — but `productive` and the D columns *are* populated.
- `posterior_d` returns `None` for organisms with no shipped generative model (rat, rabbit,
  rhesus) and for VJ loci. That is deliberate. Do not fall back to a human model.
- The shipped MMseqs2 indexes are used only when the local mmseqs **version**
  matches; otherwise arda rebuilds a private cache in `data/` on first run.
  `arda build-index` (re)builds the shipped indexes for your version.
- `map_d=True` on synthetic/partial input with no real junction simply finds no
  D — harmless; pass `map_d=False` to skip the search.
- IgBLAST is needed only to build references, never at annotation time.
- `arda correct` uses `seqtree` (a core dependency since 2.5.5);
  without it every `correct` test **skips silently** rather than failing.
