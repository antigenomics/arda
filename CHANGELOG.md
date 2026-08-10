# Changelog

Notable changes per release. Earlier releases are described by their git tags
(`git tag --sort=-v:refname`); this file starts at 2.5.0.

## 2.15.0

### Fixed — contig assembly (⛔ read counts and contig sequences change)

A second audit went over SHM calling and contig assembly. Eight defects, each with a regression test
verified to fail without it.

**⛔ Contigs were READ-ORDER DEPENDENT.** The extension tie-break was `len(ext) > len(best_ext)` with
a strict `>`, so equal-length candidates were resolved by the order the posting list happened to be
in — which is AIRR row order, which comes from a threaded MMseqs2 search. The contig *sequence*, and
every junction derived from it, could differ between runs on the same input. Now ordered on
`(length, sequence)`, a total order, in both the 3′ and 5′ passes.

**A contig's junction is now only attributed to members that COVERED it.** Membership is granted on
a `min_overlap` match, and once the extension passes have accumulated germline at the contig ends
that match can be *pure germline V* — the 5′ pass says so itself: "that region is shared germline, so
any V-read of the gene extends it correctly". Every incomplete member was nonetheless stamped with
the contig's junction, so a read of a **different clone of the same V gene** was credited to this
clonotype's `duplicate_count` on no clone-specific evidence. Members now carry their span in contig
coordinates and must cover the junction by ≥10 nt.

⚠ **This reduces reported abundance.** On the committed 1,035-read example: rescued members 162 →
145, reads counted 384 → 367, clonotypes 19 → 18. Those 17 reads are *genuinely uncounted, not
moved* — a germline-only read matches no junction, so coverage assignment cannot re-place it. It is a
deliberate precision-over-recall trade, unlike the error-correction path where read conservation is
an invariant. ✅ The clonotype that went is the case *for* the change: `CAASMAGGGNKLTF` under
`TRAV13-1*01` with 2 reads against the same junction under `TRAV13-1*02` with 28 — two alleles of one
gene, i.e. de-duplication rather than loss.

**A rejected contig now releases its reads.** `used` was set as reads were recruited but a contig
dropped for having <2 members never gave them back, so a seed that failed to extend was permanently
consumed and could not join a later contig even as an ordinary extension member. Seeds are tried
longest-CDR3-tail first, so this stranded exactly the short-tailed reads that most need a contig.

**The assembly k-mer index is bounded at insert.** It posted every k-mer position of every mapped
read of the locus with no cap; `scan_cap` bounded only how many postings were *consumed*, and
`--assemble` is on by default. ⚠ Proved equivalent before shipping, since a cap that changes the
candidate set is a behaviour change and not an optimisation: both consumers already slice
`[:scan_cap]`. Measured on 20,000 reads sharing a 60 nt germline prefix — candidate sets identical
for every k-mer, postings 1,600,000 → 784,000, ~112 MB → ~82 MB.

### Fixed — clonotype reporting

**Isotype is one vote per FRAGMENT, at last.** `_dominant_ccall` deduplicated its read list down to
fragments and then re-expanded to one entry *per row*, so a fragment whose two mates both carried a
`c_class` voted twice, and an assembly-rescued fragment voted again. A one-fragment minority could
outvote a two-fragment majority. ⛔ The tally was also order-dependent — `Counter.most_common(1)`
breaks ties by insertion order — so a two-way isotype tie could report a different class run to run.
Now lexicographic.

**Under `--all-junctions`, an assembled row outranks the read's own truncated junction**, and a
clonotype left with no reads is not emitted. Coverage pass 1 walks the concatenated frame mapped-rows
first, so the read's own *truncated* junction won the race and the contig's clonotype was emitted with
`duplicate_count` 0 — measured on a 2-read fixture: `CAKGALQ` dup=2 beside `CAKGALQKW` dup=0. Neither
half touches read conservation: the reads move, and a dropped row has none by definition.

**`reads_with_junction` no longer double-counts assembly-rescued reads** (once for the incomplete
Stage-1 row, once for the assembled row) while being documented as a Stage-1 statistic.
`reads_from_assembly` and `empty_clonotypes` report the two new quantities separately.

### Fixed — `arda resolve-ties` degraded silently without a reference

It caught `OSError` per locus, continued with an empty germline set, and returned output
byte-identical to its input with `expanded: 0` — indistinguishable from a library that genuinely had
no ties. The same "raise, never degrade" rule `--min-junction-q` and the quality rescue already
follow. ⚠ This had kept CI red since tie lists were added.

### Added — the junction boundary in GERMLINE coordinates

`v_anchor_nt` and `j_anchor_nt` (0-based Cys104 / [FW]118 offsets in the called allele) are emitted
per read, appended last. They existed only inside the reference before, so a consumer holding the TSV
could not tell a framework mutation from a junction one. Now:

* a `v_mutations` entry at 1-based `p` is junction-internal iff `p > v_anchor_nt`;
* a `j_mutations` entry at `p` is junction-internal iff `p <= j_anchor_nt + 3`.

This is what makes allele re-assignment and SHM correction a *downstream* job — see
`docs/shm.rst` for the framework-only identity recipe and the measured split.

### Changed

**`scripts/pack_reference.sh` rejects an empty reference file**, not just a missing one. A
`build-db` run without IgBLAST emits 0-byte `alleles.*`/`markup.*`, which passed the presence check,
shipped, and crashed every fresh install with `NoDataError: empty CSV` — v2.5.7 shipped exactly that.
The smallest genuine file is 24,372 B against a 1,000 B floor (#83).

**Correction to the 2.14.0 notes:** the coverage k-mer cap's under-assignment is **library-dependent**,
not a general property. It tracks how much germline the junctions share: **1.6 %** on a TRA amplicon
(short, near-germline junctions over 36,741 roots) against **0.013 %** on a TRB amplicon of comparable
depth (72,339 roots, D+N-bearing junctions that are far more specific).

### Known defects, documented not fixed

`v_mutations` / `j_mutations` **and** `v_identity` include junction-internal positions — the V
germline's 3′ tail and the J germline's 5′ head lie inside the junction, and both are scoped by
segment rather than by the junction boundary. Measured on a TRA amplicon, where TCRs cannot
hypermutate so every entry is spurious by construction: **1.046 V and 1.658 J entries per read**, and
`v_identity` reporting **0.8723 on a TCR** whose framework-only identity is **1.0000**. ⚠ Frequency
alone does not separate error from allele from junction diversity — frequency **and** position against
the anchor does. Fixing it inside arda changes every published SHM number, so it gets its own round;
`docs/shm.rst` carries the recipe and the measurements meanwhile.

## 2.14.0

### Fixed — a full-depth read leak in coverage assignment (⛔ read counts change)

**`--ec-mode accurate|amplicon|rnaseq` lost reads out of the clonotype table**, and the amount
scaled with library depth, so no local test saw it. Measured across the 16-sample golden set at full
depth: **7 samples leaked, −107,440 reads in total**, worst on the TRA amplicons (SRR5233636
−25,468 of 1,838,213, i.e. 1.39 %; SRR5233635 −7,383) and present on **both MIGEC published-truth
libraries** (−755 IGH, −1,351 TCR). The cell lines all *gained* (Jurkat +44, Raji +46,
MouseSpleen +5), which is exactly why a monoclonal QC panel never caught it.

The cause is not the error model. `_assign_coverage` bounds each k-mer's postings at `cap` and fills
them in descending abundance; an **alias** — a junction the quality gate vacated, kept in the index
only so partial reads that covered it still reach the parent — is ordered by its **parent's** count,
which is high by construction. Aliases therefore sorted to the front and evicted genuine
low-abundance roots, and every partial read whose only home was such a root went unassigned. On
SRR5233636, every arm emits an **identical 36,587-row clonotype table**, so this is purely read
assignment:

| index                                | reads assigned | vs `fast` |
|--------------------------------------|---------------:|----------:|
| cap 64, aliases ordered by abundance |      1,812,740 |   −25,473 |
| cap 64, aliases off                  |      1,838,181 |       −32 |
| **cap 64, roots first (this release)** | **1,841,624** | **+3,411** |
| cap 1024, aliases by abundance       |      1,869,556 |   +31,343 |

The alias mechanism was added to rescue 5 reads of 9,208 on Ramos and was costing 25,441 here. It is
worth keeping — the same aliases *gain* 31,343 once the index can hold both — it simply must not
outrank the roots. Every mode now gains rather than loses: on SRR5233636 `accurate` +3,411,
`amplicon` +2,827, `rnaseq` +3,393.

**`CorrectReport.reads_assigned`** is new and is the reason this was invisible. `reads` is the
spanning-read count taken *before* correction, so it cannot move; comparing it across `--ec-mode`
returns 0 on every sample and reads exactly like conservation holding. `reads_assigned` is
`sum(duplicate_count)` over the emitted table — the quantity the invariant is actually defined on.

**The quality rescue no longer merges across loci.** It grouped candidate parents by junction
*length* alone, and `--ec-mode amplicon` opens the radius to 12 substitutions, so on SRR5233636
3 of 9,025 rescues were 1-read **TRB** clonotypes absorbed into abundant **TRA** clonotypes at 11–12
substitutions — misattributed expression, not lost reads. The search now runs once per locus, which
also bounds it: it is quadratic within a length bin and single-threaded, and a 397 k-clonotype
library was still running after 39 minutes of CPU against ~5 for a 49 k one.

⚠ **V and J stay ignored by the rescue, deliberately.** The abundance model defaults to
`--require-vj` because a true error keeps the germline call — true for the 1–3 substitution
neighbours it collapses. The rescue targets the opposite class, where the *whole junction window* is
unreliable (median mean Phred 16.5–20.1), and a read that bad has an unreliable V/J call for the
same reason. 50.9 % of `amplicon` rescues cross V. What protects a genuine clone is the two gates
that *are* trustworthy: its reads must be measurably bad, and the parent must be
`lowq_min_ratio` times more abundant.

**The quality gate's parent must now be able to have PRODUCED the read.** Its only test was
`count(parent) > count(child)` — one extra read made anything within `--max-subs` a parent, and at
3 substitutions nothing supports that. The plausibility test is computed from **this read's own
Phred**: the product of 10^(−Q/10) over the discriminating positions, times the parent's count, must
reach the child's. ⚠ Using the global `--error-rate` instead was tried and is wrong the other way —
it makes the gate a strict subset of the abundance model, which is what the gate exists to reach
past (MIGEC at 1e-5 went 1,633 → 1,633 error clonotypes, from 158).

**`--ec-mode amplicon|rnaseq` now raises without `junction_quality`** instead of silently skipping
the rescue and reporting counters that are indistinguishable from a clean library — the same
"raise, never degrade" rule `--min-junction-q` already followed.

**The Nextflow module emitted a flag `arda rnaseq run` did not have.** `main.nf` appended
`--clonotype-key` to the `run` command line, where only `correct` accepted it, so
`arda_clonotype_key = 'junction'` failed the process outright. `run` now accepts it and plumbs it
through to `correct`.

### Changed

**Orphon V/J genes are excluded from the reference.** IMGT `/OR` genes sit outside their locus
(`TRBV20/OR9-2` is on chromosome 9, TRB is on 7) and cannot rearrange, so they were pure false-call
surface. ⚠ This affects `arda build-db`; an auto-fetched prebuilt reference needs regenerating to
pick it up.

**D germline loading is cached.** `_d_germlines` re-opened and re-parsed `d_germlines.fasta` on
every call and `_clonotype_d` calls it once per clonotype — 54.4 µs each, so **2.0 s** wasted on a
36,741-clonotype amplicon and **21.6 s** on a 397,305-clonotype library, on the default `--map-d`
path. Now 0.056 µs.

### Added

**Tie lists (`--tie-lists`, off by default) and `arda resolve-ties`.** A read aligned over a span
that several germlines share carries no base that separates them, so naming one is a claim the data
does not support. Membership is decided per read from the span it already aligned — a string
comparison against the reference, not a new alignment — then ranked library-wide so the allele the
whole library supports leads, with only unambiguous reads voting.

### ⛔ Known defect, newly measured — `v_mutations` / `j_mutations` include junction positions

`docs/shm.rst` claimed the mutation lists were scoped to V and J **structurally**, so "a junction
position has no germline coordinate to be filed under and cannot enter the list by any code path."
**That claim is false and is retracted in this release.** The N-pad is excluded, but the V
germline's 3' tail and the J germline's 5' head lie *inside* the junction, and mutations are scoped
by segment (`t <= t_vend`, `[t_jstart, t_vjend]`) rather than by the junction boundary — so
exonuclease chew-back and non-templated N/P bases are reported as substitutions against a germline
that does not template them.

Measured on a **TRA amplicon** (SRR5233636, 500,000 reads), where T-cell receptors do not
hypermutate so every entry is spurious by construction: **1.046 V and 1.658 J mutations per read**,
with **86.2 % of J entries at J germline position <= 10**. Splitting the load by the frequency of
each `(allele, position, alt)` across reads carrying that allele separates three superimposed
populations:

| frequency | share of J entries | share of V entries | what it is |
|---|---:|---:|---|
| < 0.01 | 13.0 % | 39.0 % | sequencing error |
| 0.01–0.5 | 80.8 % | 59.2 % | junction diversity misattributed as SHM |
| >= 0.5 | 6.2 % | 1.9 % | a genuinely wrong **allele** call |

The high-frequency tail is real and separate — `TRAV8-6*01` positions 281/282 at 0.88, `TRAJ8*01`
position 1 at 0.67 — i.e. the called allele is not the one the reads carry. `v_identity` has the
same scope defect: it runs to `t_vend`, so it is depressed by junction diversity rather than by
mutation load.

**Not fixed in this release.** The fix needs the scaffold's CDR3 boundary threaded into the C++
markup and the Python reference implementation together, and it changes every published SHM number,
so it gets its own round with its own measurement. Documented in `docs/shm.rst` with the workaround
in the meantime: separate the three populations by frequency.

### Known limitation, measured

The coverage k-mer `cap` under-assigns reads in the default path, independent of the alias fix —
but ⚠ **the size is library-dependent, not a general property**, and the release notes originally
overstated it. It is driven by how much germline the junctions share: on a **TRA** amplicon
(SRR5233636, short near-germline junctions, 36,741 roots sharing k-mers heavily) it is **1.6 %**,
while on a **TRB** amplicon of comparable depth (SRR5233642, 72,339 roots, D+N-bearing junctions
that are far more specific) the same sweep moves **+236 reads at cap 128 and +404 at cap 256 —
0.013 %**. Same biology that made TRA the read leak's worst case.

On SRR5233636 with `--ec-mode fast` (no aliases exist), the clonotype table is identical at every
cap:

| cap | reads assigned | wall |
|----:|---------------:|-----:|
|  64 |      1,838,213 | 207.4 s |
| 128 |      1,844,359 | 312.1 s |
| 256 |      1,850,917 | 485.9 s |
| 512 |      1,867,867 | 740.7 s |
|1024 |      1,867,904 | 865.0 s |

It saturates at 512: the whole prize is +29,654 reads and 1024 adds 37 more. Peak RSS is flat, so it
is CPU in the alignment inner loop. The default stays at 64 — 3.6× on this stage for 1.6 % of
coverage-based abundance is a per-workload judgement. ⛔ Guaranteeing every root a posting does
**not** help: built, measured, byte-identical on an amplicon and a bulk sample at two caps. No root
is ever fully unreachable; the reads are lost because a root's surviving postings sit at k-mer
positions the partial read does not cover.

## 2.13.2

### Added — the denoising framework is reachable from Nextflow, and the cluster path is documented

Audit of "is everything shipped and usable" found two real gaps. Both are integration, not
behaviour: **no arda output changes in this release.**

**The Nextflow module could not reach the framework.** `main.nf`, `nextflow.config` and `meta.yml`
had zero mentions of `--ec-mode`, `--clonotype-key` or `--junction-quality`, so a pipeline user had
no way to it short of `task.ext.args` surgery. Now:

```groovy
params {
    arda_ec_mode       = 'amplicon'   // fast (default) | accurate | amplicon | rnaseq
    arda_clonotype_key = 'junction'   // full (default) | junction
}
```

Both are validated against their allowed sets, and the module **warns** when `arda_ec_mode` and
`regime` disagree — the amplicon and rnaseq presets are tuned to opposite clonotype-size
distributions, so the mismatch is a real cost rather than a no-op. Declared in `nextflow.config`
with the measurement behind each, and read through `params.getOrDefault` so the module stays
correct when included without its config.

**SLURM was implemented but undocumented.** `arda slurm` / `split` / `merge` and
`arda.cluster.render_submit_script` have shipped for releases with no page describing them. New
`docs/cluster.rst`: the one-command chain, why `split` and `split_pairs` are not interchangeable
(one writes FASTA and round-robins *records*, splitting a fragment's mates across shards), why
Stage 2 must run **once globally** rather than per shard, aligner pinning, regime choice, and
resource sizing — including that prefilter threads saturate at 16 and regress at 32.

New `docs/use_cases.rst`: bulk RNA-seq, amplicon, monoclonal QC, negative controls, low-frequency
variants, SHM, and the rules for comparing arda against another tool (name the stage; benchmark
every tool at its best config; compare at gene level; give each call metric its own denominator).
It also states plainly that a 200 k subsample of a bulk library yields 1–3 clonotypes and nothing
computed on it means anything.

README gains the monoclonal-QC and per-regime rows, the read-conservation invariant, and why the
modes are off by default.

## 2.13.1

### Fixed — a Q1 base in `junction_quality` made the whole AIRR unreadable

`junction_quality` is a Phred+33 string, and **chr 34 is `"`, i.e. Q1** — a legitimate score that
any low-quality base produces. polars' CSV reader treats it as a quote character, so **one such
base collapsed the parse of the entire file**:

```
ComputeError: CSV malformed: expected 1 rows, actual 155 rows
```

Found on a real Raji run in the round-23 benchmark: exactly one row of the file contained a `"`,
and every `correct` leg for that sample died. It is not rare — it needs one Q1 base anywhere in one
junction, so any library with a low-quality tail hits it.

Underneath was a second problem: **arda wrote the format two ways.** The streaming writer
(`_markup.format_rows`, which produces every `map` output) emits raw fields and a truly empty
string for a missing value; polars' `write_csv` quotes, rendering an empty string as the two
characters `""`. Reading unquoted is right for the big files and would have turned every empty
field of an older polars-written one into a literal `""`.

So both sides are fixed: every AIRR reader now reads unquoted **and normalises a literal `""` back
to empty** (an AIRR field is never legitimately that two-character string), and arda's polars
writers emit `quote_style="never"` so there is one dialect going forward.

⚠ `examples/rnaseq/clones.tsv` changes accordingly — 18 rows where `""` becomes genuinely empty.
No value changes.

Pinned by `test_a_q1_base_in_the_quality_string_does_not_break_the_parse`.

## 2.13.0

### Added — a quality-aware denoising framework: `--ec-mode amplicon|rnaseq`, `--clonotype-key`

Stage 2's abundance model asks whether the parent could have produced this many misreads by
chance. That question is **only answerable for near neighbours**, and it is now measured exactly
where it stops being answerable. If a clonotype *k* substitutions from its parent is accumulated
independent error, the 1- and 2-substitution intermediates must also be observed — the same process
generates them at far higher rate. On Jurkat TRB (dominant clone 9,932 reads, 48 nt):

| k | clonotypes | expected at `1e-3` | median mean-Q | has an observed *(k−1)* neighbour |
|---|---|---|---|---|
| 1 | 108 | 476.74 | 31.4 | — |
| 2 | 82 | 11.20 | 25.2 | **82 / 82** |
| 3 | 28 | 0.17 | 24.2 | 3 / 28 |
| 4 | 13 | 0.0019 | 24.0 | **0 / 13** |
| 5 | 18 | 0.0000 | 20.1 | **0 / 18** |
| ≥6 | 14 | 0.0000 | 16.5–18.2 | **0 / 14** |

**k ≤ 2 is a ladder** and the abundance model is valid on it — which is what `--max-subs 3` has
always been. **k ≥ 4 is a cliff**: zero intermediates, and 13 observed where the model predicts
0.0019. Those are single bad reads, and quality says so independently (median mean junction Phred
falls monotonically 31.4 → 16.5; the k ≥ 5 class is 100 % sub-Q30 against the dominant clone's
5.9 %).

⛔ Widening `--max-subs` to 10 *does* clean that library up (53 → 11 clonotypes) **for the wrong
reason** — the abundance test it applies there has probability 0 to every printed digit, and
`--error-rate` is inert from 1e-3 to 1e-1 on that class for the same reason. So the new modes reach
it on the evidence that actually distinguishes it:

* `--ec-mode amplicon` — rescue at mean-Q < 25, radius 12 subs, ratio 50×. A targeted library is
  deep, so a 1-read neighbour of an abundant clone is almost always error.
* `--ec-mode rnaseq` — rescue at mean-Q < 20, radius 6 subs, ratio 200×. Bulk RNA-seq is sparse
  (0.02–3 % receptor), singletons are the norm and mostly real, so the rescue stays narrow.

⛔ **Nothing in the framework discards a read.** A candidate with no qualifying parent keeps its
reads and is reported as an orphan. That is not caution, it is the measured requirement: on a
polyclonal hypermutated repertoire a whole-junction mean-Q floor at Q30 strands **3.70 %** of all
junction-bearing reads with no parent to inherit them. Read conservation is pinned over every
regime × key by `test_no_regime_or_key_ever_loses_a_read`.

`--clonotype-key junction` canonicalises V/J to the junction's majority before grouping, so **call
splits collapse** — a junction byte-identical to an abundant clone's under a different V or J call,
which no error model can see because there is no discriminating base. On Jurkat that is the largest
error class by reads (130 of 14,531, including an allele-level TRG split). Cost on a polyclonal TRA
amplicon: 132 of 19,956 clonotypes merge (0.66 %), and the minority call there carries one read
against 4–10 on a short 30–39 nt junction.

Jurkat, `--clonotype-key junction`, 14,531 reads throughout:

| `--ec-mode` | clonotypes | reads | TRB purity | rescued | orphans |
|---|---|---|---|---|---|
| `fast` | 89 | 14,531 | .99562 | 0 | 0 |
| `accurate` | 53 | 14,531 | .99696 | 0 | 0 |
| **`amplicon`** | **10** | **14,531** | **.99990** | 43 | 2 |
| `rnaseq` | 40 | 14,531 | .99800 | 13 | 16 |

MIGEC spike-ins: all three published clonotypes survive every mode (error clonotypes 1,630 → 108
under `amplicon`) at exactly 310,559 reads.

New C++ extension `arda._denoise` (per-read mean junction Phred; the wide-radius, length-bucketed,
early-exit parent search — 31,943 clonotypes is 5.1e8 ordered pairs on an IGH repertoire). Both
have Python reference implementations asserted identical in the tests.

### Fixed — the coverage k-mer cap kept an arbitrary 64 roots, not the abundant ones

`_assign_coverage` bounds how many roots one germline-shared k-mer may name, and kept the first
`cap` in target order — so the survivors depended on the root list's order and any change to the
root set silently reshuffled them (merging call splits cost 3 reads of 43,475, none of which had a
junction of their own to fall back on). Now inserted in descending abundance, then sequence: the
roots a partial read is most likely to belong to are the ones that survive. On the TRA amplicon
this also **recovers 35 reads** that were previously unassigned (43,475 → 43,510).

⚠ Residual, named rather than swept: with the quality gate on, 9 reads of 43,510 (0.021 %) are
still not placed — 7 carrying a frameshift-N junction and 2 with no junction at all, i.e. reads
that never formed a clonotype and were only ever placed by alignment.

## 2.12.1

### Fixed — the quality gate MOVES a read onto its parent, it never discards it

A read that reached a complete junction came off a real rearrangement of that locus.
`--min-junction-q`'s evidence is about **one base** — it says that base is a miscall — and says
nothing about whether the molecule existed, so discarding the read understates the parent's
expression by exactly the reads the correction claims belong to it.

2.12.0 filtered gated reads out of the frame. Coverage assignment then happened to realign most of
them back onto the parent — 325 of 330 on a TRA amplicon — so the rule held *by accident*, resting
on an alignment succeeding, and 5 reads were lost. The gated read's clonotype key is now rewritten
to the parent's and it routes through the exact-key pass. Three leaks had to be closed:

* coverage re-derives every read's clonotype from the **unfiltered** Stage-1 frame by the key still
  on the read, so a moved read whose old clonotype survived (its other reads passed the gate) was
  handed straight back and the move silently did nothing;
* a clonotype the gate empties **vacates its key**, and any other read carrying it — an
  incomplete-junction read, which never reached Stage 2 — loses its only anchor;
* the vacated **junction** left the alignment index, so a partial read whose only ≥20 nt overlap was
  with it went unassigned.

Measured. Jurkat `ERR3003543`, reads assigned **14,531 at every gate** while clonotypes go 90 → 54
and purity TRA .99055 → .99528, TRB .98963 → .99096. MIGEC spike-ins at `--error-rate 1e-5`, reads
assigned **310,559 at every gate**, all three published clonotypes kept, error clonotypes
1,630 → 79: error reads fall 8,398 → 2,026 and the parent gains **+6,545**.

⛔ Scoped to gate-vacated junctions, **not** to every collapsed child. Aliasing all of them was
built and measured: it moves *default* output (Ramos 9,208 → 9,234 with the gate off) and still
loses 14 reads at Q20. Rejected.

## 2.12.0

### Added — `--d-max-evalue`: the D call is a dial, and `d_support` ranks it

The E-value that accepts a D call was a fixed operating point. It is now a knob on `annotate` and
on every `rnaseq` stage that maps D (`map`, `correct`, `assemble`, `run`, `reduce`), and on
`dmap.map_d_junction`. Measured against IgBLAST at gene level, the shipped 0.2 is deliberately the
**weakest** band in both libraries tested — it is the highest-recall setting, which is the right
default for a repertoire tool and the wrong one for a call you will act on:

| `--d-max-evalue` | TRB amplicon: called / rate / agreement | bulk IGH: called / rate / agreement |
|---|---|---|
| **0.20** (shipped) | 18,362 / .4026 / .9765 | 648 / .3610 / .9417 |
| 0.10 | 14,690 / .3221 / .9842 | 557 / .3103 / .9494 |
| 0.05 | 10,749 / .2357 / .9911 | 461 / .2568 / .9746 |
| **0.01** | 5,355 / .1174 / **.9985** | 347 / .1933 / **1.0000** |

(TRB: SRR5233641, 45,604 reads with a projected V..J interior against 31,608 IgBLAST D calls at
`v_score >= 70`. IGH: SRR5233639 at full depth, 1,795 reads with an interior, 1,056 truth calls.)

⛔ The CLI default is `None`, **not** `0.2`. The shipped operating point is alphabet-dependent —
0.2 for nt, 0.05 for aa — so a literal `0.2` would have silently loosened `--seqtype aa` by 4×
while looking like a no-op. Pinned by `test_the_d_evalue_cli_default_does_not_loosen_the_aa_gate`.

⚠ A composition-preserving shuffle null puts the false-D rate at **.0884 in TRB** against a real
call rate of .4025 at the shipped gate (IGH .0095, its germlines being longer). A TRB D at 0.2 is a
ranked hypothesis, not an identification.

### Fixed — a tandem D-D must run in genomic order

D-D fusion is a rearrangement: the upstream D joins the downstream D and everything between is
deleted, so the product carries them in genomic 5′→3′ order. `transfer._dd_orientation_ok` refuses
a pair that does not, over `_D_GENOMIC_ORDER` (`TRBD1 < TRBD2`, `TRDD1 < TRDD2 < TRDD3` — the
architectures pinned independently of species). Applied **after** the positional sort, because the
rule is about order on the read, not about which segment scored higher; a refused pair collapses to
the single higher-scoring D rather than reaching down the score list for a producible partner.

On the TRB amplicon this takes tandem calls **15 → 5**, deleting **only** the impossible ones —
TRBD2→TRBD1 7→0, TRBD2→TRBD2 3→0, TRBD1→TRBD2 **5→5** — with the single-D count identical either
way (18,362). The shipped `examples/dd.airr.tsv` record (TRD, `TRDD2*01 → TRDD3*01`) is genomic
order and survives untouched.

⛔ **This does not make TRB tandem D-D real.** Under a flank-only shuffle (100 permutations,
conditioned on a real first D, D1's span fixed): 5 observed against 2.71 expected, Poisson
`p = .139`. The pre-fix `p = .0031` was not evidence either — its "excess" *was* the 10 impossible
calls, which such a shuffle cannot generate and therefore under-counts. What the gate buys is that
the residual signal is composed only of producible pairs.

⛔ **IGH is deliberately absent from the table.** In *human* IMGT the second number of
`IGHD<family>-<position>` is the genomic position; in *mouse* it is a family-member index with no
locus meaning, and the two vocabularies collide on real gene names (`IGHD1-1`, `IGHD2-15`,
`IGHD5-5`, `IGHD5-12`, `IGHD6-6` exist in both). `_map_d` is handed sequences, not an organism.
IGH tandem D-D is already 11 → 2 on real bulk IGH from the `/OR` orphon exclusion.

### Added — the clonotype table can now be cut up: D-D markup columns

`correct --map-d` named a `d2_call` and gave no way to find it. It now also emits `np1`/`np2`/`np3`
and `v_sequence_end`, `d_sequence_start`/`d_sequence_end`, `d2_sequence_start`/`d2_sequence_end`,
`j_sequence_start` — 1-based closed, in **junction** space (the clonotype table has no read), with
`-1` for "not located". `DCall.markup(junction_nt)` returns the same cut as labelled parts.

The partition closes: `np1 + D1 + np2 + D2 + np3 == junction[v_sequence_end : j_sequence_start-1]`
on **every** record carrying a `d2_call` — 5/5 read-level and 4/4 clonotype-level on the TRB
amplicon, 0 broken. ⛔ The boundaries *inside* the junction are one consistent reading, not ground
truth: chew-back and N/P addition make the V-end / np / D / J-start partition non-identifiable from
sequence, hardest for D, which is trimmed at both ends. The tests assert on **calls** and on the
partition **closing**, never on an NDN-internal boundary.

⚠ Real IGH tandem D-D is **0 in every local dataset at every stage** (0 of 15,070 mapped IGH reads,
0 of 4,783 assembled contigs, 0 of 1,742 IGH clonotypes) — 100 bp reads do not span a D-D. IGH
coverage of this path is synthetic.

### Added — quality-aware error correction: `map --junction-quality`, `correct --min-junction-q`, `--ec-mode`

`correct` used no quality information at all, so a sequencing miscall and a real low-frequency
variant were the same object to it — both are "a rare neighbour of an abundant clonotype" — and
`--error-rate` could only trade them off globally. Phred does separate them, because it is a
different measurement. Measured at the mismatching base over 310,559 real MIGEC spike-in windows:
the parent clone's own bases sit at median Q 38 (5.1 % below Q30), the two published spike-in
variants at Q 35 and Q 34 (17.6 % / 16.7 %), and the 1-substitution error cloud around them at
**median Q 24, 54.3 % below Q30**.

Three pieces, all **off by default**; the shipped output does not move.

* `arda rnaseq map --junction-quality` adds a `junction_quality` column — the read's Phred+33
  string over exactly the bases of `junction`, same orientation. Stage 1 is the only place the
  FASTQ quality is still in hand (it was read solely for `merge_pair`'s tie-break and discarded).
  +2.2 % wall and +4.4 % bytes on 100 k amplicon reads. ⛔ For a `rev_comp` hit the quality belongs
  to the read as submitted while every coordinate is on the coding strand, so it is reversed and
  then **verified against the junction it claims to describe** — a same-length slice off the wrong
  strand is a corruption nothing downstream can detect. Verified on 4,370 junction-bearing reads
  (4,365 of them reverse-strand) against an independent FASTQ extraction: 4,370 exact, 0 wrong.
  Refused with `--reconstruct`, whose merged fragment has no single input quality string.
* `arda rnaseq correct --min-junction-q Q` drops a read whose junction differs from its putative
  parent at any base below Q. Matching bases are not evidence and are never read, so this is not a
  min/mean quality over the junction. A clonotype with no more-abundant neighbour is never gated, a
  read with no quality string is kept, and a missing `junction_quality` column **raises** rather
  than silently not gating.
* `arda rnaseq correct --ec-mode fast|accurate` presets it. `fast` is byte-identical to today's
  default; `accurate` is `--min-junction-q 20`, the low end of the plateau (the effect is flat over
  Q20–32 and eats real variants by Q35). `binom`/`betabinom` are deliberately **not** in a mode:
  measured on the same 302 k-read library at `--error-rate 1e-5`, `betabinom` returns byte-for-byte
  the same 1,633 clonotypes as `simple` at 98× the wall, and `binom` reaches 23 clonotypes but at
  57× — and on a monoclonal precision arm both are exact no-ops (301 clonotypes either way).

What it buys: keeping both published MIGEC variants used to cost Jurkat TRB purity .99540 → .96034.
With the gate on it costs nothing — 2/2 variants at TRB purity **.99530 (Q20) to .99600 (Q35)**,
at or above the shipped default's purity, which keeps neither variant. Jurkat's spurious load falls
from 297 to 62 distinct junctions with the true clone untouched. ⛔ It cannot rescue the 0.0072 %
variant: that one sits below the RT template-error floor, whose competitors are high-Q by
construction.

### Added — `v_mutations` / `j_mutations`: the SHM record, in germline coordinates

Two new AIRR columns, `G45A,C112T` — germline base, 1-based position **in that segment's own
germline allele**, read base. Built in the walk `segment_cigars` already makes, so the marginal cost
is **36 ms per 100,000 mapped reads** (measured on 35,825 real bulk IG alignments: 26.9 → 39.7 ms,
cigars byte-identical) and **+2.25 %** on the output TSV (21,461,770 → 21,944,646 bytes).

The information was not missing. `sequence_alignment` / `germline_alignment` already carry every
column, and the germline they report matches the shipped allele on **28,365 of 28,365** mapped reads
of a real bulk IG library (66,526 V mismatches, zero disagreements). What was missing is that
recovering it needs arda's scaffold geometry — a consumer that does the obvious thing and diffs the
two alignment strings gets 100,091 mismatches on that library, of which **20,140 (20.1 %) are N-pad
or constant-region columns**: it attributes junction positions to a germline.

⛔ Which is why the scoping is structural rather than a filter. A mutation inside the V..J interior
is not attributable to any germline — recombination chews the segment ends back and adds
non-templated N/P bases, so the V-end / NDN / J-start partition of a junction often is not
identifiable from the sequence at all. The lists are built only for the V and J segments; the pad is
not a segment, so an NDN position has no germline coordinate to be filed under and cannot enter the
list by any code path. Substitutions only: indels stay in the CIGAR (`I`/`D`), and germline
coordinates past an indel are still correct because the walk tracks the target position across the
gap columns.

An AIRR/SAM extended CIGAR (`=`/`X` in place of `M`) was the alternative and is 34 % smaller
(+272,524 bytes against 411,202 on the same run), but it does not carry the germline base — so a
consumer still has to fetch the allele and re-index it — and it changes the meaning of a shipped
column instead of adding one.

### Fixed — conflicting CDR3-anchor rows were resolved by file order

IMGT ships two accessions under one allele name (mouse `IGKV10-96*01` is both AF441451/287 nt and
M15520/286 nt; `IGLV2*01` is J00599 and M17529), so `cdr3_anchors.tsv` can carry two rows for one
`(segment, allele)` with **different `germline_nt` and `status`**. `load_anchors` was last-wins, so
which junction germline the Cys104 gate scored against was decided by row order — on **3 mouse
alleles**, invisibly. `IGLV3*01` could resolve to a `truncated` row over an `ok` one.

Resolution is now explicit and logged: prefer `status == "ok"`, then the longer templated germline.

⛔ The conflict test compares only the fields that **decide the junction** (`anchor_nt`,
`germline_nt`, `templated_aa`, `status`). A TRAV/DV allele legitimately appears twice — once from
the TRA pass, once from TRD's `v_shared` — differing only in `locus`; treating those as conflicts
would emit 15 warnings per human load and train the reader to ignore the 3 that matter.

### Fixed — the flagship speed path was unreachable on a plain `pip install`

`segments.fasta` is **generated, not shipped**, and the auto-fetched reference tarball does not
carry it. So on a plain `pip install arda-mapper`, every `--two-pass` / `--fast-segments` run
degraded to the one-pass search behind a single log line: correct output, exit 0, and none of the
speed — the configuration that makes arda faster than MiXCR, silently unavailable out of the box.

`_cached_segment_db` now generates it (~0.3 s, once, under the same build lock) instead of
returning `None`. Verified: with `segments.fasta` deleted, a `--two-pass --fast-segments` run
rebuilds it **byte-identically** and produces **byte-identical** output.

⛔ Generation and stale-format *re*generation are separate functions on purpose, because they ask
different questions and need different done-predicates. `_has_jc_targets` is false for a missing
file, so reusing the stale-format predicate would make a missing file read as *already regenerated*
and the lock would skip the build — silently, in the same direction as the bug. `_has_jc_targets`
also used to raise `FileNotFoundError` on a missing file, which would have thrown inside the lock;
it is now total.

### Fixed — `arda rnaseq run` could reach only the dominated config

The pipeline entry point — the one a Nextflow or SLURM user actually calls — exposed `--two-pass`
and nothing else. And `--two-pass` **alone** is a loss: 0.762× on bulk and 0.87× on an IGH
amplicon. The single tunable it offered was the one configuration that makes arda slower.

`--fast-segments`, `--v-only-on-segment`, `--prefilter` and `--indel-rescue` are now on `run` too,
with the regime rule in the command's own help: `--two-pass --fast-segments --v-only-on-segment`
is the amplicon configuration, `--prefilter` is the bulk one, and they do not compose.

### Added — `arda export-ref`: the reference, with its markup, out of the CLI

The reference is arda's most valuable offline artifact — every in-frame V·J germline scaffold with
IgBLAST-quality FR1–4 / CDR1–3 coordinates, plus the per-segment markup and the per-allele CDR3
anchors — and until now it was only reachable by reading the build's TSVs by hand and re-joining
them against the FASTAs. That join is exactly the kind of thing that goes wrong quietly:
coordinates are **1-based closed**, a `J + C` scaffold has no V at all, and the aa reference has
three frames per D allele, so an off-by-one produces plausible nonsense.

```sh
arda export-ref --locus TRB                          # 2,112 TRB scaffolds, TSV, regions as columns
arda export-ref --kind segments --format fasta       # the 924-target segment reference
arda export-ref --kind anchors  --locus TRB          # per-allele CDR3 anchors
arda export-ref --format gff3 -o trb.gff3            # regions as features for a browser
arda export-ref --format airr                        # scaffolds shaped as AIRR Rearrangements
```

Three kinds (`scaffolds`, `segments`, `anchors`) × four formats (`tsv`, `fasta`, `gff3`, `airr`),
with `--locus` filtering and `--seqtype nt|aa`. GFF3 is 1-based closed like arda, so coordinates
pass through unchanged; the TSV states the convention in a header comment.

The tests do not check that the exporter ran — they check that what it wrote **round-trips**: every
region's `*_seq` equals the slice its own `*_start`/`*_end` imply, the declared junction equals the
slice its CDR3/FR4 coordinates imply, and regions are ordered and non-overlapping.

## 2.11.1

### Changed — `transfer_hit` no longer scans the junction anchor twice

`_anchored_vj_bounds` computes the longest prefix of the junction that the called V's germline
templates, and `v_anchor_prefix` then recomputed **the same scan over the same slice** to feed the
Cys104 junction gate (`_junction_nt` cuts `query_seq[cs-3-1 : f4+2]`, byte-identical to the window
`_anchored_vj_bounds` already built). It now returns that number and the gate reuses it, so a read
with an ambiguous V call no longer walks every called allele a second time.

Amplicon 100 k reads, 8 threads: **5.55 s -> 5.35 s**. Output byte-identical on an amplicon and a
bulk library; 624 tests pass.

Cumulative over 2.11.0 + 2.11.1, same workload: **5.96 s -> 5.35 s (10.2 %)**, against MiXCR's
5.90 s — 1.10x faster on wall at 12.7 s of CPU versus 45.2 s, and 631 MB versus 3,027 MB.

## 2.11.0

### Fixed — a junction was emitted even when the V was trimmed past its own Cys104

A rearrangement can trim V back beyond the conserved Cys104. The scaffold projection still lands
somewhere, so arda emitted a junction opening on bases the V germline never templated. Measured
against IgBLAST on 100,000 TRA amplicon reads (arda-benchmark `results/round18`): **1,396 of 46,785
junctions disagreed, and every single one was a pure 5′ boundary offset** — the 3′ [FW]118 end was
correct in 100 % of them, 1,376 contained the true junction and 20 were contained by it, none was a
different region. 648 were exactly +9 nt and 643 of those were one V gene, `TRAV25`, whose Cys104
sits at germline 264 while those reads' V ended at 258.

arda already computed the signal and threw it away: `_anchored_vj_bounds` measures the longest
prefix of the junction window that the called V's own `germline_nt` explains, and returns
`(0, 0)` — *without* suppressing the junction — when that prefix is 0. `transfer_hit` now refuses
the junction (and `junction_aa`, `cdr3_aa`, `productive`, `vj_in_frame`) when the prefix is under
`MIN_V_ANCHOR_PREFIX`. Precision among emitted junctions **.96953 → .99919**, for 92 of 44,414
correct junctions; on a bulk IG library the cost is **zero** — every emitted junction already
cleared the bar. The cut is 2 nt, not 3: a synonymous TGT→TGC is the one SHM event that reaches the
conserved codon, and it leaves two bases.

### Fixed — `j_call` was copied from the scaffold whether or not the read carried a J

`transfer_hit` set `j_call = ref.j_call` unconditionally, so a V-only read inherited the J of
whichever V×J scaffold it landed on. On a bulk RNA-seq library that was **1,823 of 2,737 mapped
reads**, 1,776 of which had an empty `j_sequence_start` and an empty `junction` in arda's own
output — the record already said there was no J alignment. `j_call` is now blanked when neither the
scaffold alignment reached J germline nor the junction anchor found a J start. `j_call` precision
against IgBLAST: bulk **.1129 → .7842**, amplicon **.9685 → .9953**.

Gated on the scaffold declaring `j_sequence_start`/`vj_end`: the aa markup does not, and blanking
on a reference that cannot say deleted every protein-input `j_call`.

### Changed — `transfer_hit` walks the alignment ONCE instead of four times

Besides the seven regions, `transfer_hit` needs three single scaffold positions projected onto the
query: the V germline end, the J germline start, and the V coding-frame anchor. Each went through
`_project_point`, i.e. its **own** `_markup.transfer_regions` call -- a fresh 6-argument binding
crossing, two fresh `std::string` copies of the *same* alignment, and a fresh forward walk. Measured
at ~443 ns each against ~822 ns for the real multi-region call.

Projecting a point is the degenerate region `[p, p]`, so all three now ride along in the single
multi-region call and are read back by index. The coding-frame anchor is knowable before the walk
(it depends only on `ref.starts[0]` and `hit["tstart"]`).

**`transfer_regions` calls per 100 k-read amplicon run: 207,007 -> 54,178.** Output is byte-identical
on both an amplicon and a bulk library.

`_project_point` is kept -- only positions known *before* the walk can be folded in.

### Changed — AIRR TSV formatting moved into the `_markup` C++ extension

`airr_out.format_rows` did, per record, 52 `dict.get` calls, 52 `None` tests, 52 exact-type tests,
a 52-element list build and a `str.join` -- 2.8 M dict lookups and 2.8 M list appends per 54 k-record
chunk. It was the largest single block of Python left on the per-read path.

`_markup.format_rows` does the same work with the column-name objects hashed once and the chunk
accumulated into one buffer. **In-run A/B on 100 k amplicon reads: 0.365 s -> 0.171 s (2.1x), about
3.4 % of end-to-end wall.** Output is byte-identical, asserted by `tests/unit/test_airr_out.py`
across filled records, all-`None` records, missing keys, non-string values and non-ASCII.

⚠ The microbenchmark says 3.4x and cProfile attributes 0.914 s of own time to this function; both
overstate it. The uninstrumented in-run measurement (0.365 s, 6.3 % of wall) is the one to quote.

The Python version stays as `_format_rows_py` -- the reference implementation, the fallback when the
extension is not built, and the thing the equivalence test compares against.

### Added — `build-db --allow-chimeras`: the TRDV x TRAJ scaffolds the default refuses

The default reference declines to build `TRDV x TRAJ` on the grounds that TRDV1/2/3 are dedicated
delta V genes, so the pairing is a chimera (`refbuild/loci.py`, pinned by
`test_tra_does_NOT_share_the_trdv_stem`). That is a biological claim, and the external evidence
disputes it: on 48,030 TRA amplicon reads, IgBLAST calls `TRDV1` + a `TRAJ` on **530** of them
(1.10 % of the library, median v_score 93.8, every one carrying a junction) and **MiXCR
independently agrees**, while arda calls the same J as both tools and emits no `v_call` at all.
That single class is **83 % of arda's entire remaining `v_gene` gap** on that library.

Which side is right is a domain judgement, so it is now a flag rather than a silent default.
`--allow-chimeras` gives TRA `v_shared=("TRDV", "")`; everything else is untouched.

⛔ **Measured, and the flag does not deliver the whole class.** Of 22 TRDV alleles, 15 are
`TRAV/DV` genes already present under TRAV, so exactly **7 dedicated TRDV alleles** are new
(V 102 -> 109). They imply 483 scaffolds, of which **476 are dropped for incomplete IgBLAST region
markup** and 7 survive -- all `TRDV1*01`, against TRAJ13/16/24/39. Human scaffolds go
15,414 -> 15,421.

End-to-end on the same 100 k amplicon reads:

| | default | `--allow-chimeras` |
|---|---|---|
| v_gene | .9867 | **.9874** |
| v_allele_exact | .9455 | **.9479** |
| TRDV1-class reads recovered | 0 / 530 | **34 / 530 (6.4 %)** |
| of which junction byte-exact vs IgBLAST | — | **34 / 34** |
| whole library v_gene newly right / newly wrong | — | **+37 / -3** |

Every junction it recovers is correct, and the ceiling is the markup pipeline, not the flag: the
truth's TRDV1 reads use TRAJ52 (x165), TRAJ8 (x21), TRAJ54 (x19) and other Js that have no
surviving scaffold. Off by default. Assert the **scaffold count** after building, never the flag.

### Added — `--v-only-on-segment`: a J-less read is aligned against its V segment, not 15,414 scaffolds

Off by default; requires `--two-pass`. A `v_only` read carries no J — that is the class, not a
search failure — so searching it against the whole V×J reference asks a question its sequence
cannot answer. It is 77 % of the amplicon rescue set, at **338 µs/read against 31 µs** for a
named-target alignment. MMseqs2 still does the alignment and still returns a real bit score, over
exactly the nucleotides a whole-scaffold alignment of a J-less read would have covered, so
`--min-score` keeps its meaning; anything that fails falls through to the full-reference rescue, so
no read is lost.

⛔ The class is gated by **geometry**, not by the shortlist reason. `v_only` means "no J segment
hit", which on a 100 nt bulk read carrying SHM is not "no J in the read": the segment pass misses
short hypermutated IGHJ and the full reference then finds it. Only reads whose V alignment stops
before their own Cys104 are routed. The separation is total — of the reads whose rescue *did*
produce a junction (77 bulk, 8 amplicon), **every one** reaches Cys104. Ungated this cost 77 of 213
bulk junctions; gated it costs **zero reads and zero junctions** in either regime.

Measured, 100,000 reads, 8 threads, M3, against `--two-pass --fast-segments` alone:

| library | wall | user CPU | reads lost | junctions lost |
|---|---|---|---|---|
| TRA amplicon | 6.81 s → **5.96 s** | 20.4 s → **13.3 s** | 0 | 0 |
| bulk RNA-seq | 3.19 s → **2.60 s** | 11.1 s → **5.4 s** | 0 | 0 |

## 2.10.0

### Fixed — a target-inverted alignment row silently produced a phantom clonotype

arda detects a reverse-strand nt hit only from the **query** side (`rev = qstart > qend`,
`annotate/mapper.py`). When MMseqs2 expresses the minus strand on the **target** side instead —
`tstart > tend`, with `taln` reverse-complemented — the row reads as forward. `transfer_regions`
then walks the target strictly forward from `tstart` (`_markup/markup.cpp:185-201`), which slides
the whole scaffold markup by `(tlen + 1 - tstart) - tstart` nt and takes the junction window off
Cys104 onto whatever codon lands there.

Nothing downstream catches it. The only gate on a re-annotated contig is `assemble._CANON`'s
`^C…[FW]$` — a conserved-motif test, not an anchor — so a window that happens to open on a
spurious `TGT` and close on `TGG` passes, in frame, with no stop codon.

Measured on a delivered **Jurkat** run (`ERR3003543`, arda 2.5.6 — Jurkat is a monoclonal T-cell
line and must yield essentially one TRB and one TRA):

| | junction_aa | reads |
|---|---|---|
| true clone | `CASSFSTCSANYGYTF` | 15,380 |
| **phantom** | `CVLLCQQFLDLFGSLW` | **7,408** |

`tlen` 349, true `tstart` 170, reported `tstart` 180 → the window slid **exactly 10 nt** into V
framework 3. Ablating the two phantom contig rows and re-running `correct` moves the true clone
**15,380 → 21,138**: the phantom had taken 5,758 reads from it. The same shape appears in the mouse
sample (`CNVFLCSMVQHPF`, 1,239 reads, shifted +13 nt from `CALWYSTHYVF`) — two organisms, one bug.

**These rows are not recoverable minus-strand hits to be reflected into forward coordinates.** They
are internally inconsistent: on that read `germline_alignment` is the reverse complement of the
scaffold while the query matches the plus strand at 91.4 % identity — the row's own reported
`pident` — and `identity(qaln, taln)` is 0.232, which is why `v_identity` came out 0.216 against
~0.98 for every normal record. They are now dropped in `_best_hits`, the one place that already
knows a row is unusable, so the read routes to the full-reference rescue or stays unmapped. That is
the only outcome that cannot ship a well-formed junction that is wrong.

⚠ Emission is **MMseqs2-build dependent**: 120 such rows across the six delivered samples, **0** on
the build in the local env at the same arda version. So this is a robustness gate rather than a
regression fix, and comparing two arda versions on one machine cannot catch it. Blast radius when
they do occur: 13 of 108,069 Jurkat reads (0.012 %), 3 of which emitted a junction.

Regression tests: `tests/unit/test_best_hits.py::test_a_target_inverted_row_is_not_a_usable_hit`
and `::test_an_inverted_row_never_beats_a_forward_one`. No reference DB, no MMseqs2 — the row is
expressible directly as a TSV line, so the test cannot be skipped away.

### Changed — `--max-subs` default 2 → 3

`max_subs` is a seqtree **search radius**, not a threshold: the accept/reject decision is the
length-scaled probability model (`count[parent] * (error_rate * L)**n >= count[child]`). A radius of
2 truncated the search below what the model would already accept on a deep clone, leaving a
sequencing-error trail uncollapsed in exactly the samples where the right answer is known.

Measured on four delivered libraries, clonotypes out:

| sample | kind | `max_subs=2` | `max_subs=3` | `max_subs=4` |
|---|---|---|---|---|
| Jurkat | monoclonal T line | 74 | **57** | 57 |
| Raji | monoclonal B line | 91 | **58** | 58 |
| GM12878 | oligoclonal B-LCL | 13 | **13** | 13 |
| MouseSpleen_WT | polyclonal spleen | 7,942 | **7,942** | 7,942 |

The polyclonal and oligoclonal repertoires are **unchanged** at 2, 3 and 4 — the model refuses those
collapses on abundance regardless of radius — so this is not a diversity-destroying change. It
saturates at 3.

### Added — `annotate.project`: the junction placed by arithmetic, not by alignment

AIRR `junction` runs from the Cys104 codon to the [FW]118 that opens FR4. Today arda learns those
two positions by aligning the read against a `V + pad + J` scaffold, and the **15,414-scaffold
reference exists for that and essentially nothing else** — the V and J *gene* calls already come out
of the 924-target segment pass at .9997 / .9998 agreement.

But the segment pass already returns `(target, tstart, qstart)` per read per side, and
`cdr3_anchors.tsv` already records `anchor_nt`. `refbuild.segments` writes each segment target so
its germline starts at offset 0, so target coordinates *are* germline coordinates:

```
offset = (anchor_nt + 1) - tstart
pos    = qstart + offset                    # forward hit
pos    = (qlen - qstart + 1) + offset       # reverse-complement hit
```

Two integer adds. Measured against IgBLAST truths at `v_score >= 70` on **254,867 scored reads**:

| sample | locus | n | accuracy |
|---|---|---|---|
| IGH_naive | IGH | 78,394 | **.99977** |
| IGH_repertoire (91.77 % median V identity) | IGH | 90,663 | **.99634** |
| SRR5233635 | TRA | 41,881 | **.99947** |
| SRR5233641 | TRB | 43,520 | **.99949** |
| IGH_repertoire | IGK | 143 | **1.00000** |
| IGH_repertoire | IGL | 17 | **1.00000** |

Accuracy is >= .993 in **every** V-identity stratum, including `<90 %` at n = 38,007.

**Not wired into the output path.** It yields `junction` and its coordinates; it does *not* yield
`v_identity`, `sequence_alignment`, `germline_alignment`, the per-segment CIGARs or the `mmseqs2_*`
block, all of which `annotate.transfer` derives from the alignment's `qaln`/`taln`, and
`_align_implied` also decides the allele. Whether removing junction placement from the critical path
recovers wall time is a separate measurement.

**It refuses rather than degrades**, on the eight conditions in `project.REFUSALS`: missing or
non-`ok` anchor (`no_anchor`), an unvalidated locus (`unvalidated_locus`), an indel gate that was
never run (`indel_unchecked`), `segmap`'s two-diagonal indel signature (`indel_split`), V and J on
opposite strands (`strand_mismatch`), an order violation (`order`), a projection landing off the
read (`off_read`), and a junction length that is not a multiple of 3 (`bad_codon`). A well-formed
junction that is wrong is the worst output this codebase can produce — the reference-geometry bug
shipped junctions that started `C`, ended `[FW]`, passed `--complete-only` and were short by exactly
the allele's truncation.

⛔ **TRD is declined, because it has ZERO coverage.** The per-locus bar was ">= .99 at n >= 2,000, or
the locus goes on the refusal list". Across two TR amplicons the segment pass never handed a single
TRD read both anchors, so all 767 TRD junctions in the truths fell through to the aligner and TRD
never appears at all. *Absent* is not *validated*, and the only TRD number that exists is 43/51 =
0.843. `UNVALIDATED_LOCI` declines it, keyed off the **J** anchor — TRAV/DV alleles rearrange to
either TRAJ or TRDJ and the J decides the locus.

⚠ Scope is set by yield, not accuracy: **87.0 % of hit TRA-amplicon reads and 77.3 % of TRB carry
both anchors, against 7.2 % of bulk reads** (bulk is 53.0 % `V_only`, 31.4 % `C_only`). Bulk reads
mostly do not span a junction at all, so this is an amplicon lever.

### Added — two off-by-default measurement hooks

`ARDA_MMSEQS_SEARCH_OPTS` appends flags to every `mmseqs.search` argv, and `ARDA_PROJECT_DUMP` /
`ARDA_VONLY_DUMP` write per-read diagnostic tables. None is a supported interface; they exist
because every recorded flag measurement in the benchmark was taken against a reference that has
since been rebuilt twice, so re-measuring had to become cheaper than editing source. The first
immediately earned its keep: `--exact-kmer-matching 1` turns out to flip allele-level calls on
IGH_naive at an *identical* clonotype count, which no count-based metric can see.

### Changed — `project_junction`'s `split_checked` is required, with no default

An inert indel gate must be loud. `split_checked` tells `project_junction` whether the caller has
already run `segmap`'s two-diagonal split check; a default would let a caller skip the check and
still get a projection, which is the failure mode the `indel_split` refusal exists to prevent.
Making it a keyword with no default turns that into a `TypeError` at the call site instead. The new
`indel_unchecked` refusal covers the case where the caller passes `False`.

## 2.9.0

### Added — `--indel-rescue`: an ungapped extension cannot score a read that has an indel

An ungapped extension follows ONE diagonal, so a read carrying an insertion or deletion relative to
germline scores only up to the indel. The unit test measures it directly: a 120 nt read scores
**240** against its germline target, and **120** — exactly half — with a single base deleted. Its
scaffold then gets chosen on truncated evidence.

That is not a corner case. Measured on **341,294 real IGH mates** (two IGH RepSeq amplicons,
IgBLAST truth, `v_sequence_alignment`/`v_germline_alignment` gaps — no cigar, nothing inferred):

| V identity | `>=98` | `95-98` | `90-95` | `<90` |
|---|---|---|---|---|
| reads carrying a V indel | 0.74 % | 1.63 % | 3.56 % | **8.00 %** |

3.18 % pooled, and the rate tracks SHM load because AID produces indels, not only substitutions.

The signature is in the seed votes before any extension runs: two well-supported diagonals on the
**same** target, offset by the indel length. Votes are already sorted by `(target, diagonal)`, so
one pass reads it.

**A flagged read is REROUTED to the gapped rescue, never dropped** — so a false positive costs a
little speed and cannot cost a read. The demotion happens after `shortlist()` has asserted its
partition is total, so the no-read-is-lost invariant is untouched. Counted as `indel_rescued` in
the run report.

**What it changes, measured.** Recall is *identical* with and without it (174,066 of 174,226 real
amplicon fragments either way), because recall asks whether a read was found and these were found
regardless. What moves is which allele and which junction they get — which is load-bearing, since
the clonotype key is `(locus, v_call, j_call, junction)` at allele level. Adjudicated against
IgBLAST at gene level on exactly the reads whose call moved:

| sample | v_call moved | correct without | correct with |
|---|---|---|---|
| IGH_repertoire (hypermutated) | 586 | 53.24 % | **84.13 %** |
| IGH_naive (low SHM) | 100 | **77.00 %** | 63.00 % |
| **pooled** | **686** | 56.71 % | **81.05 %** |

**+167 reads corrected, +24.3 points.** The sign flip is the mechanism confirming itself: where
real indels are common the flag fixes calls truncated at the indel; where they are rare its false
positives — repeats read as two diagonals — dominate. `locus` never moves.

⛔ **Off by default, and it is `--fast-segments`-only.** On 13 bulk RNA-seq datasets it demotes
**zero** reads and the AIRR output is byte-identical, which is correct: bulk TR carries ~0 indels.
Turn it on for hypermutated IG work; it does nothing elsewhere.

### Measured — `--fast-segments` on real IGH amplicon

**1.87x at 41 % less memory** (100,000 pairs, ~90 % receptor: 319.74 s → 170.77 s, RSS 4,016 →
2,382 MB), because it raises `fast_fraction` from **0.052 to 0.5018** on the same reads.

⛔ **Every `fast path` figure in the `--two-pass` documentation is MMseqs2-specific.**
`fast_fraction` is a property of *(reads x segment mapper)*, not of the reads: MMseqs2 misses the
short IGHJ on 95 % of these 5'RACE reads and `_segmap`'s ungapped extension finds it on half.
`v_only` rescues fall 169,004 → 85,933. Read the old table as-is and you leave the fast path off on
exactly the library where it is worth 1.87x.

Recall cost, against an independent IgBLAST truth on 344,554 real mates: **.99525** vs the
default's .99582 — 197 reads, of which **200 are below 90 % V identity and 0 are above 90 %**.

## 2.8.0

### Added — `--fast-segments`: the segment pass answered structurally, not by homology search

The two-pass path searched a 924-target segment reference purely to learn, per read, its best V
allele and its best J allele with coordinates — no cigar, no backtrace. That is a **structural**
question about a fixed 236 kb germline reference, not a homology search. `arda._segmap` (new C++
extension) answers it by seeding, voting by diagonal and extending ungapped. Measured on 100,000
amplicon reads against the shipped reference, 8 threads:

| | wall | agreement with `mmseqs search` |
|---|---|---|
| `mmseqs search` | 2,770 ms | — |
| `_segmap` | **74 ms** | V allele .9997, J allele .9998, C allele 1.0000 |

End to end on 50,000 pairs, `--two-pass` with and without the flag: **9.10 s → 6.12 s (1.49×)**,
`locus` identical on every read, `v_call` .999794, `junction_aa` .999938, 6 reads lost of 48,620.

⛔ **Off by default, and the residual delta is why.** Six reads and ten V calls of ~48,600 is small
but is not zero, and the shipped path does not move them at all. It only *nominates*: every
candidate is still aligned against the full V+pad+J scaffold and scored by MMseqs2.

What made it equivalent were two constants, both wrong at first and both caught by measurement:

* **k = 12**, the `-k` arda already passes MMseqs2 — not 16, which was inherited from `prefilter`
  where it is calibrated for *rejection*. Seed length sets sensitivity to mismatches: at k=16 the
  mapper seeded 53,048 reads against mmseqs' 53,121, and a read with no segment hit is assumed
  non-receptor and is **never rescued**, so those 73 were simply lost. At k=12 `no_segment_hit`
  matches mmseqs exactly and the AIRR delta fell from 30 lost reads to 6.
* **a significance floor of 40.** mmseqs applies `-e 1e-3`; this scheme has no e-value, so without
  a floor 43,010 reads pick up a constant-region hit against mmseqs' 473. Half of those score
  exactly 38 — a bare seed plus a couple of flanking matches.

### Fixed — IgBLAST was run without its J-frame table, so every truth had NO junctions

`optional_file/<organism>_gl.aux` tells `igblastn` each J allele's reading frame; without it there
is nothing to place the Phe/Trp 118 anchor against, so it emits no `cdr3`, no `junction` and no
`junction_aa` **on any read** — while still calling V and J, reporting a normal `v_score` and
exiting 0.

Both callers looked for the file under `paths.bin_dir()`; it lives beside the executables under
`igblast.igblast_root()`. Those are **the same directory in a source checkout** (`setup.sh`
installs IgBLAST into `<repo>/bin`) and different on every auto-fetched install, so it worked
everywhere it was developed. And `auxiliary_data=aux if aux.exists() else None` made the miss
silent.

Measured cost: a 10,000-read amplicon truth with `j_call` on 9,070 of 9,300 reads and `junction_aa`
on **zero**, written up as an IgBLAST limitation at 151 bp before it was traced here. With the file
passed, 437 of 449 reads carry one, 4–19 aa. `igblast.auxiliary_data()` now **raises** rather than
degrading: "no junctions" is indistinguishable from a truth that genuinely has none.

⚠ Any IgBLAST truth built with an auto-fetched install before 2.8.0 has this defect. Check
`junction_aa` fill before scoring against it.

### Fixed — a misleading FASTQ diagnostic

`_read_pairs_dnaio` reported every dnaio `FileFormatError` as "R1 and R2 differ in length; one file
is truncated". dnaio raises that type for malformed *records* too — a real case here is a `+` line
that kept its original SRA description after the `@` line was renamed to carry a mate suffix. It
now only claims a truncation when dnaio actually reported a pairing or length problem.

### Reference — the TRA/TRD V-gene sharing is ASYMMETRIC, and that is deliberate

`TRD` declares `v_shared=("TRAV", "/DV")`; `TRA` declares nothing, and a symmetric-looking
"fix" for that was built, measured and **reverted**. The sharing is not symmetric because the
biology is not:

* **TRDV1/2/3 are dedicated δ V genes** and rearrange to **TRDJ**. `TRDV1 + TRAJ` is not a
  rearrangement that occurs, so building that scaffold invites reads onto a chimera.
* **TRAV/DV genes pair with either, and which J they took is what defines the locus.** Both real
  directions are already covered: `TRAV/DV × TRAJ` comes free with the TRA build (IMGT files those
  genes under TRAV) and `TRAV/DV × TRDJ` is what `TRD`'s `v_shared` adds.

What made the wrong version look justified was an IgBLAST truth calling 147 amplicon reads
`TRDV1*01` + a TRAJ with a real junction, against which arda scored 0.0952 on that stratum. arda
declining to call a V there is **arda being right about the biology**; the truth is wrong.

The shipped reference is therefore **unchanged** — a full `build-db` + `build-index` from IMGT
reproduces the committed `markup.tsv`, `alleles.fasta` and `combinations.tsv` byte for byte
(15,414 human scaffolds) — and `_REFERENCE_TAG` stays at 2.5.7. `test_tra_does_NOT_share_the_trdv_stem`
pins the asymmetry so the symmetric version cannot be reintroduced by inspection.

### Added — tests for FASTA input and `--limit`

Both features already worked and neither was tested. 17 cases, covering the paths that actually
differ: `read_pairs` sends the *unlimited* case to dnaio and the *limited* case to a pure-Python
reader. Includes format detection by content rather than extension, paired FASTA agreeing with
paired FASTQ record for record, a truncated FASTA mate still raising, `--limit` counting **pairs**
rather than reads, and a limit not failing on a divergence beyond it.

### Changed — the segment reference collapses the J×C product too

`refbuild/segments.py` exists to remove a cross-product: it replaced 15,414 V×J scaffolds with 775
V + 124 J targets. It did that on the V side and **copied the 345 J+C scaffolds through verbatim**,
which are themselves a J×C product (IGH 14 J × 11 C, IGL 9 × 7, TRB 16 × 2). Every J+C scaffold of
a locus ends in the *same* constant sequence, so a read reaching C was aligned against all of them
to learn one `c_call`, at a redundancy factor equal to the locus' J-allele count — **69× on TRA**.

Measured on a TRA amplicon: those 345 targets were **27.7 % of the database and 76.4 % of its
alignments**, 4,977 alignments per target against 603 for a V target.

The constant region is now its own target, exactly as V and J are — one `C|<allele>` target per
distinct C allele, **345 → 25**, whole reference **1,244 → 924**. A read spanning J into C hits its
`J|` and its `C|` target separately and names its J+C scaffold through
`Reference.jc_combinations()`, the C-side twin of `combinations.tsv`.

**This is faster and more accurate**, the same way the V×J collapse was — a J call decided by a
whole-scaffold bit score whose *constant* half is arbitrary is a worse J call. On 50 k TRA amplicon
pairs, scored against the one-pass output on the same reads:

| | before (1,244 targets) | after (924 targets) |
|---|---|---|
| segment search | 6.00 s | **3.18 s** (1.89×) |
| whole `map --two-pass` | 11.87 s | **8.92 s** (1.33×) |
| segment alignments | 2,248,847 | **531,703** (4.23×) |
| `v_call` disagreements | 147 | **85** |
| `j_call` disagreements | 401 | **296** |
| `c_call` disagreements | 13 | 13 |
| `junction_aa` disagreements | 18 | 18 |
| reads lost | 2 of 48,627 | 7 of 48,627 |

The five extra losses are V-less J→C reads carrying **no junction**, so no clonotype is affected: a
read spanning the J/C boundary can fall below the search threshold on each half where the
concatenated scaffold cleared it. They are counted as `no_segment_hit`, not silently dropped.

⚠ **A C target is kept for every locus, including the six where the C call carries nothing.** Only
IGH's constant genes separate anything reportable — its 11 alleles are 7 classes, i.e. the isotype;
TRA/TRD/IGK have one C allele each, TRB/TRG two, and IGL's seven IGLC are all one class. Dropping
those 14 targets is measurably faster (2.89 s vs 3.18 s) and costs no fast-path read. They are kept
because a C target does a second job unrelated to information content: it is the only segment
target a read lying wholly inside the constant region can hit, and without it such a read hits
nothing, never enters `seen`, and is never rescued — **14 of 453 reads vanish** on the real-read
fixture, every one a V-less J→C read.

⛔ The J+C contest is nominated **from the J, not from a C hit**. Requiring C evidence is not
equivalent, and it let the exact bug the contest exists to prevent back in: on
`SRR5233639.12648/1` it invented `TRBV12-3*02`, destroyed `c_call` TRBC2*01, and fabricated
`junction_aa` CASSFAGLVNIDEQFF on a read the one-pass calls V-less.

⛔ `_segment_rows`' polars filter and `_segment_best_hits`' own kind guard are two statements of
one rule in two languages. Adding `C|` to the loop alone made the reduction discard every C row, so
`best_c` was always empty and 15 J→C reads vanished with **`no_segment_hit` not even moving** — the
rows were dropped before anything counted them.

### Fixed — six silent two-pass defects, found by audit

All six produce correct-looking output and none raises. Each was reproduced against the shipped
human reference before being touched.

* **`_segment_rows` had no unusable-row filter.** polars sorts nulls FIRST under
  `descending=True`, so a row with an empty `bits` field became rank 0 for its `(query, side)` and
  **evicted the read's real best hit** — reproduced with `V|A*01 <empty bits>` beside `V|B*01 120`,
  where the reduction returns `V|A*01` and the 120-bit row is gone. The null then reached
  `float(row["bits"])` and raised mid-chunk, after earlier chunks were written: a partial AIRR file
  that looks complete. `_best_hits` has had this filter since 2.7.2; the segment path never got it.
* **Composite allele names could never match `combinations.tsv`.** A segment target inherits its
  scaffold's `v_call`/`j_call` verbatim and those are sometimes ambiguity lists, while
  `load_combinations` registers only individual members. Measured: 23 of 775 `V|` and 2 of 124 `J|`
  targets are composite, and **all 2,852** (composite V × any J) pairs were absent — zero hits.
  `IGHV3-23*01,IGHV3-23D*01` and `IGKV1-39*01,IGKV1D-39*01` are on that list, i.e. the most-used
  human IGHV and IGKV genes: every such read fell to rescue, reported as a chimera the reference
  does not contain.
* **`jc_combinations` keyed on the comma-joined group string**, and the two sides group alleles by
  *different* rules — a J+C scaffold collapses identical J sequence, a `J|` target inherits the V×J
  collapse of identical assembled `V+pad+J` plus reading frame. **24 J+C scaffolds** were
  unreachable from any J hit, including every IGLJ2/IGLJ3 read. Now keyed per allele on both sides.
* **A pre-2.8.0 `segments.markup.tsv` poisoned that lookup.** It loads into the same `entries` dict
  *after* `markup.tsv`, and its `JC|` rows collide with their base scaffolds on `(j_call, c_call)`;
  the later insertion wins, so every value became an id the full target DB does not contain and the
  contest was off for the whole run. Reachable on the first run after any upgrade.
* **The J+C contest accepted a walkover.** When the V×J alignment produced nothing, its row was
  taken unconditionally and the read dropped from the rescue set — keeping a V-less answer (no
  `v_call`, no junction, no clonotype) the full reference was never asked about. The read now stays
  in `failed`.
* **`build_segment_reference` truncated in place on the concurrent runtime path.** Now under the
  same `arda._locking.build_lock` the mmseqs DB build uses, with the format re-checked after
  acquiring it, and a read-only reference tree degrades to a warning instead of killing the run.

### Fixed — an exact J tie was broken lexicographically

Found by the test suite, not the audit. J alleles of one gene are short and differ by a base or
two, so a read routinely ties EXACTLY between two `J|` targets — and the tie was resolved by target
name, which silently decided which V×J scaffold the read was aligned against.
`SRR5233639.3589/2` ties `J|IGLJ2*01,IGLJ3*01` and `J|IGLJ2A*01` at **54 bits each**; the comma
sorts before `A`, so the composite won and the read was seated on a scaffold scoring **93** while
its true home scores **96**. `_MAX_TIED_J` mirrors `_MAX_TIED_V`: tied J alleles now contribute
candidates and `_best_hits` decides on whole-scaffold bit score, as the one-pass does. Candidates
are the two axes (tied V × best J, best V × tied J), bounded by the sum of the caps, not the
product.

Measured on 50 k TRA amplicon pairs against the one-pass output: `j_call` disagreements
**296 → 214**, `v_call` 85 → 90, `c_call`/`junction_aa` unchanged.

**Blast radius: `--two-pass` only.** `segments.fasta` is reached solely under `if two_pass:` in
`rnaseq.map`, so the shipped default path — the one-pass search against the full V×J reference — is
byte-for-byte unaffected by everything above. `--two-pass` is off by default and is documented as an
amplicon lever.

`SegmentStats.jc_targets` is now `SegmentStats.c_targets`. `segments.fasta` is generated by
`build-index`, not shipped, so no reference asset changes.

## 2.7.2

### Fixed — a second `--two-pass` crash, on alignment rows with no score

2.7.1 dropped rows with a blank query id; the amplicon run then died at the next site instead:
``TypeError: float() argument must be a string or a real number, not 'NoneType'``. ``bits`` is cast
with ``strict=False``, so an unparseable score becomes null rather than raising at parse time, and
that null reached ``float(row["bits"])`` in the J+C contest. Same malformed row, second crash site.
Both are now filtered in one place and both route the affected read to the full-reference rescue —
the guarantee ``--two-pass`` is built on.

With the two fixes, ``--two-pass`` completes on amplicon for the first time. Measured on 1 M pairs
(PRJNA371303, `map` stage): **455.33 s vs 804.69 s (1.77×)** on TRA and **479.35 s vs 825.83 s
(1.73×)** on TRB, losing **10 and 16 reads** of ~1 M, with `junction_aa` concordance **.9995** and
**.9984** against the one-pass output. It is the only amplicon lever measured so far that preserves
calls — `--adaptive` is 3.15× but moves 34 % of V-gene calls in that regime, and a
one-allele-per-gene reference costs 42 % of allele-level V calls on TR data by construction.

⚠ These are defensive fixes, not a root cause. 16–30 % of reads still fall through to the
full-reference rescue on amplicon, which is most of what ``--two-pass`` gives back in speed, and
the malformed rows suggest the same origin. The sub-DB construction still needs an audit.

### Fixed — a stale claim about `createsubdb`, and the code written on top of it

``_subset_db``'s docstring said ``createsubdb`` does not write a ``.lookup`` for the destination.
It does (verified against mmseqs 18-8cc5c, which also copies it WHOLE rather than subsetting it).
Writing the file by hand on that assumption overwrote mmseqs' own and broke the rescue path.
New ``tests/unit/test_subset_db.py`` pins the real behaviour — four invariants, each of which has
shipped broken at least once and every one of which fails **silently**.

## 2.7.1

### Fixed — `--two-pass` crashed on amplicon after writing a partial output

`arda rnaseq map --two-pass` died with `KeyError: None` roughly 90 % of the way through a 1 M-pair
TCR amplicon run, having already written ~872 k of ~969 k rows. A truncated AIRR TSV plus a
non-zero exit is the worst shape of failure here: the output file looks like a finished one, and a
pipeline that only checks for the file's existence sees a completed run.

Cause: `_align_implied` searches a hand-built prefilter sub-DB, and `mmseqs createsubdb` does not
carry header/`.lookup` entries for everything in it, so `convertalis` emits rows with a **blank
query field**. `_best_hits` turned that into a `None` dict key and `seqs[None]` raised.

Such a row cannot be attributed to any read, so it is now dropped and logged. That is the correct
repair rather than a workaround: a read absent from the hit mapping is realigned against the full
reference by the caller, which is exactly the exactness guarantee the two-pass path is built on.

Found by running `--two-pass` on amplicon for the first time since the J->C contest fix, the
tied-V fix and `top_hit`.

### Performance — the record builder and the chunker

* `segment_cigars` and `_aln_identity`, the two per-alignment-column Python loops in
  `transfer_hit`, move into the C++ extension: **23.03 -> 0.69 us** (33x) and `transfer_hit`
  **130 -> 27 us** per mapped read (4.8x). Both keep their Python originals as the reference the
  ports are differentially tested against, and as the fallback when the extension is not built.
* `format_rows` no longer looks up each column twice (1.22x on that term).
* `chunked_fragments` computes `frag_stem` only at a possible cut point instead of once per
  record -- 1.2 M calls became ~10 on a 1.2 M-read run. Cut decisions verified identical over 300
  randomised trials including the `--reconstruct` parity case.

Context for the sizes: on a receptor-rich library the record builder is ~6 % of a run and
`mmseqs search` is ~78 %, so these are real but bounded. On amplicon the search is **96.9 %**.

## 2.7.0

### Added — `--prefilter`: an exact k-mer gate in front of MMseqs2

On bulk RNA-seq the receptor fraction is 0.02–3 %, so `mmseqs search` spends essentially all of
its time proving that reads are *not* receptor reads. Measured on 4 M reads of SRR10611239 (947
map, 0.024 %): the search alone is 48.9 s. The fitted cost model says why —
`wall ≈ reads/46,353 + hits/350`, so the dominant term is the **read count**, not the answer.

A read can only align to a V(D)J scaffold if it shares an exact 16-mer with one. Testing that is a
lookup; proving it is Smith-Waterman. `src/_prefilter/` (pybind11, sibling of `_markup`) does the
cheap test first: exact 16-mers, both strands, ≥ 1 hit, two-level lookup (a 4^12 bitset — 2 MB and
L2-resident — then a binary search into a sorted `uint64` array), `std::thread` fan-out, early
exit on the first hit. Against the real human reference: **157,838 seeds / 1.3 MB**, built in
0.21 s, scanning at **36.4 M reads/s** on 8 threads.

`arda rnaseq map --prefilter`, **off by default**. Measured on aldan3 (`map` wall from the run
report, 32 threads), against TRUST4 on the same reads:

| sample | receptor % | before | after | speedup | TRUST4 | ratio |
|---|---|---|---|---|---|---|
| SRR10611239 | 0.024 | 82.58 s | 7.79 s | **10.6×** | 23.46 s | **0.332** |
| SRR6926533 | 0.123 | 39.11 s | 5.26 s | **7.4×** | 8.71 s | 0.603 |
| SRR8363894 | 0.772 | 46.86 s | 13.19 s | **3.6×** | 11.48 s | 1.149 |

**Why not an MMseqs2 flag.** Its own prefilter cannot be tuned into this and cannot act early
enough. A 15-setting sweep found **no lossless candidate above 1.05×**; the informative row is
`--min-ungapped-score 30`, which is free (0 reads lost) and gives *no speedup* — so the cost is
the k-mer stage, which no flag exposes. And MMseqs2 can only prefilter reads already in a DB, so
the FASTA write and `createdb` (12.0 s, 19.6 % of the 4 M-read profile) are paid for every read
regardless. This runs before both.

**What it costs.** ~0.5 % of real reads, and the shape matters more than the number:

* every lost read scored **75–79 bits** against the `--min-score 75` cutoff — max lost score 79,
  so nothing confident is ever dropped;
* **`junction_aa` moved on zero shared reads** — the clonotype key is
  `(locus, v_call, j_call, junction)`, so a moved junction splits a clonotype, and none moved;
* the loss is **entirely IG, zero across all four TR loci** — the predicted signature, since one
  substitution destroys k consecutive exact windows and SHM supplies substitutions.

The index is built from `Reference.target_fasta`, the **same FASTA MMseqs2 searches**. That is
load-bearing: against a `V+pad+J`-only reference the design measured 16.29 % loss, **69.27 % of it
J→C reads**, and indexing the constant region took it to 0.53 %. Deriving the index from the
search target makes that hole structurally unreachable. `prefilter_stats` in the run report gives
`seen`/`passed`, the only number that says whether it earned its keep on a given library.

### Changed — search batches are no longer read chunks

Every `mmseqs search` call costs ~0.7 s of fixed setup whatever it is handed (10 chunks = 42.73 s
of search against 36.45 s for one), and a prefiltered 400 k-read chunk holds ~1,900 reads. Reading
stays chunked so memory is bounded; survivors now accumulate until they amount to a chunk's worth
of work. Raising `--chunk-size` to 4 M reaches the same wall at **741–1642 MB** RSS — batching
gets it at a quarter of that. Flushing only ever happens on a chunk boundary, never inside one:
`chunked_fragments` keeps a fragment's mates together, and splitting them is a bug this pipeline
has shipped once already under `--reconstruct`.

### Changed — FASTQ parsing moved to C (new dependency: `dnaio`)

Once the prefilter removed the search, **reading became the largest single cost** of a bulk run —
65 % of a 0.024 %-receptor library, where before it was 3 % and explicitly not worth touching.
Reworking the pure-Python loop (`zip(*[iter(fh)] * 4)` rather than four `readline()` calls) bought
1.43×; dnaio is 2× on top of that including mate tagging.

It is adopted because it makes the **same two assertions** the pure-Python reader makes, and both
are load-bearing — a truncated mate file and a shuffled one each produced a published false
discovery in arda's benchmark that had to be retracted. Its errors are re-worded so callers still
see "mate mismatch" and "truncated" as distinct failures, and a truncated gzip still surfaces as
`ValueError` rather than the bare `EOFError` Typer prints as "Aborted.". It is **not** used when
`--limit` is set: limit is a head, so a truncation beyond it must never be reached, and dnaio
validates pairing while filling its own buffers. The pure-Python reader remains as a fallback, and
a test asserts the two produce identical streams.

### Added — `--adaptive` on `rnaseq map`

`_extend_uncertain` has implemented the align-term lever since 2.6.0 (`--max-accept 40`, then
re-search uncapped only the reads scoring under 90 bits; 2.17× on 1 M bulk reads, zero reads lost)
but it was reachable only from the Python API. Off by default: on the real-read fixture it moves
`junction_aa` on 3 of 453 reads, two of them at 128 and 131 bits, so the trigger cannot be
calibrated on score alone. **Measured worse than plain `--prefilter` in the prefiltered regime**
(11.13 / 16.01 / 8.19 s against 10.20 / 14.66 / 6.42) — once the scan term is gone the re-search
is pure overhead.

### Measured and rejected

`--two-pass` on top of `--prefilter`. The theory was that prefiltered survivors are amplicon-like
(34,582 of 41,145 hit — 84 %) and so land in its 3.51× regime rather than the 0.762× bulk one.
They do not: 14.69 / 22.84 / 9.17 s against 10.20 / 14.66 / 6.42. Its cost is not the hit rate.

## 2.6.3

### Fixed — the mmseqs auto-fetch could publish a half-written binary

arda runs concurrently against one cache **by design** — a Nextflow process per sample, a SLURM
array task per shard — and on a cold cache every one of them calls the fetch in the same moment.
That path installed with `shutil.copy2` straight onto `bin/mmseqs`, which is this repo's
signature bug in its purest form. Two silent failures:

* `copy2` writes **in place**, so two processes doing it at once interleave their bytes and
  produce a corrupt binary that nothing reported an error about;
* `dest.exists()` — the gate deciding the work is already done — goes true on the **first byte**
  `copy2` writes. A third process therefore finds an mmseqs, trusts it, and executes a truncated
  file. The `chmod` came *after* the copy as well, leaving a second window in which the file
  existed but was not executable.

Now there is one downloader under `build_lock`, the binary is made executable while still under
a private name, and a single `os.replace` publishes it. A rename is atomic within a filesystem,
so `bin/mmseqs` only ever exists as a complete, executable binary — which is what makes
`dest.exists()` a legitimate gate instead of a bug. (`os.replace` is a rename only within one
filesystem, which is why 2.6.2's move of the staging directory next to the destination is a
prerequisite, not a separate fix.)

The two new tests were checked against the old implementation and **both fail there** — a
concurrency test that passes on the broken code pins nothing. The decisive one interrupts the
copy *after* it has written bytes, which is what a full disk, an OOM kill or a SLURM timeout
does, and asserts no partial binary is left for the next process to run.

The reference fetch (`_database_fetch`) and the new IgBLAST fetch already worked this way; this
brings the third one into line.

### Fixed — `pip install 'arda-mapper[mmseqs]'` was a hard failure

The `mmseqs` extra depended on `arda-mmseqs`, a companion distribution that is built by CI but
has never been uploaded to PyPI. An extra naming a package the index does not have does not
degrade gracefully — it is a resolution error:

```
ERROR: No matching distribution found for arda-mmseqs<19,>=18; extra == "mmseqs"
```

And that exact command was what the README, the installation docs, and `mmseqs_binary()`'s own
"binary not found" error all told you to run. The one actionable line in the error message was
itself broken.

`mmseqs` is now an empty alias, like `rnaseq` already was, so existing pins keep resolving.
Nothing is lost: the bundled wheel was only ever a convenience over the runtime auto-fetch,
which pulls the same static binary, is what a plain install already does, and needs the network
exactly once — as `pip` itself just did. The four places that recommended the extra now describe
auto-fetch, and the resolution test asserts the advice names a route that exists.

There is no install-time alternative to look for: **a wheel has no post-install hook.** `pip`
unpacks an archive and runs nothing, by design, so "download it during `pip install`" is not
something a package can opt into. Fetching at first use is the standard answer to that
constraint, and it is what arda already does.

## 2.6.2

### Fixed — auto-fetch staged its download in the system temp dir, and died on a cluster

Both auto-fetches unpacked into `tempfile.TemporaryDirectory()`, i.e. `/tmp`. On a cluster `/tmp`
is routinely a small node-local disk: aldan3's is **2.0 GB with 29 MB free**, and the IgBLAST
fetch — a ~400 MB archive that extracts to more — died there on the first real run with
`OSError: [Errno 28] No space left on device`.

Both now stage on the **destination filesystem** (`dir=dest.parent`), which also makes the
download → extract → lay out → `os.replace` sequence single-filesystem by construction rather
than by assumption. MMseqs2 is fixed the same way; it had the identical bug and would have hit
it next.

Invisible on any developer machine with a roomy `/tmp`, so a test asserts *where the archive is
written*, not merely that the install succeeded.

## 2.6.1

### Fixed — `arda igblast` did not work from a `pip install`, at all

IgBLAST resolved only through `<project>/bin`, which is populated by `setup.sh`. A plain
`pip install arda-mapper` has no source checkout and never runs `setup.sh`, so **every**
`arda igblast` invocation failed. It is the command that produces the gold standard the
benchmark scores against, so on a fresh install the accuracy arm of a 27-dataset panel run
returned nothing on all 27 — the tool runs themselves were fine.

The error also pointed at the wrong thing:

```
IgBlastError: IgBLAST ships no internal annotation for organism 'human'
```

That names a missing data file inside a present install. The actual state was that IgBLAST had
never been installed, and `has_internal_annotation` was resolving `internal_data/human/` under a
directory that did not exist. Reading it as a broken or unsupported reference costs an hour.

IgBLAST is now fetched from NCBI on first use, exactly as MMseqs2 already was:

```sh
pip install arda-mapper
arda igblast -i reads.fq -o truth.tsv     # fetches the release once, then runs
```

Resolution order is `$ARDA_IGBLAST` → `<project>/bin` if a checkout already has one →
`<cache>/igblast`, auto-fetched. `$ARDA_IGBLAST_ASSET` overrides the platform asset and
`$ARDA_NO_AUTO_FETCH` refuses the download with an error that names the fix. Windows resolves
an asset now too (`x64-win64`), where the old script raised `SystemExit`.

The install is **atomic and concurrency-safe**, because arda runs concurrently against one
cache by design. The release is laid out in a *sibling* staging directory and moved with a
single `os.replace`, and the readiness gate is a marker written last rather than
`igblastn.exists()` — this repo has shipped that exact bug twice, and here it would surface as
executables present with `internal_data/` still copying, i.e. as the misleading error above.

`scripts/fetch_igblast.py` is now a wrapper over the packaged module instead of a second copy
of the logic; keeping two is what let the runtime path go missing. `igblast_version()` reports
the installed NCBI release (1.22.0 at time of writing) so a benchmark can record it.

## 2.6.0

### Fixed — a run was not reproducible

`arda rnaseq correct` did not produce the same output twice. Three runs over one unchanged
200k-read AIRR, same flags, gave three different `clones.tsv`. polars' `group_by` is a
multithreaded hash aggregation and the clonotype fold leaned on it three times without pinning
an order: group order (so `_parents` collapsed an error child onto whichever parent it met
first — one identical junction flipped between `IGHV3-11*06` and `IGHV3-21*08`), read order
within a group (so `_assign_coverage`, which is first-with-longest-overlap-wins, moved
`duplicate_count` between 11 and 9 for one IGK clonotype), and the order of equal-`count` rows.

None of it was visible. The row *count* was stable at 449 every time, so the table looked
reproducible while its V calls and abundances were not.

Two more coin flips in the same class, upstream: `_best_hits` resolved exactly-tied bit scores
with an unordered sort and an unordered `unique` — paralogous scaffolds tie routinely, and
**all 25** tied queries in a fixture flipped their V call when the input rows were reversed —
and chunk boundaries could fall between a fragment's two mates, costing that fragment the
isotype donation `_apply_constant_rule` makes per chunk. The default paired path got away with
the latter (records arrive strictly two per fragment) but `--reconstruct` breaks that parity.

After the fix, on 200k reads: all three artifacts byte-identical across 1, 4 and 8 threads;
`correct` byte-identical over five repeat runs and six input-row permutations. The shipped
example is unchanged and the realworld IgBLAST-concordance suite still passes — this is a
reproducibility fix, not an accuracy change.

### Added — `arda rnaseq slurm`, byte-identical to a single-node run

`arda slurm` could not run the RNA-seq path, and aiming it there would have been silently
wrong: it writes FASTA (dropping quality), round-robins *records* so a fragment's mates land in
different shards, and has no Stage 2/3 at all.

New `arda rnaseq split | reduce | slurm`. Only Stage 1 is distributed — `correct` counts
distinct fragments globally and `assemble` grows contigs *across* reads, so sharding either
would count a clone once per shard and never build the long-CDR3 contigs Stage 3 exists for.
Shards are **contiguous blocks of read pairs**, so concatenating the per-shard AIRR in shard
order reproduces the single-node row order exactly, and `run` and `reduce` call the *same*
Stage-2/3 function rather than two copies that could drift.

Verified on 60k real read pairs, single-node vs 5 shards: `airr.tsv`, `assembled.airr.tsv` and
`clones.tsv` all byte-identical.

The amplicon `arda split` / `arda slurm` are unchanged, and now say in their help that they are
not for paired RNA-seq.

### Added — resource reporting for every stage, and provenance in the report

`assemble` and `correct` now report `wall_seconds`, `peak_rss_mb` and `rss_gain_mb`; `map` was
moved onto the same helper so all three mean the same thing. This matters because Stage 3 is
the expensive one: mapping is flat at ~300–400 MB at any depth, but the clone set scales with
repertoire richness (2.7 GB at 28k clonotypes vs 314 MB for a colder sample with *more* reads),
so a `--mem` directive sized from the mapping number alone gets OOM-killed.

`peak_rss_mb` is documented as what it actually is: the whole-process high-water mark **as of
that stage's end**, monotone, because `getrusage` offers no per-stage reset. That is the number
a memory directive has to cover; `rss_gain_mb` gives the per-stage attribution.

The merged report also gained `mmseqs_version` and a `reference` fingerprint — when two
delivery modes disagree, a different aligner build or reference is the usual cause and neither
is otherwise visible.

### Added — mmseqs that just works, and is the *right* mmseqs

`mmseqs_binary()` took whatever was on `PATH`. An index is only reusable by the release that
compiled it, so a mismatched binary made every run discard `database/`'s precompiled DBs and
rebuild a private cache — no error, just a slow start and results not comparable with anyone
else's. Seen on a cluster whose `~/bin/mmseqs` shadowed conda's version-matched one.

Candidates are now version-matched, with auto-fetch of a known-good build as the fallback and a
warning that names the consequence if even that fails. `$ARDA_MMSEQS` stays unchecked (an
explicit override is the user's call), and no filtering happens when no index ships — a plain
`pip install` has nothing to be compatible with.

`versions_compatible()` also handles a spelling `version_key` could not: the static release
asset prints its **full commit hash** where bioconda and the index marker print
release+short-commit. Release 18 *is* commit `8cc5c…`, so those are one build — and since the
static asset is what arda auto-fetches, the stricter comparison would have rejected arda's own
downloaded binary on every macOS install.

New `pip install 'arda-mapper[mmseqs]'` ships the binary in a companion wheel
(`packaging/arda-mmseqs`), so there is no first-run download and no network needed. A separate
distribution because pip cannot be asked to *prefer* a wheel — build tags are ranked, not
chosen; the same shape `cmake`, `ninja` and `ruff` use. MMseqs2 is MIT and its notice ships
with the wheel.

### Changed

* The Nextflow module pins `mmseqs2 =18.8cc5c` exactly (was `>=15`). A floating pin lets conda
  resolve a different aligner than the CLI and SLURM paths use — an accuracy-differs-between-
  modes hazard, not just a caching one.
* The module gained the `meta.yml` it never had.

### Added — a segment reference, and a fast path that cannot lose a read

Searching 15,414 V×J scaffolds re-pays for the shared V region once per J that V pairs with —
median 13, and 67 for TRA. `arda.refbuild.segments` emits every V, J and J+C allele as its own
target instead: **1,244 human targets, 12.4x fewer**. A read's best V and best J then name
exactly one scaffold through `combinations.tsv`, so the second alignment is one target rather
than ~277. Measured on a TRA amplicon: **2.20x end to end, 85 % of reads on the fast path.**

The reference is derived from `alleles.fasta` + `markup.tsv` in under a second, so `build-index`
generates it rather than shipping it — the same argument that already excludes the mmseqs
indexes from the 3.2 MB release asset.

`arda.annotate.shortlist` exists to keep this from becoming a different tool. arda's claim is
near-zero Stage-1 false negatives, so every read is partitioned into `implied` (fast path) or
`rescue` (realigned against the full reference), and the partition is asserted total. V-only,
J-only, an unknown V×J pair, a failed second alignment — all rescued, none dropped. Measured: of
5,278 hits on real amplicon reads, **one** was lost, at bit score 56 — below the default
`--min-score 75`, so nothing arda would have reported. Gene-level agreement with the one-pass
path: locus 99.96 %, V 99.09 %, J 99.24 %.

Two mistakes in this path produce *correct output that is silently no faster*, so both are
now documented and tested: `JC|` targets are keyed by scaffold id rather than allele (getting it
wrong collapsed the fast path from 85.3 % to 0.1 %), and `top_hit` on the segment pass destroys
the V+J pairing the whole scheme depends on.

Strand is handled per read, not per library. The segment pass runs `--strand 2` and mmseqs
reports reverse hits with `qstart > qend`, which makes a hand-built prefilter diagonal
meaningless; reverse-strand reads are aligned against their own reverse complement and the
coordinates flipped back. A fixed-orientation assumption would have been wrong on the first
library tested — 99.1 % of its reads were reverse-strand, and library strandedness is not a
constant across protocols.

There is deliberately **no alpha/delta ambiguity rule**. TRD *is* TRAV/DV + TRDJ — the J (and C)
decides the locus — and the reference already encodes it: of 1,050 scaffolds built from a TRAV/DV
segment, 1,005 are locus TRA and 45 are locus TRD. An earlier draft rescued these as ambiguous,
discarding real rearrangements the reference contains. `combinations.tsv` is the arbiter.

### Fixed — the two-pass fabricated junctions on J->C reads

A read that runs through its J into the constant region has two plausible homes: the V×J
scaffold its best (V, J) pair names, and the J+C scaffold the segment pass actually hit. Among
775 V alleles a 100 nt read always has *some* V above threshold, so the shortlist always found a
pair — and forcing that choice took a read scoring **141** on a J+C scaffold, re-seated it at
**99** on a V×J one, destroyed its `c_call` (the isotype) and **invented a `junction_aa`** out of
the spurious V. 4 of 453 mapped real reads, 3 of them gaining a fabricated junction.

Both scaffolds now compete on bit score, exactly as they do in the one-pass search. Two targets
per read is still ~138x fewer alignments than the full reference. On real reads `locus`,
`c_call`, `junction_aa` and `productive` are now byte-identical to the one-pass output.

Two further corrections came out of wiring it up. V alleles of one gene differ by a nucleotide or
two, so a short read routinely cannot separate them — **a measured mean of 3.45 tie exactly** on
segment bit score. Picking one of them by segment score alone diverged from the one-pass V *call*
on 38 % of TRA amplicon reads, and since the clonotype key is `(locus, v_call, j_call, junction)`
at ALLELE level, that silently splits and merges clonotypes. Every tied allele's scaffold is now a
candidate (capped at 8; they share a diagonal, so it is one extra prefilter line each, not an
extra search).

That made the alignment TSV 2.88x larger — and it was being written with the full 17-column
format, `cigar`/`qaln`/`taln` included, rebuilding a fraction of the 194 MB -> 877 MB peak-RSS
regression `top_hit` exists to prevent. Reducing to one row per query first fixes the memory and,
unexpectedly, the calls: both paths now break exact ties by mmseqs' own ordering rather than by
two different rules, so **allele-level `v_call` agreement went 0.9735 -> 0.9956** and `junction_aa`,
`locus` and `c_call` are byte-identical. Deterministic across 2, 4 and 8 threads.

`--two-pass` is now reachable: it is wired into `arda rnaseq map` and `arda rnaseq run` (and so
through the SLURM and Nextflow paths, which call the same function), with the segment DB and
`combinations.tsv` built and parsed once per run rather than per chunk, and the shortlist
accounting reported under `segment_search` in `arda.json`. A missing `segments.fasta` falls back
to the one-pass search with a warning rather than failing.

**Off by default, and it is an amplicon optimisation rather than a general one.** The scheme
needs a read to hit both a V and a J so the pair names one scaffold. Primer-anchored amplicon
reads do (85 %); bulk RNA-seq reads land anywhere in a transcript and do not (5 %). Measured:

| library | receptor % | fast path | speedup |
|---|---|---|---|
| TCR amplicon | 48.47 | 85.13 % | **3.51x** |
| bulk RNA-seq | 2.74 | 4.98 % | 0.762x — 31 % **slower** |

On bulk the segment search is overhead on top of a rescue that is nearly the whole set. Neither
regime loses a read arda would report (0 lost at or above `--min-score`, 0 gained, both). The
bulk case is a *scan-term* problem and the segment reference structurally cannot address it.

### Fixed — `productive="F"` claimed 75 % of mapped reads are non-productive

`productive` is a property of the V-J junction, and `phase` is only computed when a read reaches
both CDR3 and FWR4. A read with a V but no junction fell through to `"F"` — a bare V fragment
reported as a *confirmed* non-productive rearrangement. On the real bulk fixture that is **342 of
453 mapped reads (75 %)**, since most bulk reads lie wholly inside V and never reach CDR3.

The module already had the rule right one branch over (a V-less read leaves `productive` empty —
"not non-productive, it is unevaluable"). `productive` and `vj_in_frame` are now empty unless a
junction was observed; `stop_codon` deliberately is not, because a stop in the V-side regions is
directly observed either way. The test asserts the invariant both ways — `productive` is set iff
`junction_aa` is — so neither an over- nor an under-report can pass.

### Fixed — the isotype vote counted mates, not fragments

`_dominant_ccall` documents itself as the dominant resolved class over a clonotype's *fragments'*
constant mates, but iterated per-mate `sequence_id`s. A fragment with both mates assigned voted
twice while a fragment with one assigned mate voted once, so a 1-fragment minority could outvote
a 2-fragment majority. One vote per fragment now, order preserved so the tie-break stays
deterministic. Only affects samples with uneven mate assignment.

### Fixed — `--error-method binom|betabinom` could hang forever

`_root` walks parent pointers in an unbounded `while`, so a 2-cycle is not a wrong answer — it is
a run that never returns. `_parents` cannot make one (it requires `count[parent] * p_err >=
count[child]`, forcing counts to increase along the chain), but `_error_pileup` re-decided
parentage from per-position depth with no ordering condition, and that test is symmetric at low
coverage: `_binom_sf(1, 2, 0.001) = 0.001999 > alpha`, so each of a pair passes as the other's
error child. Reproduced on real Stage-1 output (16,157 rows → 1,190 clonotypes): mutual pairs at
indices 77↔143 and 857↔861 under both methods. Fixed by restoring the invariant — a parent must
be strictly more abundant than its child. Not on the default path (`error_method="simple"` is
acyclic by construction), and the regression test lives in the DB-free suite, because a stage
that never terminates must not be gated behind an mmseqs skip.

### Fixed — tied V candidates shared one alignment diagonal

A V×J scaffold is `V + pad + J` with the V left-aligned, so a read sits at the same offset only
in scaffolds whose V allele has the same length — and alleles of one gene do not. On the human
reference **11 V genes carry alleles of differing length, 8 differing by ≥ 35 nt and one by 72**,
and `mmseqs align` returns nothing once the diagonal is off by more than ~35 nt, so the shared
diagonal silently dropped the sibling scaffold. Each candidate is now shifted by its V-end
difference. A candidate whose geometry is unknown keeps the unshifted diagonal and is no worse
off than before.

### Changed — each stage materialises only the AIRR columns it reads

`_error_pileup`, `_assign_coverage` and `assemble_contigs` each built Python strings for all 83
columns of the Stage-1 AIRR, including `sequence_alignment`, `germline_alignment` and the seven
region sequences. They subscript 6, 8 and 11 of them. Measured on 181,200 real rows: **2.42 KB/row
for all columns vs 0.41 KB/row for the used ones — `correct` peak RSS 1099 → 714 MB**. At full
depth (~3.6 M mapped rows) that is ~7 GB. `clones.tsv` is byte-identical across the full
three-stage pipeline.

### Added — `--adaptive`, measured and off by default

`--max-accept` is unbounded by default, so arda aligns every hitting read against all ~300 of its
prefilter candidates and keeps one. Capping it is the largest single lever on the align term
(75 % of bulk search wall). The cap alone is lossy, so `--adaptive` caps everything and then
re-searches *uncapped* only the reads whose capped score falls below a trigger. Measured on 1 M
real bulk reads: **2.17× with zero reads lost**, the uncertain set being 0.5 % of the library.

Off by default because read preservation is not the whole guarantee: it also changes
`junction_aa` on 3 of 453 reads on the real-read fixture — two of them scoring 128 and 131
against a 90-bit trigger — and moves **~23 % of the clonotype table** (allele-level Jaccard
0.7706 on 1 M bulk reads). A high score does not certify that the best alignment was found, so a
score-only trigger may not be calibratable at all. Shipped rather than deleted because the
measurement bounds how much of a bulk run is alignment work that never reaches the output; the
test pins the junction-move count so it cannot silently grow.

### Fixed — segment targets could not be resolved through `Reference`

`Reference` is built from `markup.tsv`, which describes scaffolds, so `ref.get("V|IGHV3-7*02")`
returned `None` and `_annotate_chunk` dropped every segment hit as unmapped: searching the
segment reference produced **0 annotated reads against the scaffold reference's 278** on the same
input, even though the search found them. `segments.markup.tsv` shares the schema and the key
space and now loads through the same path. Optional by design — a reference built before 2.6.0
has no segment entries rather than failing.

After this the segment reference annotates the same 278 reads and finds all 126 V-only reads. All
775 V segments agree with their scaffolds on FR1–FR3 region coordinates (0 differing), because a
scaffold is `V+pad+J` with the V at position 1. The two paths still disagree on the V *gene* for
13 of 126 V-only reads, because the segment search ranks paralogues by V alone while the scaffold
search ranks by `V+pad+J`; running `arda igblast` on exactly those reads agrees with the segment
call **9 times to 3**. Segment-as-primary remains unreleased pending the junction and clonotype
gates.

### Fixed — the two-pass tests were testing the one-pass search

Generating the segment reference in `build-index` and gitignoring it left CI — which never runs
`build-index` — with no `segments.fasta`, so `map_rnaseq` fell back to the one-pass search as
designed and the comparison tests compared the one-pass output *to itself*, and passed. Third
instance of this repo's signature failure: silent success over nothing. A module fixture now
builds the segment reference, and `test_two_pass_is_actually_engaged` guards the guard — it
asserts the fast path resolved something and says in its failure message that every other
two-pass assertion is vacuous otherwise.

### Changed — a measured cost model, and a bigger chunk

Fitting 13 full-depth cluster runs (60,252 s) gives `wall_map ~= total_reads/44,470 +
mapped_reads/681`: **a read that hits costs ~65x one that does not**. That splits `map` into a
scan term and an align term whose ratio is set by the library's receptor fraction, and it is why
the segment reference helps most on amplicon and least on a 0.0003 %-receptor negative. The
RNA-seq chunk moved 200k -> 400k on the same evidence. `scripts/bench_cost_model.py` refits it.

Refit directly against `mmseqs search` rather than against whole-run wall (round 5): `wall_search
~= reads/46,353 + hits/350`, where *hits* is reads with at least one prefilter candidate. It fits
three independent points to 0.0 s and carries **no fixed per-call term** — the earlier
`mapped_reads/681` was a whole-run figure imported into a per-hit decomposition, which forced a
spurious intercept. Four consecutive 400k searches take 35.0/34.8/35.2/35.9 s, i.e. flat, so
chunking costs 1.021x and no more.

Unit tests now run against **real reads** (660 pairs from public SRR5233639, 88 KB) rather than
synthetic ones only.

### Fixed — every cold-cache install of a bumped version 404'd

`reference_url` was derived from `__version__`, so releasing 2.5.7 pointed first-run users at a
`v2.5.7` asset that does not exist — the reference tarball only changes when the reference does.
Pinned to `_REFERENCE_TAG` with a fallback. Found by the cold-cache experiment, not by a test.

### Changed — lint gates CI, and the docs Makefile matches it

`ruff check src/` reported 16 errors and CI never ran it. Fixed and gated. `docs/Makefile` built
without `-W` while CI built with it, so a local docs build could pass on something CI rejects.

## 2.5.7

### Fixed — a stale mmseqs cache was reused instead of rebuilt

`_cached_target_db` detects a stale target DB by mtime (older than `alleles.fasta`) and calls
`_createdb_atomic` to rebuild it — but that build guard gated "already built" on bare **existence**
(`done=db.exists`, from the 2.5.5 concurrency work). So it found the stale file, called the work done,
and skipped the rebuild: every run then searched the previous scaffolds, projecting the current
`markup.tsv` coords through an out-of-date alignment and sliding the junction off Cys104. Surfaced as
shifted mouse markup — a 352-nt TRB scaffold cached under a reference since rebuilt to 346 nt gave
`junction_aa='LCASSFHRDYNSPL'` instead of `CASSFHRDYNSPLYF`. "Done" now means a *current* db (exists
**and** at least as new as its source fasta); the same skip also affected `build-index` after a
`build-db`.

### Fixed — the committed precompiled indexes were never used

The precompiled mmseqs DBs shipped in `database/vdj/*/mmseqs/` are gated on a version marker, compared
with an exact `==`. But `mmseqs version` prints the same release+commit with different punctuation
across builds — the official static binary says `18-8cc5c`, the bioconda build `18.8cc5c`. When the
toolchain moved to the static binary, every committed marker (written under conda's mmseqs) stopped
matching, so the shipped indexes were silently ignored and every run rebuilt a private cache. The
comparison now uses a separator-insensitive `mmseqs.version_key` (folds `[\s._-]+`→`-`, lowercases),
which bridges the cosmetic difference while still rejecting a genuinely different version.

### Changed — dev/CI toolchain moved from conda to uv

mmseqs2's bioconda binary was the only thing tying the dev environment to conda; it also ships an
official static build, which arda already auto-fetches. `setup.sh` now builds a `uv` `.venv` (portable
under bash and zsh) and CI installs via `uv`. Conda stays only for the Nextflow/ISP integration, which
ships its own `environment.yml`. Does not affect `pip install arda-mapper`.

## 2.5.6

### Fixed — a concurrent fetch could publish a half-extracted reference

2.5.5 made the mmseqs **index** build safe under concurrency. The **reference fetch** that runs just
before it had the same disease, and it was the more dangerous of the two.

`paths.database_dir` decides "is the reference here?" by testing whether `<cache>/database/vdj/` is a
directory. Its mere existence is the gate. So anything that lets `vdj/` become visible before it is
complete does not produce a slow download or a crash — it makes every *other* arda process search an
**incomplete reference** and report success over nothing.

Two things let that happen:

* **No lock.** arda is routinely run concurrently against the same cache — the Nextflow module
  launches one process per sample, a SLURM array one per task — so on first use in a fresh
  environment all of them fetched into the same path at once.
* **A cross-filesystem move.** The old code extracted into `/tmp` and called `shutil.move`.
  `shutil.move` is a rename only *within* one filesystem; across a boundary it silently degrades to a
  recursive **copy** — and `/tmp` and `~/.cache` normally are two different filesystems. So `vdj/` was
  populated file by file, in full view of everyone else. It also `rmtree`'d the existing tree before
  writing the new one, deleting files a concurrent reader could be mid-search on.

The fetch now holds a build lock, extracts into a staging directory **inside the destination** (so the
swap is a rename, not a copy), and publishes with a single `os.replace`. `vdj/` is either absent or
complete — never in between. Under `--force` the old tree is *renamed* aside rather than deleted, so a
reader already holding its files keeps them.

The lock is now one implementation (`arda._locking.build_lock`) shared by the reference fetch and the
index build, rather than two copies of the same loop.

### Changed

* `seqtree>=0.4` (was `>=0.2`). Clonotype output is unchanged — `examples/rnaseq/clones.tsv`
  reproduces bit-for-bit.
* Docs: corrected the memory guidance. Peak RSS tracks **repertoire richness**, not read depth —
  mapping is flat (~300–400 MB at any depth), but Stage 3 holds the clone set, so a B-cell-rich tumour
  peaked at 2.7 GB (28k clonotypes) while a colder sample with *more* reads used 314 MB. The old
  "< 400 MB, independent of read depth" was measured before the Stage-3 assembler existed and would
  have OOM-killed anyone sizing a SLURM or Nextflow memory directive from it. Budget ~4 GB.

## 2.5.5

Two ways a correct arda install could still fail to work. Both hit the delivery path.

### Fixed — concurrent runs corrupted the mmseqs index

**arda is routinely run concurrently against the same reference** — the Nextflow module launches one
process *per sample*, a SLURM array one per task — and the index build was not safe against that.

`mmseqs createdb` creates its `db` file the instant it starts writing. Every *other* arda process
therefore saw `db.exists() == True`, skipped building, and searched a **half-written database**. The
failure was silent and total: **`0/200000 reads mapped, loci={}`**, no error, clean exit code. A whole
27-dataset benchmark run came back as zeros. (`build_index` was worse: it unlinked the files a
concurrent reader was mid-search on, and mmseqs died with `Cannot open index file db_h.index.1`.)

The build now holds a lock, writes into a private temp dir, and moves the finished files into place
with `db` **last** — readers test `db.exists()`, so it must not appear until its siblings are all
there. A killed builder leaves a temp dir, never a half-built DB that looks complete.

*Shipping a prebuilt index would not have fixed this*: the index is keyed to the mmseqs **version**,
and arda uses whatever `mmseqs` is on `PATH`. Any user with their own build (as on a cluster) rebuilds
locally and races anyway.

### Changed — `seqtree` is now a core dependency, not an optional extra

`pip install arda-mapper` is all the bulk RNA-seq pipeline needs.

As an extra it produced this package's worst failure mode: a plain install would map and assemble a
100 M-read sample for 45 minutes and *then* die on a bare `ModuleNotFoundError` — after all the
expensive work, before writing a single clonotype. It slipped past every smoke test because
`arda --version` succeeds without it, and it shipped to a collaborator's Nextflow module, whose own
comments called `-profile conda` "works out of the box" while it could never emit a clonotype table.
Patching the error message (2.5.1) treated the symptom; the dependency should simply not have been
optional.

`seqtree` ships the same 12 wheels on the same platforms as arda (py3.10–3.13, linux/mac/win) and
pulls no runtime dependencies of its own, so requiring it costs nothing.

* `[rnaseq]` survives as an **empty alias**, so existing `arda-mapper[rnaseq]` pins (the Nextflow
  module, collaborator instructions, older docs) keep resolving without a warning.
* **Removed every `pytest.importorskip("seqtree")` gate** (10 of them). With seqtree required, a skip
  there can only *hide* a failure — and a skipped test is exactly how two reference bugs shipped.
* The Dockerfile keeps its `import seqtree` build check.

## 2.5.4

**Completes the 2.5.3 reference fix.** 2.5.3 dropped germlines that are *truncated* into the
junction, but left in a second class that is just as wrong and turned out to be doing most of the
damage: germlines whose junction anchor arda **cannot find at all**.

### Fixed

* **Unanchorable V alleles kept building scaffolds — and were read magnets.** When arda cannot locate
  Cys104 in a germline (`status=no_anchor`), IgBLAST still annotates a CDR3 on its scaffold, at a
  position arda has no way to verify. It is demonstrably the wrong one: `TRAV23/DV6*04` yields a
  scaffold junction `CTTSGTYKYIF` (11 aa) while every other allele of that gene templates `CAAS` and
  yields 14 aa.

  On a real tumour this was not a rounding error. After 2.5.3 dropped `TRAV20*03`, its reads did not
  come out correct — they **moved to `TRAV23/DV6*04` carrying the identical 11 aa junctions**, and
  that one allele then held **55 % of all TRA reads**. The gene label changed; the corrupt junction
  did not.

  2.5.3's mitigation (keep the scaffold, blank its `junction_aa`, set `productive=F`) does not work:
  the runtime projects the scaffold's `cdr3_start`/`cdr3_end` onto the read and re-derives a junction
  from those coordinates. And the check it used — "junction does not start with C" — never fired
  here, because `CTTSGTYKYIF` *does* start with C. It is simply the **wrong C**.

  Such an allele is now **dropped from the scaffold set when its gene keeps at least one anchored
  allele** — the sibling is what proves it wrong, and the drop is then free: **0 genes lost, in any
  organism** (46 human V alleles, 44 mouse, 2 rhesus, 1 rabbit). A gene whose alleles are *all*
  unanchorable is left alone: nothing contradicts its annotation, and dropping those would delete
  **41 *functional* mouse IGHV/IGLV genes**.

* **New CI invariant**: a scaffold's junction must agree with its own allele's anchor — no scaffold
  may be built from an unanchorable allele that has an anchored sibling, and no junction may be
  shorter than the V alone templates. (Length, not byte-identity: IMGT's gapped and ungapped records
  for the same allele can differ — mouse `IGLV2*01` disagrees at 12 of 294 bases — and the anchors
  table reads one file while the scaffold reads the other.)

Human reference: 15,690 → **15,069** V-J scaffolds. Mapping is unchanged; clonotypes and V calls move.

## 2.5.3

**Reference fix — germlines that are truncated INTO the junction no longer build scaffolds.**
Clonotype output changes; mapping does not. Re-run rather than diff against 2.5.2 clonotypes.

### Fixed

* **3'-truncated V alleles produced junctions that were short but looked canonical.** Some IMGT V
  records stop before the canonical 3' end, so the `V + pad + J` scaffold built from them **lacks
  nucleotides a real rearrangement has** — and every read projected onto such a scaffold came out
  with a junction short by exactly that much. `TRAV20*03` templates 3 nt into the junction where its
  siblings template 13, so every clonotype built on it was **3 aa short**. On a real tumour it took
  1,039 reads (48 % of all TRA reads) across 33 clonotypes, and arda's TRA CDR3-length distribution
  peaked 3 aa below every other tool's.

  Nothing downstream could catch it: the junction still starts with Cys and ends with [FW], so it
  looks canonical and `correct --complete-only` passes it. The build's completeness gate only checks
  that region *coordinates* exist, which they do.

  An allele is now dropped from the **scaffold set** if its germline reaches at least one codon less
  far into the junction than the longest allele of its own gene. Orphan-safe by construction — the
  threshold is relative to the gene's own best allele, so that allele always survives and no gene can
  lose all of its alleles. 63 human V alleles, 53 mouse (`TRAV20*03/*04`, `TRBV4-3*02/*03/*04`,
  `IGKV3-20*02`, …). They remain in `cdr3_anchors.tsv`, now flagged **`status=truncated`** — a status
  `docs/reference_build.rst` has documented from the start and which was never implemented.

* **Scaffolds whose V has no findable Cys104 claimed to be productive.** 570 human and 179 mouse V-J
  scaffolds emitted a `junction_aa` that does not begin at Cys104 — i.e. not a junction at all — and
  549 of the human ones were flagged `productive=T`. Such a scaffold now **keeps its place in the
  reference** (the read still maps, so filter recall is untouched) but carries **no junction** and
  `productive=F`, so `correct --complete-only` drops it instead of minting a bogus clonotype.

  Dropping those alleles outright was considered and rejected: it would delete **41 *functional*
  mouse IGHV/IGLV genes**, and arda is the only tool in the benchmark that runs on mouse.

* Effect on the shipped example: mapping is **bit-identical** (925/1035 reads, same per-locus
  counts); clonotypes go 21 → 20 because reads on a dropped allele reassign to a surviving sibling of
  the same gene, merging rows that shared a junction and differed only in `v_call`.

### Added

* **`arda build-db --one-allele-per-gene`** — build scaffolds from one representative allele per gene
  (`*01` where it exists, else the lowest-numbered). Human: **16,035 → 4,443 scaffolds (3.6x smaller)**
  with **all 290 V genes retained**. Deliberately not a literal `*01` filter: 19 human V genes
  (`IGHV2-70D`, `IGHV3-25`, `IGHV3-43D`, `IGHV3-62`, `IGHV3-64D`, …) have no `*01` record and would
  vanish silently. Off by default.

* **Reference invariants now run in CI**, without a built DB (`tests/unit/test_reference_invariants.py`):
  every scaffold junction starts at Cys104; no scaffold is built from a truncated allele; no gene loses
  every allele; TRD still carries the shared `TRAV*/DV*` genes; human still has its 345 `J+C` scaffolds.
  Both reference bugs found so far were invisible to CI because the only tests covering them needed a
  built database and were skipped — that is exactly how the 2.5.2 defect regressed.

## 2.5.2

**Annotation fix (γδ / δ chains).** TRA and TRD are interleaved on chr14 and share V genes filed under
`TRAV` as `.../DV...` (e.g. `TRAV14/DV4`). Whether those landed in the TRD scaffold set silently
depended on IMGT functionality flags in the `TRDV` file — a clean `arda build-db` could produce TRD
scaffolds with **no** `/DV` genes, so a δ rearrangement on such a V gene had no TRD scaffold to match
and was **miscalled TRA** (the locus follows the J/D/C gene, never the shared V).

* `refbuild`: `Locus.v_shared=("TRAV", "/DV")` now pulls the shared genes into TRD deterministically.
* Added a DB-free regression test — the existing `test_locus_disambiguation` needs a built DB, so it
  was skipped in CI, which is how this regressed.
* Human reference rebuilt; TRD scaffolds now carry `TRAV/DV`. No API change.

On real δ amplicon data each library now maps ~99.6–99.9 % to TRD.

## 2.5.1

Packaging and integration fixes. No change to mapping, assembly, correction or any output
column — 2.5.0 results reproduce exactly.

### Fixed

* **`pip install arda-mapper` could not run `arda rnaseq`.** `arda rnaseq correct` imports
  `seqtree`, which ships only in the optional `rnaseq` extra, so a plain install mapped and
  assembled and then died with a bare `ModuleNotFoundError` — after the expensive stages, and
  before writing any clonotype table. The extra is now documented in the README and in
  `docs/installation.rst`, and `correct` raises a `ModuleNotFoundError` that names the fix
  (`pip install 'arda-mapper[rnaseq]'`). A genuinely broken `seqtree` build still surfaces its
  own `ImportError` rather than being misreported as "not installed".

* **The Nextflow module never installed that extra**, so `-profile conda` — the path its own
  comments advertise as "works out of the box" — could not produce a clonotype table.
  `integrations/nextflow/arda/environment.yml` now pins `arda-mapper[rnaseq]`. The `Dockerfile`
  smoke test gained `python -c "import seqtree"`; `arda --version` succeeds without it, so the
  image built green and failed mid-pipeline.

* **The module's `airr` output channel emitted two files.** `path("*.airr.tsv")` also matches
  `<prefix>.assembled.airr.tsv`, which `rnaseq run` has written since 2.4.0, so `ARDA.out.airr`
  carried a 2-element list and the assembled AIRR was silently mislabelled as the mapped one.
  The mapped AIRR is now matched by exact name, and the contigs get their own
  `assembled_airr` channel (`optional: true`, for `--no-assemble`).

## 2.5.0

### Behaviour changes for existing users

**`d_call` output moves.** Three independent fixes changed which records get a D and which D
they get. If you depend on 2.4.0's D calls, re-run rather than diff.

* The V..J interior is now bounded by the **per-allele junction anchors** shipped in
  `database/vdj/<org>/cdr3_anchors.tsv`, not by the mmseqs scaffold projection. A scaffold has
  a 9 nt N-pad where a read has 20–40 nt of N-D-N, so the projection collapsed the very window
  the D lives in — on real human IGH, an 11 nt interior where the truth is 37 nt. This, not
  "paralogous germlines + SHM", was the cause of the old 46–69 % IGH D concordance. IGH D-gene
  agreement with IgBLAST is now **86–98 %** across the five organisms (human 97 %, rat 98 %,
  rhesus 94 %, mouse 94 %, rabbit 86 %).
* The call is gated on a **Karlin–Altschul E-value** (`_D_MAX_EVALUE = 0.2`), replacing four
  hand-tuned per-locus score floors. The E-value ships as the AIRR `d_support` column so a
  consumer can re-threshold. Some records that got a D at 2.4.0 correctly get none now: the
  best hit simply does not clear the gate.
* **TRBD2 is never called with a TRBJ1.** TRBD2 lies 3′ of the entire TRBJ1 cluster and V(D)J
  joining deletes the intervening DNA, so the pair cannot be produced. Unenforced, TRBD2
  (16 nt) outscored TRBD1 (12 nt) on noise and took 17 % of real human TRB J1-cluster D calls.

**Output schemas gain columns** (additive; consumers reading by name are unaffected):

* `rnaseq correct` clonotype table: `d_call`, `d2_call`, `d_support`, `d2_support`.
* `rnaseq assemble` AIRR: the four above plus `np1`, `np2`, `np3`.
* `--seqtype aa` annotation now populates the D columns, which were previously always empty.

**`rnaseq correct` row order is now deterministic.** Clonotypes tied on
`(duplicate_count, consensus_count)` used to fall back to read order, which comes from a
threaded mmseqs search — the same FASTQ produced the same rows in a different sequence each
run. Ties now break on `(junction, v_call, j_call)`.

### Added

* **`arda.cdr3fix` + `arda markup`** — mark up and repair a bare `(CDR3 aa, V, J, species)`
  record with no read behind it: which residues each germline templates, where the submission
  disagrees, how far the disagreement extends, and a conservative repair. Everything is
  *junction space* (Cys104 … Phe/Trp118, both anchors included), which is what VDJdb's `cdr3`
  column holds and is **not** arda's `cdr3` field. Emits a VDJdb-compatible `cdr3fix` JSON
  blob and a human-readable `--report`.
* **`arda.annotate.dmap.map_d_junction`** — D and tandem D-D on a bare nucleotide junction,
  with no mmseqs pass.
* **`arda.dpost.posterior_d`** and `arda markup --d-posterior` — the D gene and its position
  inferred from the junction *length*, which pins `insVD + |D surviving| + insDJ`. Places the
  D to a median 1–3 nt even when the protein shows none of it. Shipped for human IGH/TRB/TRD
  and mouse TRB; every other pair returns nothing rather than guessing.
* **D mapping on amino-acid input**, against each D germline's three translated frames. A call
  on ~36 % of real IGH records, agreeing with the nucleotide call on 98 % of them; mostly
  silent for the TR loci, whose D is too short to survive trimming into protein.
* **New reference artifacts**, per organism: `cdr3_anchors.tsv` (per-allele conserved Cys104 /
  Phe-Trp118 anchors, with a `status` that is flagged, never guessed) and `d_prior.tsv`
  (generative-model summaries used by `dpost`; derived, not measured — see `SOURCES.md`).
* **`examples/`** — a runnable tour, every artifact derived from data already in the repo and
  rebuilt by `python examples/regenerate.py`: one mRNA per locus; the two human reads (of
  7,341 across five organisms) that carry a tandem D-D; seven VDJdb records covering every
  junction-repair outcome; and a 1,035-read FASTQ that runs the whole bulk RNA-seq pipeline in
  about six seconds.
* `--map-d/--no-map-d` on `rnaseq run`, `rnaseq assemble` and `rnaseq correct`.

### Fixed

* **`cdr3fix` repairs always land on a canonical junction.** `cdr3_repaired` is accepted only
  when it opens with Cys104 and closes with Phe/Trp118; otherwise the submission comes back
  untouched. Trimming flanking framework (`_MAX_TRIM = 3`) is now budgeted separately from
  inventing germline (`_MAX_FIX = 2`), and a trimmed residue costs something — a free trim
  could tie the untrimmed alignment, win the tie-break, and eat the conserved Phe of clean
  short IGK junctions. On the committed 250-row VDJdb fixture: repair reproduced **100/100**,
  zero novel rewrites, and `good`/`vCanonical`/`jCanonical` agree with VDJdb on every row.
* `Cdr3Error.applied` claimed an edit had been written to the output when a `Failed*` side had
  in fact discarded the whole repair.
* `Cdr3Markup.v_end` counted residues of the *submitted* junction, so a V-side trim or insert
  shifted it off `cdr3_repaired` — and `dpost` slices the non-templated middle with exactly
  those coordinates.
* `Cdr3Markup.v_canonical` / `j_canonical` described the submission rather than the repair,
  disagreeing with VDJdb's own booleans on 76 of 250 fixture rows.
* `dpost` reported `entropy = -0.0` for a degenerate posterior.
* `examples/example.airr.tsv` had been stale for four release rounds — written with 49 AIRR
  columns while the schema grew to 83. It is now regenerated and tested.

### Removed

* The private per-locus D score floors (`_D_MIN_SCORE`, `_D2_MIN_SCORE`), superseded by the
  E-value gate.
* The `family*01` rung of the allele-resolution ladder: reachable for **0** genes across all
  five shipped organisms.

### Docs and tests

* `arda.cdr3fix`, `arda.dpost` and `arda.annotate.dmap` are in the Sphinx API reference; the
  zero-warning gate still passes.
* Corrected two false claims: `rnaseq assemble` is implemented and on by default, and
  `rnaseq run` is three stages / four files.
* Test suite: 132 → 228 test functions, **268 passing** (1 skipped: optional `airr` dep).
