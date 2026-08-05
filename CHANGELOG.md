# Changelog

Notable changes per release. Earlier releases are described by their git tags
(`git tag --sort=-v:refname`); this file starts at 2.5.0.

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
