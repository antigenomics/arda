# Changelog

Notable changes per release. Earlier releases are described by their git tags
(`git tag --sort=-v:refname`); this file starts at 2.5.0.

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
