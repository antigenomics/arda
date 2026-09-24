# Quality control in arda: what the table could not say, and what a cohort needs

2026-09-24. Shipped in 2.24.0. Supersedes nothing -- `stats.py` (2.14.0) is the thing this
extends, not replaces.

## The premise this starts from, and it is not "arda has no QC"

**"arda needs a QC report."** It has had one since 2.14.0. Every `arda rnaseq` / `arda amplicon`
run writes `<prefix>.stats.tsv` unconditionally -- `pipeline.write_stats_for:140`, not behind a
flag -- with about 60 metrics over five scopes, and a `docs/qc.rst` page describing them.
The gap is narrower and more specific than "no QC", which is why the work is small.

**"and then someone joins the tables across samples."** Nobody does. `stats.py:16-23` says in its
own docstring that long format was chosen because "the one thing a QC table must support is
`grep` / `join` / a per-metric plot across samples", and `docs/qc.rst:76-78` repeats it. Grep for
a consumer: no code globs `*.stats.tsv`, `*.arda.json` or `*.clones.tsv` across samples;
`cli.py`'s `_mode_run:896` loops samples and emits no combined artifact;
the Snakemake workflow's `rule all` only `expand`s per-sample paths. The word
"cohort" appears once in `src/`, in a comment at `pipeline.py:344` justifying a column name. The
format was built for a join that was never written.

## What is actually missing

1. **Shape.** The table carries min, max and mean. A bimodal junction length, a read-length cliff
   and a clone-size distribution with no singletons all read as an ordinary mean, and all three
   are the diagnosis.
2. **Single cell.** `arda cells` writes no `.stats.tsv` at all, and its `.arda.json` shares zero
   keys with the bulk one. A cohort mixing the two has nothing to join on.
3. **The join**, plus two defects that would corrupt it rather than fail (below).
4. **Anything that renders it.** A bulk run produces no figure of any kind.

## The scope line, and why it is where it is

`~/vcs/code/vdjtools` (v2, Python) already ships `stats/{diversity,rarefaction,inext,spectratype,
usage,functional,shm}`, `overlap/{metrics,similarity,track,tcrnet}`, `preprocess/decontaminate.py`
-- the cross-sample index-hopping filter -- and `io/cohort.py`, hive-partitioned Parquet for
cohorts larger than RAM. **It takes `arda-mapper` as a base dependency**, and its own CLAUDE.md
says "Delegate rather than reimplement ... annotation/markup -> arda". Building diversity here
would be a library reimplementing its own dependent.

arda's QC answers **"did this run work, and is this sample like its batch"**. Never "what is this
repertoire".

| not arda's | whose | why |
|---|---|---|
| diversity, clonality, rarefaction, D50, iNEXT | `vdjtools.stats` | estimators, not run QC, and already shipped downstream |
| overlap, public clonotypes, TCRnet | `vdjtools.overlap` | needs cross-sample clonotype matching |
| contamination / index-hopping removal | `vdjtools.preprocess.decontaminate` | same, and it *filters*; arda flags |
| spectratype as an analysis | `vdjtools.stats.spectratype` | arda emits the raw junction-length histogram as a QC scope only |
| per-base quality, GC, primer and barcode QC | pRESTO, VDJPipe, fastp | arda never sees a pre-filter read |
| the MiAIRR read-QC fields | the sequencing facility | see below |

The line: a histogram from one `group_by` over a column arda already wrote is QC. An estimator
over that histogram is biology.

**Never: the MiAIRR QC fields cannot be filled honestly from here.**
`total_reads_passing_qc_filter` is the facility's usable-read count, upstream of annotation.
arda's `mapped_reads` is reads that matched a receptor -- 0.024 % of a bulk RNA-seq library is
normal and is not a QC failure, and `total_reads` is whatever arda was handed, with no knowledge
of what was filtered before. Emitting either under the MiAIRR name would assert something arda
does not know. `software_versions` and `germline_database` do have honest sources and are already
in the `run` scope under arda's own names; renaming them buys nothing. So no MiAIRR aliases ship,
and a full `Repertoire` / `DataProcessing` document waits for someone actually submitting to a
repository.

## S0 -- distributions, and one address per fact

Five scopes, all keyed `locus:bucket` so one parser reads them all, all **sparse** -- only
occupied buckets get a row, the rule `v_gene` already used. An amplicon occupies ~30 of ~200
possible junction lengths.

`junction_aa_len` (spanning reads only), `read_len` (10-nt buckets of `mmseqs2_qlen`),
`clone_size` (powers of two, weighted by clonotypes *and* by reads), `isotype` (`c_call`), and
`chain_support` from S1. Plus `reads_rev_comp` and `reads_junction_completed`.

Note what that costs: **nothing new is computed.** `mmseqs2_qlen`, `rev_comp`,
`junction_completed_nt` and `c_call` were all already being read into `stats.py` and discarded --
`_READ_COLS` listed `mmseqs2_qlen` at line 62 and no aggregate ever used it.

Never: **`mmseqs2_qlen` is a float string.** mmseqs writes `407.0`, and `cast(pl.Int64,
strict=False)` on `"407.0"` is null, not 407. Every bucket came out 0, `_hist` filtered them all
out, and the distribution vanished silently rather than failing. Cast via `Float64`.

Three fixes to the run scope, each of which corrupts a join rather than breaking it:

- **`_merge_map_reports` renames the facts a lane-split sample reports.** `wall_seconds` becomes
  `wall_seconds_max`/`_sum` and `peak_rss_mb` becomes `peak_rss_mb_max`. So a join on
  `run/map/wall_seconds` drops exactly the samples delivered as several FASTQ pairs, with no
  error anywhere. The canonical name is restored in `_flatten_report` -- the one place a report
  becomes QC rows -- leaving `.arda.json` untouched. `wall_seconds_sum` is a different fact and
  keeps its name.
- **`segment_search.reasons` never arrived.** `_flatten_report` flattened one level and dropped
  any value that was itself a dict. `reasons` is the only evidence for why `fast_fraction` is low.
- **An empty filtered set was written blank, not omitted.** `chain/IGH/junction_nt_min` appeared
  with an empty value when no IGH read spanned both anchors, and an empty cell casts to 0 --
  "the shortest IGH junction was 0 nt" instead of "there wasn't one". `docs/qc.rst:80-81` already
  promised omission.

`write_stats_json` emits the same rows as typed, nested JSON, built from the same `rows` list so
the two cannot diverge.

Acceptance: a multi-read-group sample and a single-file sample produce the same `run/map` metric
names, and no row anywhere has an empty value.

## S1 -- single cell in the same scopes

`collect(cells=<prefix>)` reads the single-cell report and `.chains.tsv` into the **bulk** scopes:
`chain` per locus, `v_gene`/`j_gene` with coverage, `junction_aa_len`. "Which loci did this sample
yield, at what junction lengths" is the same question either way. What is genuinely single-cell --
cells, molecules, knee, contig N50 -- comes through `run` from its own report, where it cannot be
mistaken for a read count.

Never: **a chain the extra-chain gate marked `extra` is not yield.** `min_extra_molecules` exists
because an extra chain supported by one molecule is ambient 96-97 % of the time; counting those
would report the contamination as signal, which is the one thing the gate is for.

`pairing_rate`, `doublet_rate` and `molecules_placed_fraction` are derived rather than left to the
reader: the first two are the AIRR Community's chain-pairing QC, they are what a batch is compared
on, and `cell_summary:789` has already assigned every cell the status they count.

Also here: `summary["figures"]` was assigned at `singlecell.py:1130`, **after** the JSON was
written at `:1127`, so the one key naming the panels a run drew never reached disk.

Acceptance: a single-cell prefix and a bulk prefix produce overlapping scope sets, and
`chain/<locus>/chains` excludes every `extra` row.

## S2 -- the roll-up

`arda qc batch -d results/ -o batch [--samples sheet.tsv]`. `src/arda/qc.py`: `collect_batch`,
`pivot`, `outliers`, `write_batch`.

**It reads only `*.stats.tsv`.** Never an AIRR, never a clonotype table. A 1,000-sample cohort is
~20 MB and one `pl.concat`; it runs on a laptop against results copied off a cluster with the bulk
data left behind, and it cannot become the slow step. Same split as Stage 1 (per read, shards)
against Stages 2-3 (global, do not): the expensive reduction already happened once, per sample,
where the data was.

The sheet gains optional `project` and `batch` -- labels, read from a sample's first row, that
nothing in the pipeline looks at. They were previously warned about as unknown columns.

**No shipped threshold.** Within each `(project, batch)` group, per metric: `median`, `mad`, and
`0.6745·(x − median)/mad`, flagged at `|z| ≥ 3.5` (Iglewicz-Hoaglin). Reported, never
applied. A constant saying what a good `mapped_fraction` is would depend on the library, the
organism and the depth, and one that looked calibrated would be worse than none -- the same
reasoning that shipped `chain_fraction` unset in `design-singlecell.md`.

Never: **a group under 5 samples gets no z.** MAD over four points is not a scale estimate, and a
z derived from one flags whichever sample happens to sit furthest from the middle. Metrics whose
MAD is 0 get none either, and non-numeric values never acquire a median.

Note: **an empty `key` round-trips as null.** An unkeyed scope writes an empty cell and polars
reads an empty cell as null, so `sample`'s key comes back as `None` unless filled -- and a join on
`key` then matches nothing for the scope carrying most of the metrics. `collect_batch` fills it.

Acceptance: every per-sample row in the batch table equals that sample's own `stats.tsv` row,
field for field.

## S3 -- one HTML file

`arda qc report -i batch.qc.json -o batch.qc.html`, or a single sample's `.stats.json` -- one run
is a cohort of one, lifted into the batch shape so there is one renderer and one input format.

Never: **no CDN, no bundle, no dependency.** The data is inlined and ~350 lines of plain
JavaScript draw SVG. A chart library would be richer and would cost the one property the file
exists for: it opens on an air-gapped login node, off a USB stick, or as an email attachment, and
keeps working after the results directory is gone. Same stance as `scplot.py`, whose gnuplot
script is written whether or not gnuplot exists. Palette and ink are `scplot`'s, not a second one.

The robust z is recomputed **in the browser** rather than shipped per cell: it keeps the JSON
small and guarantees the shading agrees with the numbers beside it.

Never: **`createElementNS(null, 'tbody')` is not an HTML tbody.** It is an element in the null
namespace, which the layout engine draws as nothing. The table filled correctly, held six rows in
the DOM, and rendered as a bare header. Found by opening the page in a real browser, which is why
that step is in the verification list and not optional.

Acceptance: the rendered page contains no `<script src>`, no `<link href>`, no `@import`, no
`url(http…)` and no `fetch(` -- the SVG namespace URI excepted, which browsers never retrieve.

## Order

S0 -> S1 -> S2 -> S3, each mergeable alone. S0 is the only one with a defect in it; the other
three are additive and touch nothing an existing run reads.
