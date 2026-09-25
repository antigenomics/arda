---
name: arda
description: >-
  Annotate TCR/BCR sequences with FR1-FR4 / CDR1-CDR3 regions, V/D/J gene calls, constant-region
  isotype and junction boundaries, emitting spec-valid AIRR. Use for nucleotide or amino-acid
  input, single sequences or whole FASTQ; for extracting a repertoire from bulk RNA-seq, amplicon
  or single-cell libraries; for marking up and repairing a bare (CDR3 amino acid, V, J) record
  that has no read behind it, as in VDJdb; for pulling germline FR/CDR subsequences per allele;
  for run and batch quality control over one sample or a whole cohort; for inferring which V
  alleles a donor carries and restricting calls to a personalized germline; or for rebuilding the
  reference database from IMGT.
---

# arda

arda does the expensive IgBLAST work **once, offline** — a reference of every in-frame V·J
germline scaffold with FR1–4 / CDR1–3 markup — then at runtime maps queries to it with MMseqs2
and projects the markup through the alignment in C++. That makes it embeddable and ~4–8× faster
than IgBLAST, with 98–99.7 % region concordance on real GenBank mRNA across all five organisms.

It also handles records with **no read behind them** — a CDR3 amino acid plus a V and J call, as
in VDJdb — marking up which residues each germline templates, repairing the junction, and
inferring the D gene from the junction's length.

## Install and environment

`pip install arda-mapper` (imports as `arda`). Nothing else is needed: the curated `vdj/`
reference auto-fetches into `~/.cache/arda` on first use, and a static `mmseqs` binary fetches
itself if one is not already resolvable. A source checkout uses the committed `database/` instead.

- `ARDA_NO_AUTO_FETCH` — air-gapped runs against a pre-populated cache.
- `ARDA_MMSEQS` — pin a binary; overrides everything and is not version-checked.
- `ARDA_HOME` — a source checkout; **disables auto-fetch**. To simulate a pip user, set
  `XDG_CACHE_HOME` instead (the cache is `$XDG_CACHE_HOME/arda/database`).
- `arda info` prints resolved paths and tool availability — run it first when anything looks wrong.

Bulk RNA-seq needs no extra: `seqtree` and `dnaio` are core dependencies. The `rnaseq` and
`mmseqs` extras are empty aliases kept so old pins still resolve.

## Core API

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],  # raw strings or (id, seq) pairs
    seqtype="nt",          # "nt" or "aa"
    organism="human",      # human | mouse | rat | rabbit | rhesus_monkey
    map_d=True,            # map D segments for VDJ loci — works on aa input too
)
```

Each record is a dict of AIRR fields in query space, 1-based closed: `locus`,
`v_call`/`d_call`/`d2_call`/`j_call`, `c_call`/`c_class`, `productive`/`stop_codon`/`vj_in_frame`,
`rev_comp`, `v_identity`, `sequence_alignment`/`germline_alignment`, `{v,j,c,d}_cigar`,
`*_germline_start`/`*_germline_end`, `v_sequence_end`, `j_sequence_start`, `np1/np2/np3`,
`d_support`/`d2_support`, `junction(_aa)`, and per region in
`(fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4)`: `{r}_start`, `{r}_end`, `{r}`, `{r}_aa`. Ambiguous
D and C calls are comma-joined allele lists, as V and J are. The TSV is a spec-valid AIRR
Rearrangement file and passes `airr.schema` validation.

`d_support` is the Karlin–Altschul **E-value the D call was gated on** (accepted at `<= 0.2` for
nt, `<= 0.05` for aa). It ships so a consumer can re-threshold: keeping rows with
`d_support <= x` for `x < 0.2` reproduces exactly what a stricter arda would have called. A
missing `d_call` on a VDJ locus usually means the best hit did not clear the gate, not that
mapping was skipped. For explicit strand / sensitivity / streaming control use
`arda.annotate.mapper.annotate_records` and `annotate_file`.

→ **[references/annotation.md](references/annotation.md)** — full field list, parameter
semantics, AIRR column order, the D E-value gate and the genomic-order constraint.

## Never loop — mmseqs2 is the parallel layer

**Gather every sequence, make ONE `annotate_sequences` call, then work on the batch output.**
Each call pays a fixed ~825 ms mmseqs2 process + index-load cost; 300 sequences cost the *same*
~930 ms total because mmseqs2 threads internally.

```python
recs = arda.annotate_sequences([(cid, seq) for cid, seq in all_chains], organism="human")
by_id = {r["sequence_id"]: r for r in recs}     # map back per-item downstream
```

Never wrap a per-item `annotate_*` in a loop, a `ProcessPoolExecutor` or a `ThreadPoolExecutor`:
a process pool that forks after mmseqs2/BLAS have spawned threads **deadlocks**, a thread pool
just serialises on the same overhead, and either way the fixed cost is paid N times instead of
once.

## Region and junction semantics

- Coordinates are projected through the alignment, so they are correct for truncated, mutated and
  reverse-strand queries alike.
- There is **no coverage filter**: a partial read — or a bare germline V or J — maps to its
  scaffold and returns only the regions inside its coverage. A bare V → `fwr1..fwr3`; a bare J →
  `fwr4`. That is how to get per-allele germline FR/CDR subsequences without inventing a
  rearrangement.
- `junction` spans Cys104 through the [FW]118 that opens FR4; `cdr3` is J-anchored. Out-of-frame
  junctions are reported with an N-bridge (`_`).

→ **[references/region-segments.md](references/region-segments.md)** — the bare-germline recipe,
junction/CDR3 detail, coordinate round-trip rules.

## Organisms and loci

| Organism | Loci with full markup |
|----------|-----------------------|
| human, mouse | TRA, TRB, TRG, TRD, IGH, IGK, IGL |
| rat, rabbit, rhesus_monkey | IGH, IGK, IGL (IG only) |

D segments are mapped for the VDJ loci (IGH, TRB, TRD), D-D fusions sought in all three, on
protein input as well as nucleotide. rat/rabbit/rhesus have no TR loci in IMGT, so those build
`EMPTY` and say so in `loci_manifest.tsv` rather than failing silently.

## CLI

```bash
arda info                                   # resolved paths + tool availability
arda annotate -i reads.fastq.gz -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta -o out.tsv --seqtype aa --no-map-d
arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report -   # bare (CDR3aa, V, J) records

arda rnaseq   --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # bulk: map+assemble+correct
arda amplicon --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE -d out/   # targeted RepSeq, other preset
arda cells    asm/PBMC.consensus.fq.gz -p out/PBMC            # single cell (UMI consensus in)

arda map      --r1 R1.fq.gz --r2 R2.fq.gz -o mapped.airr.tsv  # the stages, separately
arda assemble -i mapped.airr.tsv -o assembled.airr.tsv
arda correct  -i mapped.airr.tsv --extra-airr assembled.airr.tsv -o clones.tsv
arda shm      -i mapped.airr.tsv -o rescoped.airr.tsv         # recount SHM outside the junction
arda stats    -i mapped.airr.tsv -c clones.tsv -r SAMPLE.arda.json -o SAMPLE.stats.tsv
arda qc batch  -d out/ -o out/batch --samples sheet.tsv   # every sample's QC as one table
arda qc report -i out/batch.qc.json -o out/batch.qc.html  # ...as one self-contained page
arda scenarios -i clones.tsv -o d_prior.tsv   # EM over recombination scenarios -> a generative model

arda genotype     -i mapped.airr.tsv -o donor.genotype.tsv --loci TRB   # which V alleles this donor has
arda resolve-ties -i mapped.airr.tsv -o narrowed.airr.tsv --genotype donor.genotype.tsv

arda igblast    -i reads.fastq -o truth.airr.tsv              # gold-standard IgBLAST, all loci
arda export-ref --kind segments --locus TRB --format fasta
arda build-db    --organism all             # offline reference build (needs IgBLAST)
arda build-index --organism all             # rebuild mmseqs indexes for the local mmseqs
arda cluster submit --r1 R1.fq.gz --r2 R2.fq.gz -p SAMPLE --shards 8    # one big pair, sharded
arda cluster plan --samples sheet.tsv --work-dir work/                  # any scheduler
arda cluster submit-samples --samples sheet.tsv -d out/                 # two SLURM arrays
```

**Never: `-v` / `-q` / `--log-file` are GLOBAL and go BEFORE the subcommand**
(`arda -v --log-file run.log amplicon --r1 ...`). Progress goes to **stderr**, results to
**stdout** — a mode run prints its output paths one per line and nothing else, so `$(arda map
...)` and `arda export-ref ... > out.tsv` are safe. `--log-file` is always DEBUG whatever the
console level is, and `-q` does not silence it.

**Never: `correct` needs `--extra-airr`.** Without it the contigs `assemble` just built are
silently discarded and the clonotypes no single read spans never reach the table. The mode
commands wire it for you.

## Pick the mode, not the flags

| mode | library | configuration it carries |
|---|---|---|
| `arda rnaseq` | bulk RNA-seq, 0.02–3 % receptor | `--prefilter` |
| `arda amplicon` | targeted RepSeq / 5'RACE | `--two-pass --fast-segments --v-only-on-segment` |
| `arda cells` | single cell, **one UMI consensus per molecule** with the barcode in the record name | reference-free per-cell assembly, then pairing |

The two speed configurations **do not compose** and each is a loss in the other's regime, which is
why the mode name owns the preset. `arda cells` starts one step later than the other two: arda
does no demultiplexing, no barcode correction, no UMI collapse and no cell calling — the upstream
tool (e.g. `migec`) does all four and leaves the barcode in the record name.

→ **[references/speed-flags.md](references/speed-flags.md)** — each flag's measured regime, the
segment reference and rescue guarantee, and the internals a hot-path change must respect.
→ **[references/rnaseq-pipeline.md](references/rnaseq-pipeline.md)** — what each stage owns, the
run report, `arda stats` and `arda qc`.
→ **[references/read-groups-and-cluster.md](references/read-groups-and-cluster.md)** — a sample
split across several FASTQs, and the three cluster adapters.
→ **[references/bare-records.md](references/bare-records.md)** — `arda markup`, junction repair
and the D posterior.
→ **[references/genotype.md](references/genotype.md)** — `arda genotype` and
`resolve-ties --genotype`: the personalized germline, why it never re-aligns, and why read length
decides what it can say.
→ **[references/install-mmseqs.md](references/install-mmseqs.md)** — mmseqs env vars, the shipped
indexes, version-mismatch handling.
→ **[references/reference-build.md](references/reference-build.md)** — the `arda.refbuild`
pipeline and `build-db` / `build-index`.

## Sequence primitives

`arda.refbuild.translate` exposes fast C++-backed helpers, mirpy-API-compatible:
`translate(nt, frame=0)`, `detect_coding_frame(nt)`, `reverse_complement(nt)`,
`back_translate(aa)`, `aa_coords_from_nt(nt_start, nt_end, coding_start)`.

## Gotchas

- **`junction` is not `cdr3`.** Both conserved anchors are *in* `junction`/`junction_aa` and *out*
  of `cdr3`/`cdr3_aa`, so `cdr3_aa == junction_aa[1:-1]` always. Everything in `arda.cdr3fix` /
  `dmap` / `dpost` works in **junction** space, matching VDJdb's `cdr3` column. Mixing the two
  conventions is the most expensive mistake available here, and it corrupts Pgen, clustering and
  matching downstream.
- **An empty `d_call` is a decision, not a gap.** The call is gated on `d_support` (E-value ≤ 0.2
  nt, ≤ 0.05 aa). Human TRB gets a D on only ~47 % of junctions because an ordinary TRB interior
  is 11–21 nt and a heavily-trimmed TRBD1 scores below the gate. TRA, TRG, IGK and IGL have no D
  gene at all.
- **A conserved-motif check is not an anchor.** `TRAJ35*01`'s anchor codon decodes Cys (TGC), not
  [FW] — it is a functional IMGT `F` gene. Read `anchor_nt` from `cdr3_anchors.tsv`; a `[FW]GXG`
  motif check silently deletes the gene.
- **An overlapping V / J / NDN assignment inside the junction is not a defect.** Exonuclease
  chew-back and N/P addition mean that partition is often not identifiable from sequence alone.
  What *is* checkable: the junction's outer bounds, the gene calls, and whether arda invented a
  junction it has no anchor for.
- `posterior_d` returns `None` for organisms with no shipped generative model (rat, rabbit,
  rhesus) and for VJ loci. That is deliberate — do not substitute a human proxy. To score one of
  those pairs, fit a table with `arda scenarios` and pass it: `posterior_d(..., prior_path=)` /
  `arda markup --d-prior PATH`. Nothing in the installed database is touched.
- **An IG `v_call` is only as good as the V germline the read covers**: `v_gene` recall is
  **.1170 under 60 nt and .9896 at 200 nt or more**, and .9896 is the TRA amplicon's .9867 — no
  IG-specific deficit at that coverage. Position beats length (IGHV diverges in FR1/CDR1/CDR2,
  conserved near Cys104). No threshold ships; filter `v_germline_start`/`v_germline_end` yourself.
  Tables and the rejected gate: `docs/usage.rst`.
- aa input returns region `*_aa` directly with no frame bridging, so `stop_codon` and
  `vj_in_frame` stay empty — but `productive` and the D columns *are* populated.
- `map_d=True` on synthetic or partial input with no real junction simply finds no D. Harmless;
  pass `map_d=False` to skip the search.
- The shipped MMseqs2 indexes are used only when the local mmseqs **version** matches; otherwise
  arda rebuilds a private cache on first run. `arda build-index` rebuilds the shipped ones.
- IgBLAST is needed only to build references and to run `arda igblast`, never at annotation time.
  It auto-fetches on first use like mmseqs does.
