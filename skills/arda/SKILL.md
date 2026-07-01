---
name: arda
description: >
  Fast TCR/BCR FR/CDR region annotation (Antigen Receptor Domain Annotation).
  Use whenever the user wants to annotate immune-receptor sequences with
  framework/CDR regions, V/D/J gene calls, junction/CDR3 boundaries, or AIRR
  output — for TCR (TRA/TRB/TRG/TRD) or BCR (IGH/IGK/IGL), nucleotide or amino
  acid, single sequences or large FASTA/FASTQ. Also use for: getting germline
  FR1-FR3/CDR1-CDR2 (V) or FR4 (J) subsequences for individual alleles; building
  or rebuilding the reference database from IMGT germlines via IgBLAST; or
  diagnosing mmseqs2 setup. arda runs MMseqs2 + a C++ coordinate projection
  (IgBLAST is offline/build-time only), so it is much faster than IgBLAST at
  annotation time. Load references/ files for the detailed API, region/junction
  semantics, reference-build pipeline, or mmseqs install/troubleshooting.
license: GPL-3.0
compatibility: >
  Python 3.10+; `pip install arda-mapper` (>=2.0.3) — no source checkout and **no
  `ARDA_HOME`**: the curated `vdj/` reference auto-fetches into `~/.cache/arda` on first use
  (set `ARDA_NO_AUTO_FETCH` for air-gapped runs), and the `mmseqs` binary auto-fetches a
  static build into the cache if missing — so a bare `pip install` annotates out of the box.
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
onto each query. That makes it embeddable and ~4–8× faster than IgBLAST with
~97–99% concordance.

## Core API

```python
import arda

records = arda.annotate_sequences(
    ["GACGTGCAG...", ("clone7", "CAGGTG...")],  # raw strings or (id, seq) pairs
    seqtype="nt",          # "nt" or "aa"
    organism="human",      # human | mouse | rat | rabbit | rhesus_monkey
    map_d=True,            # map D segments for VDJ loci (nt input only)
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
`v_call`/`d_call`/`d2_call`/`j_call`, `productive`, `rev_comp`, `v_sequence_end`,
`j_sequence_start`, `np1/np2/np3`, `junction(_aa)`, and per region in
`(fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4)`: `{r}_start`, `{r}_end`, `{r}`,
`{r}_aa`.

Read [references/annotation.md](references/annotation.md) for the full field list,
parameter semantics (strand/sensitivity/threads/chunking), AIRR column order, and
performance notes.

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

## Organisms & loci

| Organism | Loci with full markup |
|----------|-----------------------|
| human, mouse | TRA, TRB, TRG, TRD, IGH, IGK, IGL |
| rat, rabbit, rhesus_monkey | IGH, IGK, IGL (IG only) |

VDJ loci (D segments mapped, nt only): IGH, TRB, TRD. D-D fusions: IGH, TRD.

## CLI

```bash
arda info                                   # versions + available references
arda annotate -i reads.fastq.gz -o out.airr.tsv --organism human --seqtype nt
arda annotate -i prot.fasta -o out.tsv --seqtype aa --no-map-d
arda build-db --organism all                # offline reference build (needs IgBLAST)
arda build-index --organism all             # rebuild mmseqs indexes for local mmseqs version
arda slurm -i big.fastq -o big.airr.tsv --shards 50   # multi-node: split → array → merge
```

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
organism (needs IgBLAST, fetched by `setup.sh` into `bin/`).

Read [references/reference-build.md](references/reference-build.md) for the
`arda.refbuild` pipeline (IMGT germlines → V×J scaffolds → IgBLAST → markup TSVs)
and `build-db` / `build-index`.

## Sequence primitives

`arda.refbuild.translate` exposes fast C++-backed helpers, mirpy-API-compatible:
`translate(nt, frame=0)`, `detect_coding_frame(nt)`, `reverse_complement(nt)`,
`back_translate(aa)`, `aa_coords_from_nt(nt_start, nt_end, coding_start)`.

## Gotchas

- D mapping and `productive` only populate for nt input; aa input returns region
  `*_aa` directly with no frame bridging.
- The shipped MMseqs2 indexes are used only when the local mmseqs **version**
  matches; otherwise arda rebuilds a private cache in `data/` on first run.
  `arda build-index` (re)builds the shipped indexes for your version.
- `map_d=True` on synthetic/partial input with no real junction simply finds no
  D — harmless; pass `map_d=False` to skip the search.
- IgBLAST is needed only to build references, never at annotation time.
