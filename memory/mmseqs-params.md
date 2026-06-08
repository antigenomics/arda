# MMseqs2 parameters & locus handling

## One DB, all loci
The runtime target is a **single** mmseqs DB built from `database/vdj/<organism>/
alleles{,.aa}.fasta`, which contains scaffolds for *all* loci (IGH/IGK/IGL/TRA/
TRB/TRG/TRD). So one search annotates mixed **bulk RNA-seq** across every locus at
once; the per-query `locus` (AIRR field) comes from the best-hit scaffold.

## Precompiled (shipped) indexes
The mmseqs target DBs are **committed** under `database/vdj/<organism>/mmseqs/<nt|aa>/`
(createdb output + a `VERSION` marker), so annotation runs out of the box (~24 MB
total). `_cached_target_db` (annotate/mapper.py): prefer the committed DB **iff its
`VERSION` == local `mmseqs version`** (DBs are version-sensitive); else build once
into `data/mmseqs_db/<org>_<seqtype>` (private cache, never dirties git). `arda
build-index [--force]` (mapper.build_index) regenerates the shipped DBs for the
local version — a maintainer tool, deliberately NOT in setup.sh so end users don't
dirty the committed blobs. CI uses a different mmseqs build → exercises the fallback.

## Tuned defaults (annotate/mapper.py)
- nt: `--search-type 3`, `-s 7.0`, `--max-seqs 50`, `--strand 2` (both strands),
  `-a`; no coverage filter (partial RNA-seq reads must still map).
- aa: `--search-type 1`, `-s 7.0`, `--max-seqs 50`, `-a`.
- Pipeline is `createdb (cached target) + createdb query + search + convertalis`,
  NOT `easy-search`, to reuse the target DB across calls.

Rationale: queries are germline-similar (85–95% id) and short; `-s 7.0` gives
reliable best-hit recall without the cost of 8.5. `--max-seqs 50` is plenty for a
single best germline. Default `-s 5.7` also worked (~97%); 7.0 buys margin on
divergent/short reads. Measured: 98.7% region concordance vs IgBLAST on 100 real
IGH mRNA; 98.9% on 10k synthetic.

## convertalis gotcha
`convertalis` needs `--search-type 3` for nucleotide results too, or it errors
"unclear if translated or nucleotide search". Wired via `mmseqs.convertalis(search_type=)`.

## Reverse strand
`--strand 2` finds reverse hits, reported with **qstart > qend** and qaln/taln
already on the coding strand. mapper detects this, reverse-complements the query,
remaps qstart/qend onto it, and sets `rev_comp=T`. CLI `--strand both|forward`.

## Speedup vs IgBLAST (synthetic IGH, 16 threads)
~4.4× (10k) → 7.9× (100k); arda ~3.3k seq/s vs igblast ~0.42k seq/s. Throughput
rises with N as fixed DB-build cost amortizes. See `scripts/bench_vs_igblast.py`.
