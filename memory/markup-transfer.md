# Runtime markup transfer (Phase 2)

## Pipeline
`annotate_records` (annotate/mapper.py): read seqs → `mmseqs easy-search` query vs
`database/vdj/<org>/alleles{,.aa}.fasta` → best hit per query (max bits) →
`transfer_hit` projects reference region coords onto the query via the C++
`_markup.transfer_regions` → AIRR TSV.

- nt input → search-type 3, `markup.tsv` coords (nt space).
- aa input → search-type 1, `markup.aa.tsv` coords (aa space). Same projection code.

## C++ `transfer_regions` (src/_markup/markup.cpp)
Single walk over `qaln`/`taln`. Inputs/outputs 1-based closed (AIRR). Verified:
- full identity → coords unchanged.
- insertion in query within a region → region span absorbs the inserted bases;
  downstream regions shift by the insertion length (validated: +3nt in CDR3 → FR4
  shifts +3, FR4 stays `WGxG`).
- deletion → deleted ref positions contribute no query base.
- 5'-truncated query → uncovered regions return (-1,-1) → emitted as empty.
Self-hit (query == scaffold) reproduces reference coords (~209/210; rare 1-off at a
boundary, acceptable).

## AIRR record assembly (annotate/transfer.py)
- Region nt seq = `query[qs-1:qe]`; aa via translating query from FR1 query-start.
- `productive` = no stop in the translated V..J span.
- `junction` = CDR3 ± one codon (nt) / ± one residue (aa).
- `v_call`/`j_call` carry the (possibly comma-joined) ambiguous allele set from
  the deduped scaffold.

## mmseqs binary discovery
`$ARDA_MMSEQS` → `bin/mmseqs` → PATH. In normal use arda is installed *into* the
`arda` conda env where `mmseqs` is on PATH. When running from another env (e.g.
dev/base), set `$ARDA_MMSEQS` to the env binary.

## Done since v1
- **Reverse-complement** nt queries: `--strand 2`; reverse hits (qstart>qend) are
  re-oriented (revcomp query, remap coords), `rev_comp=T`. See [[mmseqs-params]].
- **Target-DB caching**: `_cached_target_db` createdb once under `data/mmseqs_db/`.
- **C++ seq primitives**: translate/detect_frame/reverse_complement/back_translate
  in `src/_markup/markup.cpp`, re-exported by `refbuild/translate.py` (mirpy-compatible).

## D-segment mapping (done)
After V/J transfer, `transfer._map_d` takes the V..J interior of the junction
(query coords `v_sequence_end+1 .. j_sequence_start-1`) and gapless-local-aligns
every locus D germline against it via the C++ `_markup.d_local_align` (match=+1/
mismatch=-1; mmseqs is unreliable at ~8-31 nt). Best hit (score≥`_D_MIN_SCORE`=6)
→ `d_call` + `d_sequence_start`/`d_sequence_end` (AIRR, 1-based closed, query
space). For D-D loci (IGH/TRD) a second non-overlapping D above `_D2_MIN_SCORE`=7
→ `d2_call`/`d2_sequence_*`; `np1`/`np2`/`np3` partition the interior between V,
the D(s), and J. D germlines ship in `database/vdj/<org>/d_germlines.fasta`
(`>locus|allele`, VDJ loci only), loaded into `Reference.d_germlines`; VJ loci get
no germlines so D mapping is skipped. `build._collect_d_germlines` writes the file
during refbuild (from `imgt.ungap_gene`). Concordance vs IgBLAST (committed
fixtures): TRB/TRD ~97% gene agreement among co-called; IGH ~46-69% (paralogous D
+ SHM). **Limitation**: a long junction can exceed what mmseqs aligns *through*,
so the projected interior collapses and D mapping silently no-ops (lowers recall,
never a wrong call). Tests: `tests/unit` (d_local_align), `tests/synthetic`
(single-D e2e + `_map_d` double-D logic + option toggle), `tests/realworld`
(per-org concordance), `tests/benchmark` (overhead).

**Optional**: `map_d` (default True) threads through `annotate_records`/
`annotate_file`/`adapter.annotate_sequences` and CLI `--map-d/--no-map-d`; gated in
`_annotate_chunk` (`dg = ... if seqtype=='nt' and map_d else None`). D mapping is a
short per-hit C++ local alignment, not a new mmseqs pass, so the cost is tiny:
benchmark (`test_d_mapping_overhead`, VDJ-only worst case, M3, 16 threads) shows
**+1.0% @ 2k, +1.3% @ 10k** wall time (~7-15 µs/seq). Negligible on real bulk
RNA-seq where most reads aren't VDJ hits.

## TODO / known gaps
- `productive` is heuristic (stop-free V..J span), not full AIRR productivity.

## Streaming I/O (done)
`annotate_file` streams the input in bounded chunks (`_CHUNK_SIZE=50k`) via a
background reader thread (prefetch queue, maxsize 2) that parses the next chunk
while mmseqs annotates the current (subprocess releases the GIL). Reference +
cached target DB load once and are reused across chunks. Memory is flat for
arbitrarily large FASTQ; output is written incrementally (`airr_header` +
`format_rows`). `annotate_records` (in-memory) shares the same `_annotate_chunk`
core. CLI `--chunk-size`. Multi-node sharding still TODO.
