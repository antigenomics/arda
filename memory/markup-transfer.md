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

## TODO / known gaps
- Reverse-complement nt queries: `rev_comp` is hard-coded "F"; mmseqs nt search
  strand handling not yet wired (detect qstart>qend / use --search-type strand).
- No mmseqs target-DB caching yet (easy-search rebuilds each call; target is small).
- Multiprocessing: mmseqs is threaded; the transfer step is serial (cheap). Chunked
  process-pool sharding deferred to the benchmark phase.
