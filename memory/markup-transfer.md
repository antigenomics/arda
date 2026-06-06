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

## TODO / known gaps
- **D-segment mapping** (NOT done): scaffolds are V×J only, so `d_call` and the D
  region inside CDR3 are not assigned. Future: after V/J transfer, align the CDR3
  interior against a D germline DB and emit `d_call`/`d_sequence_*`. Must handle
  **double D-D junctions** (D-D fusions, seen in IGH/TRD) — possibly >1 D segment
  per rearrangement (AIRR `d2_call`). See README TODO.
- Multiprocessing: mmseqs is threaded; transfer is serial (cheap). Chunked
  process-pool sharding for huge FASTQ still deferred.
- `productive` is heuristic (stop-free V..J span), not full AIRR productivity.
