# ARDA candidate and allele contracts

## Preliminary alignment candidates

`arda.alignment_candidates([(id, coding_strand_nt), ...], organism="human",
threads=0, max_seqs=300, evalue=1e-3, sensitivity=7.0)` returns flat Polars `queries` and `alignments` frames plus explicit `limits`.
Each call accepts at most `batch_size=10000` records; iterate bounded batches.
Candidates retain segment, locus, allele groups, query/reference 1-based closed
coordinates, CIGAR, mismatches, aligned strings, source and alignment score.
One batched V/J segment search retains MMseqs hits before best-hit reduction;
D retains the best positive gapless local span per allele in one native
batch. These are preliminary evidence, never fixed assignments or recombination
boundaries. MMseqs target/evalue limits and local-span limits mean this is not
exhaustive scenario enumeration; probabilistic consumers must reconsider trimming,
insertions and mismatches. Inputs must be nucleotide coding-strand sequences.
The `threads` budget applies sequentially to MMseqs and the GIL-released native
D batch. D workers share immutable germlines and own contiguous query spans;
queries preserve input order; alignments sort by sequence ID, segment, target,
query start, target start, query end and target end at every thread count. Batches below 32 queries
run serially to avoid startup cost; `limits` declares this threshold and actual
D worker count. Thread startup failure raises instead of falling back.

This API rejects `ARDA_MMSEQS_SEARCH_OPTS` overrides so declared limits match
the search. Invalid query IDs/lengths/bounds raise; reverse/inverted hits are
excluded with `limits["rejected_reverse_hits"]` reporting their count.

## Major alleles

Annotation/search defaults to `major_only=True` on `annotate_sequences`,
`annotate_records`, `annotate_file`, `alignment_candidates`, `load_reference`,
and `map_rnaseq`. CLI `annotate` and `map` default to `--major-only`; select
`--all-alleles` (or `major_only=False`) for the full frozen reference.
Searches use actual major templates; minor alignment coordinates are never relabelled.
Between-gene alternatives remain available. Full base FASTAs/markup stay unchanged.

`arda.select_major_alleles(eligible_names)` returns `gene -> allele`, preferring
`*01`, otherwise the lowest numeric suffix; starless D names work. Callers own
complete/functionality eligibility. `arda.major_allele_map(organism="human",
locus=None)` describes annotation majors: prefer available ok F/ORF anchors,
otherwise keep a preliminary template for a gene with only partial/no-anchor records.
Its partial status never establishes complete emission eligibility. Candidates
return the selected map as `major_alleles`. Models must filter incomplete templates
before calling the name-only selector. Frozen resource export/raw metadata helpers
explicitly use the full reference.

`arda.d_germline_metadata(organism="human")` returns metadata keyed by
`(allele, sequence)`: raw completeness, status and functionality. D FASTA
alignment targets remain available; complete emission templates exclude partial
records. `Anchor.imgt_completeness` preserves V/J raw completeness labels.
