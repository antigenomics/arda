# arda roadmap

Implemented: offline V·J reference build (5 organisms), MMseqs2 runtime mapping,
C++ markup transfer, reverse-complement handling, all-loci single-DB querying,
streaming/bounded-memory FASTQ I/O, out-of-frame junction translation, extended
V/J-position markup, offline GenBank-vs-IgBLAST test fixtures.

## TODO

- [ ] **D-segment mapping.** Scaffolds are V·J only, so `v_call`/`j_call` and all
      FR/CDR coordinates are assigned but `d_call` is not. Plan: after V/J transfer,
      align the CDR3 interior (between the projected CDR3 start and end) against a D
      germline DB and emit `d_call` + `d_sequence_start`/`d_sequence_end`.
  - [ ] **Double D-D junctions.** Handle rearrangements with two D segments (D-D
        fusions, seen in IGH/TRD): emit a second `d2_call` + `d2_sequence_*` and
        the intervening N regions (`np1`/`np2`/`np3` in AIRR terms).

- [ ] **Multi-node sharding.** Single-node streaming + threading is implemented;
      add a SLURM-array mode that shards a huge FASTQ across array tasks and
      concatenates the per-shard AIRR TSVs.

- [ ] **Full AIRR productivity.** `productive` is currently a heuristic (in-frame
      + stop-free V..J span); align it with the complete AIRR productivity rules
      (start codon, stop-codon scan over the whole VDJ, frame of the junction).

- [ ] **Performance.** Optional per-chunk process-pool for inputs where mmseqs is
      not the bottleneck (mostly-non-receptor bulk RNA-seq); mmseqs index reuse.
