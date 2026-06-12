# arda roadmap

Implemented: offline V·J reference build (5 organisms), MMseqs2 runtime mapping,
C++ markup transfer, reverse-complement handling, all-loci single-DB querying,
streaming/bounded-memory FASTQ I/O, out-of-frame junction translation, extended
V/J-position markup, D-segment mapping (incl. D-D fusions), offline
GenBank-vs-IgBLAST test fixtures.

## TODO

- [x] **D-segment mapping.** After V/J transfer, the V..J interior of the junction
      (between the projected `v_sequence_end` and `j_sequence_start`) is aligned
      against the per-organism D germline set by gapless local alignment in the C++
      `_markup.d_local_align` primitive — mmseqs is unreliable on ~8-31 nt D — and
      the best hit is emitted as `d_call` + `d_sequence_start`/`d_sequence_end`
      (AIRR, query coords). D germlines ship in `database/vdj/<org>/d_germlines.fasta`
      (VDJ loci only); VJ loci are skipped automatically. Concordance vs IgBLAST:
      TRB/TRD ~97% gene agreement where both call a D; IGH ~46-69% (paralogous
      germlines + SHM make IGH D inherently ambiguous).
  - [x] **Double D-D junctions.** For D-D loci (IGH/TRD) a second non-overlapping D
        is sought above a stricter threshold and emitted as `d2_call` + `d2_sequence_*`;
        `np1`/`np2`/`np3` partition the junction between V, the D(s), and J.
        (Limitation: a very long junction can exceed what mmseqs aligns *through*,
        collapsing the projected interior — this only lowers D recall, never
        produces a wrong call.)

- [x] **Multi-node sharding.** `arda split` round-robins a huge FASTA/FASTQ into N
      shards (one pass); `arda merge` concatenates per-shard AIRR TSVs (single
      header); `arda slurm` renders/submits a `submit.sh` chaining split →
      `sbatch --array` annotate → merge via an `afterok` dependency
      (`arda.cluster`). Split/merge/script are unit-tested; the cluster run is
      pending a live SLURM test.

- [ ] **Full AIRR productivity.** `productive` is currently a heuristic (in-frame
      + stop-free V..J span); align it with the complete AIRR productivity rules
      (start codon, stop-codon scan over the whole VDJ, frame of the junction).

- [ ] **Performance.** Optional per-chunk process-pool for inputs where mmseqs is
      not the bottleneck (mostly-non-receptor bulk RNA-seq); mmseqs index reuse.
