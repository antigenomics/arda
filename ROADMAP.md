# arda roadmap

Implemented: offline V·J reference build (5 organisms) with a per-organism
locus-coverage manifest (`loci_manifest.tsv`; warns on empty loci / unreachable
D germlines), MMseqs2 runtime mapping, C++ markup transfer, spec-valid AIRR output
(per-segment V/D/J/C cigars, sequence/germline alignments, read-as-submitted
orientation via `rev_comp`), reverse-complement handling, all-loci single-DB
querying, streaming/bounded-memory FASTQ I/O (optional quality retention),
out-of-frame junction translation, extended V/J-position markup, D-segment mapping
(incl. D-D fusions across all D loci), offline GenBank-vs-IgBLAST test fixtures.

## TODO

- [x] **D-segment mapping.** After V/J transfer, the V..J interior of the junction
      (between the projected `v_sequence_end` and `j_sequence_start`) is aligned
      against the per-organism D germline set by gapless local alignment in the C++
      `_markup.d_local_align` primitive — mmseqs is unreliable on ~8-31 nt D — and
      the best hit is emitted as `d_call` (a comma-separated allele ambiguity list on
      a score tie) + `d_sequence_start`/`d_sequence_end` (AIRR, query coords). D
      germlines ship in `database/vdj/<org>/d_germlines.fasta`
      (VDJ loci only); VJ loci are skipped automatically. Concordance vs IgBLAST:
      TRB/TRD ~97% gene agreement where both call a D; IGH ~46-69% (paralogous
      germlines + SHM make IGH D inherently ambiguous).
  - [x] **Double D-D junctions.** For D-D loci (IGH/TRB/TRD — every locus with a D
        germline set) a second non-overlapping D is sought above a stricter threshold
        and emitted as `d2_call` (also a comma-separated ambiguity list) + `d2_sequence_*`;
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

- [x] **RNA-seq mode (`arda rnaseq`).** Staged pipeline for bulk RNA-seq (~1-5%
      receptor reads). `map`: recall-first, paired-FASTQ (`--r1/--r2`), streams and
      writes only mapped reads (keyed by read id = read-id → junction map) + optional
      candidate FASTA + a JSON report (`arda.rnaseq.map`, reuses `annotate.mapper`
      with a `mapped_only` fast-path). `correct`: CDR3 error correction, a port of
      vdjtools `Corrector` (parent:child count ratio, ≤2 mismatch, 20×) over `seqtree`
      neighbour search (`arda.rnaseq.correct`; optional dep `arda-mapper[rnaseq]`).
      `arda igblast`: all-loci gold-standard AIRR (`refbuild.gold`). Benchmarked in the
      `arda-benchmark` repo vs assembly-based extractors (speed) and IgBLAST (accuracy).
  - [x] **Seamless pipeline integration.** A one-shot `arda rnaseq run` (map + correct
        → `<prefix>.clones.tsv` / `.airr.tsv` / `.arda.json`), a top-level `--version`,
        and a drop-in nf-core-style Nextflow module (`integrations/nextflow/arda/`:
        process + conda env + Dockerfile + per-process config + README). Bulk RNA-seq
        pipelines (nf-core/rnaseq and friends) get per-sample AIRR clonotype tables from
        the same trimmed-FASTQ channel their aligners consume, published to
        `${params.outdir}/arda/`. See `docs/pipeline_integration.rst`.
  - [ ] **Stage 3 — contig assembly.** Reconstruct full-length V(D)J contigs from the
        candidate reads (interface stub in `arda.rnaseq.assemble`; the role
        assembly-based extractors play). Deferred — a de-novo assembler, out of scope
        for the filter-first goal.
    - [x] **Contig cigars (`arda.annotate.contig`).** Once a contig exists, give it
          valid V/J/C cigars + AIRR alignments two ways that produce the *same* record:
          `reannotate_contigs` (re-align the contig through `annotate_records`) and
          `merge_contig`/`merge_contigs` (stitch the reads' existing alignments via the
          C++ `_markup.merge_alignment`, no second mmseqs pass — ~9× faster at ~10^5
          contigs/sample). This annotates a contig; the de-novo assembler above is still
          what would produce one.
  - [ ] **C k-mer prefilter (contingency).** If the MMseqs2 prefilter is
        the throughput bottleneck vs assembly-based extractors, add a parallel spaced-seed germline
        index (new `src/_vjprefilter/` pybind11 ext) that rejects non-receptor reads and
        emits V/J allele hints to prune alignment. Gated on measured need.
