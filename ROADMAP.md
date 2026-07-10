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
      (VDJ loci only); VJ loci are skipped automatically. The call is gated by a
      Karlin–Altschul E-value (`_D_MAX_EVALUE = 0.2`), which replaced four hand-tuned
      per-locus score floors. Concordance vs IgBLAST where both call a D: TRB/TRD ~97%
      gene agreement, IGH 94-98% across the five organisms (recall 52-67%).
  - [x] **Genomic order constrains D×J.** TRBD2 lies 3' of the entire TRBJ1 cluster, so
        deletional joining can never produce a TRBD2–TRBJ1 pair. Unenforced, TRBD2 (16 nt)
        outscored TRBD1 (12 nt) on noise and took 17% of real human TRB J1-cluster D calls.
        `transfer._allowed_d` masks the candidate set (holding the E-value's `n` at the full
        locus size, so a J1 record can only lose an impossible call, never gain a weak one),
        and `scripts/build_d_priors.py` zeroes the same cells before renormalising.
        IGH and TRD place every D 5' of every J: nothing is masked there.
  - [x] **Double D-D junctions.** For D-D loci (IGH/TRB/TRD — every locus with a D
        germline set) a second non-overlapping D is sought above a stricter threshold
        and emitted as `d2_call` (also a comma-separated ambiguity list) + `d2_sequence_*`;
        `np1`/`np2`/`np3` partition the junction between V, the D(s), and J. Recall is
        61% on injected IGH tandems but only 13-15% on TRB, where trimming usually leaves
        one D under the ~7 nt needed to see it; false `d2_call` on true single-D is 0-1%.

- [x] **Multi-node sharding.** `arda split` round-robins a huge FASTA/FASTQ into N
      shards (one pass); `arda merge` concatenates per-shard AIRR TSVs (single
      header); `arda slurm` renders/submits a `submit.sh` chaining split →
      `sbatch --array` annotate → merge via an `afterok` dependency
      (`arda.cluster`). Split/merge/script are unit-tested; the cluster run is
      pending a live SLURM test.

- [ ] **Nucleotide junction re-mapping → recombination-scenario counts.** `arda markup`
      and `arda.cdr3fix` currently take an amino-acid junction (the VDJdb case); add the
      `cdr3nt` path. With nucleotides the germline boundaries are exact rather than
      codon-quantised, so a `(cdr3nt, V, J)` record yields the *full recombination
      scenario*: `delV`, `insVD`, `D` + `delDl`/`delDr`, `insDJ`, `delJ`. Half of this
      already exists — `annotate.dmap.map_d_junction` derives `v_sequence_end` /
      `j_sequence_start` in nt (longest common prefix/suffix vs the shipped
      `cdr3_anchors.tsv` germlines) and calls D with an E-value gate.

      The point is the sufficient statistics: counting scenarios over a repertoire gives
      exactly the E-step tallies needed to re-estimate a generative model by EM in the
      style of Murugan et al. (PNAS 2012) / IGoR — `P(V)`, `P(D,J)`, `P(delV|V)`,
      `P(delDl,delDr|D)`, `P(delJ|J)`, `P(insVD)`, `P(insDJ)` and the insertion Markov
      matrices. arda would then generate its own priors instead of borrowing OLGA's
      (see `arda.dpost`, whose `d_prior.tsv` tables are currently OLGA/vdjrearm-derived).
      Note a scenario is not unique — several `(delV, insVD, ...)` explain one junction —
      so the E-step must sum over scenarios, not take the MAP.

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
