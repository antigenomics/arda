# arda roadmap

Implemented: offline V·J reference build (5 organisms) with a per-organism
locus-coverage manifest (`loci_manifest.tsv`; warns on empty loci / unreachable
D germlines), MMseqs2 runtime mapping, C++ markup transfer, spec-valid AIRR output
(per-segment V/D/J/C cigars, sequence/germline alignments, read-as-submitted
orientation via `rev_comp`), reverse-complement handling, all-loci single-DB
querying, streaming/bounded-memory FASTQ I/O (optional quality retention),
out-of-frame junction translation, extended V/J-position markup, D-segment mapping
(incl. D-D fusions across all D loci), offline GenBank-vs-IgBLAST test fixtures,
run QC (`arda stats`, every mode, with the junction-length / read-length /
clone-size / isotype distributions) and batch QC (`arda qc batch` / `qc report`:
a cohort roll-up with per-batch robust z, and a self-contained HTML dashboard),
single-cell support (`arda cells`: reference-free per-cell
contig assembly, chain pairing, doublet flagging, `umi_count`, and the QC surface),
multi-file samples (read groups) across the CLI, SLURM, Nextflow and Snakemake, and a
generative model of the rearrangement itself (`arda scenarios`, EM over the recombination
scenario set; `arda.hmm`, the same model read as inference).

## TODO

### Next up — ranked, 2026-09-24

Everything below this block is the full backlog, ordered by subsystem rather than by priority.
This is the short list, and each entry says what would make it *done* rather than what it is.

1. **`arda.dpost` cannot consume what `arda scenarios` writes.** `docs/scenarios.rst` calls the
   output a drop-in for `d_prior.tsv`, and it is — by *format*. But `dpost.load_d_prior` is
   `@lru_cache`d on the organism and reads one fixed path (`dpost.py:108-110`), so the only way to
   use a fitted table is to overwrite a file inside the installed database. `arda.hmm.model_for`
   already takes `prior=`; `dpost` does not. Thread a path through `load_d_prior` /
   `posterior_d` and expose it as `arda markup --d-prior PATH`. Small, and it separates *using* an
   estimate from *adopting* one — which is the decision the entry below is about.

2. **11 of the 13 shipped (organism, D-locus) pairs have no `d_prior.tsv` at all.** Not a
   regression: OLGA has no model for them, which is the whole reason the table is derived rather
   than measured. Verified coverage —

   | organism | D loci with germlines | loci with a prior |
   |---|---|---|
   | human | IGH, TRB, TRD | IGH, TRB, TRD |
   | mouse | IGH, TRB, TRD | TRB |
   | rabbit | IGH, TRB, TRD | — |
   | rat | IGH | — |
   | rhesus_monkey | IGH, TRB, TRD | — |

   `load_d_prior` returns `{}` for a missing organism and `posterior_d` then returns `None`
   (`dpost.py:197-199`), so `arda markup` silently has no D posterior for three of five organisms.
   `arda scenarios` fits exactly this table from real junctions, so the blocker is **a cohort per
   (organism, locus)**, not code. Do human and mouse first, where the benchmark repo already has
   the data, and A/B the fitted table against the OLGA-derived one before adopting either.

3. **`--error-rate`'s single default is wrong for variant preservation.** At the default `1e-3`,
   `rnaseq correct` erases both published MIGEC spike-in variants; `1e-5` recovers both exactly,
   and `1e-4` kept both while removing 72 % of real PCR errors on an independent cloud. Not a
   defect — no abundance method separates signal-to-noise ~1, which is why UMI consensus exists —
   but one constant cannot serve both regimes. Wanted: a per-library calibration *rule*, or an
   estimate off the data, not a re-tuned constant. ⚠ Whatever it becomes, it is not a QC threshold
   and must not turn into one.

4. **A per-allele-per-position SHM model**, which two separate entries below are waiting on: the
   HMM cannot be extended to IGH without one (it would explain mutated germline as N-region and do
   *worse* than exact-match anchors, which at least fail safely), and `_map_d`'s amino-acid path
   searches the three translated D frames as independent database entries, tripling `n`, when the
   prior over `insVD` already induces a prior over frame. Biggest item here by some margin.

5. **Genotyping is currently refused on purpose — revisit or close it.** `stats.py` collects
   `allele_candidate` (a recurrent high-quality V mutation carried by >= 50 % of an allele's reads,
   >= 10 reads) and says in its own docstring that it is *"a shortlist to look at, never a call"*.
   Restricting a sample's V reference to the alleles its donor actually carries is annotation, not
   repertoire biology, so it sits on arda's side of the `vdjtools` line — and the evidence is
   already being computed and discarded. But it crosses a stance the repo took deliberately.
   **Author's call; do not cross it silently.**


- [ ] **Single-cell.** Staged plan in **`project/design-singlecell.md`**; that document is
      authoritative and this entry is the index.

  - [x] **S0 — the identifier parser** (`src/arda/cell.py`, 2026-08-14). Lifts a cell barcode out of
        `sequence_id`, which is where every upstream tool already puts it. Dialects `cellranger`,
        `migec`, `prefix`, plus `--cell-regex` and an `auto` sniff. Never: `auto` never considers
        `prefix` and never believes a barcode under 10 nt: a bulk sample named `TCGA` otherwise
        parses as a one-cell library and nothing flags it.
  - [x] **S1 — `cell_id` as an AIRR column** (2026-08-14). `--cell-from` / `--cell-regex` on `map`,
        `amplicon` and `rnaseq`. Never: The hook is in `map.py`'s `flush()`, not `mapper.py:1430` —
        that line is in the unmapped branch, which `mapped_only=True` skips, so it is dead code on
        the only path `map` uses. `CELL_ID` is appended to `extra_cols` **last**.
  - [x] **S2 — `locus` and `clone_row` in `--read-map`** (2026-08-14). Never: `junction` alone is not
        the clonotype key — `correct` keys on `(locus, v_call, j_call, junction)` — so the two-column
        form did not close the per-cell join it existed for.
  - [x] **S3 — per-cell assembly, chain pairing and the QC surface** (`arda cells`, 2026-08-15).
        It went further than planned and the plan was wrong about where the junction comes from.
        Not the `--read-map` join and not the per-read junction: `arda cells` **assembles each
        cell's contigs first**, reference-free, from one UMI consensus per molecule, and
        annotates the contig. Reads of one (CB, UMI) are co-terminal so a molecule covers one
        window of the transcript, but different molecules start at different positions, so a
        cell's molecules TILE it -- 99.90% of the k-mers of Cell Ranger's 943 contigs are already
        in their own cell's molecules before any assembly runs. Measured against Cell Ranger on
        `sc5p_v2_hs_PBMC_1k`: **933/943 CDR3s recovered verbatim (0.9894), TRB 479/479**, chain
        recall 0.9777, contig N50 536 nt, 23 s for 479 cells.
        **Never: Phase a component before consensing it or a doublet is invisible by construction** --
        two chains of one locus share their constant region, so the layout puts them in one
        component and the column consensus averages their junctions into a third sequence. Worth
        0.9714 -> 0.9777 on real data, and on a synthetic doublet the difference between one
        918 nt contig with no callable junction and both true junctions.
        **Never: The extra-chain gate is PRODUCTIVITY first, count second.** Of the extra chains Cell
        Ranger agrees with 60/60 are productive; of those it does not, 20/128 are.
        Never: `doublet_candidate` on the heavy slot only; a second light chain is allelic inclusion.
        Diagnostics: knee, doublet scatter, chain support, contig lengths, filter sweep
        (`--reference`), and clustering agreement via `arda.partition`. `docs/singlecell.rst`,
        `notebooks/singlecell_qc.py`.
        **Never: The knee is Kneedle at its GLOBAL maximum, and it is guarded.** Kneedle's published
        local-maxima walk is degenerate without the paper's smoothing spline -- 287 local maxima
        on a real curve, stopping at rank 15 of 136,032. And a global maximum always exists, so
        an ambient-only library returns a rank too; `find_knee` refuses one below 10x the mean
        molecules per barcode. Agrees with migec's C++ exactly (rank 376, 308 molecules).
  - [x] **S4 — `umi_count`.** Shipped 2026-09-24 as `correct --cell-from` / `--cell-regex`.
        `duplicate_count` and `consensus_count` are untouched — they are AIRR-spec fields.
        Never: it was NOT blocked on a migec format decision, and the entry that said so was
        answering the wrong question. AIRR asks for **distinct UMIs**, and `.c<k>` / `.<m>`
        subdivide reads that ALREADY SHARE one `<umi>`, so no reading of the conditional suffixes
        can change the count; `<sample>`, `<cell>` and `<umi>` are unconditional. What is still
        blocked on migec is *molecule identity* for a `<sample>.mig.tsv` join, which nothing in
        arda needs.
        Never: measured, not reasoned — `umi_count < consensus_count` on a **saturated barcode**
        (migec splits one barcode into two molecules; differences outside the junction leave both
        records in one clonotype: 3 / 3 / 2), and **not** in contig mode, where components never
        overlap so at most one contig carries a junction and the rest are dropped by
        `complete_only`. The first draft of the rescope claimed the opposite.
        Never: OMITTED, never 0 or 1, when the dialect names no UMI or the ids do not parse; and
        `--cell-from` is a BOTH-HALVES flag in `cluster.regime_flags`. Details in
        `project/design-singlecell.md` S4.
  - [ ] **`arda singlecell` stays reserved as a MODE name** — the work lives in `arda cells`. Never: Its only sensible preset is the all-False vector,
        which is byte-for-byte what `--exact` already gives on either existing mode — a no-op mode.
        It ships when a measured speed row differs from both presets.

- [ ] **TRUST4 head-to-head on AMPLICON at full depth, plus an IgBLAST-truth accuracy leg.**
      Scheduled. What exists today is wall clock only, at 100 k / 500 k reads, same job and same
      staged input (round 20): IGH_repertoire **201.05 s** vs TRUST4 615.04, IGH_naive **136.51** vs
      359.66, migec_exp1_TCR **316.25** vs 423.96, migec_exp1_IGH 223.77 vs **225.10**. Missing: an
      hours-scale full-depth run, and **any** amplicon accuracy figure for TRUST4 — arda's and
      MiXCR's amplicon accuracy is measured against IgBLAST, TRUST4's is not. Never: Do not project the
      hours from the per-100 k walls: a projected ratio quoted as measured is the single mistake
      this project has made most often. The arm is `cluster/ampacc7.sbatch` in arda-benchmark
      (IgBLAST on 10,000 pairs at stride 100, one merged single-end FASTQ so every tool sees the ids
      arda emits); the TRUST4 leg needs writing.


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

- [x] **Samples split across files (read groups).** `--r1`/`--r2`/`--id` are repeatable and
      position-matched, and `--samples sheet.tsv|csv` reads the nf-core `sample, fastq_1, fastq_2`
      spelling with repeated ids merging in row order. `arda.samples` is the only new concept.
      Never: a multi-file sample is an ALREADY-SHARDED sample -- `pipeline.run` maps each read
      group and concatenates in DECLARED order before the existing `finish`, so nothing below
      Stage 1 changed. Byte-identical to the same reads in one file, verified on
      `tests/data/rnaseq_real` cut into four.
  - [x] **The read group is the unit of parallel work, the sample is not.** `arda cluster plan`
        emits `readgroups.tsv` + `samples.tsv` for any scheduler with the commands FILLED IN;
        `arda cluster submit-samples` renders them as two SLURM arrays; Snakemake
        (`integrations/snakemake/arda/`) and Nextflow schedule the same way. 5 samples x 4 lanes
        is 20 jobs, not 5. Never: samples run SEQUENTIALLY in one CLI call -- mmseqs threads
        internally, so N samples at cores/N is slower than N in a row.
  - [ ] **`--jobs N` for concurrent samples on one box** is deliberately absent. It can only pay
        when each sample is too small to saturate mmseqs; measure a many-small-samples sheet
        before adding it.

- [x] **Batch quality control.** Staged plan in **`project/design-qc.md`**; that document is
      authoritative and this entry is the index. The per-sample QC table has existed since
      2.14.0 and says in its own docstring that its long format is there to be joined across
      samples -- nothing performed that join, no mode but bulk wrote the table, and nothing
      rendered it.
  - [x] **The distributions, and one address per fact.** Five sparse scopes keyed
        `locus:bucket` (`junction_aa_len`, `read_len`, `clone_size`, `isotype`,
        `chain_support`), all from columns `stats.py` was already reading and discarding.
        Never: min/max/mean cannot show a SHAPE -- a bimodal junction length is two primer sets
        in one tube, a read-length cliff is an adapter left on, and each reads as an ordinary
        mean. Three run-scope defects fixed with them, all of which corrupt a cross-sample join
        rather than failing it: `_merge_map_reports` renamed `wall_seconds`/`peak_rss_mb` for
        lane-split samples only, `segment_search.reasons` never survived the one-level flatten,
        and an empty aggregation was written blank instead of omitted.
  - [x] **`arda cells` writes the same table.** Same scopes as bulk, so a cohort can mix them;
        `pairing_rate` / `doublet_rate` derived, since they are the AIRR chain-pairing QC and
        `cell_summary` had already assigned every cell its status. Never: a chain the
        extra-chain gate marked `extra` is not yield -- it is ambient 96-97 % of the time.
  - [x] **`arda qc batch` -- the roll-up.** Reads ONLY the per-sample `*.stats.tsv`, never an
        AIRR or a clonotype table, so a 1,000-sample cohort is one `concat` on a laptop against
        results copied off a cluster. The sheet gains optional `project`/`batch` labels.
        Never: no shipped threshold. Each metric carries its group's median, MAD and robust z
        (|z| >= 3.5); a group under 5 samples gets no z, because MAD over four points is not a
        scale. Flags, never filters, as `stats` already refused to decide.
  - [x] **`arda qc report` -- one self-contained HTML file.** Data inlined, SVG drawn by plain
        JavaScript, **no CDN, no bundle, no new dependency** -- it opens air-gapped and survives
        the results directory. Same stance as `scplot`'s always-written gnuplot script.
  - [x] **Why a read did NOT map.** `mapped_fraction` alone made "wrong organism", "no
        receptor content in this library" and "`--min-score` too strict" one number. The
        `run/map` scope now carries an `unmapped.*` ledger — `prefilter_rejected`, `no_hit`,
        `hit_not_in_reference`, `constant_only`, `below_min_score` — plus `accounted`, and the
        `sample` scope carries each as a fraction of `total_reads`. Never: it is a LEDGER, so it
        states its own completeness instead of leaving the reader to sum, and `accounted` is
        RECOMPUTED on a shard merge rather than summed — a bucket wired into the single-file
        path and not the merge would look fine on a laptop and be wrong on every cluster run.
        Never: `constant_only` counts READS, not the fragments `constant_only_fragments` counts;
        the rule also drops the constant-only MATE of a fragment it keeps, and using the
        fragment count left 45 of 1,320 reads on `tests/data/rnaseq_real` in no bucket at all.
    - [ ] **Splitting `no_hit` further** — no V / no J / V and J on different mates, as MiXCR
          does — is not done. It needs the reason a hit was rejected to survive out of
          `_best_hits`, which is inside the per-read hot path; measure the cost before adding
          it, and only if a real library makes `no_hit` ambiguous.
  - [ ] **Repertoire biology stays out**, deliberately. Diversity, clonality, rarefaction,
        overlap and cross-sample clonotype matching are `vdjtools`', which takes arda as a base
        dependency. arda's QC answers "did this run work, and is this sample like its batch".
  - [ ] **The MiAIRR Repertoire/DataProcessing document** is absent because arda cannot fill
        the read-QC fields honestly -- `total_reads_passing_qc_filter` is the facility's count,
        upstream of annotation, and `mapped_reads` is not it. Unblocks when someone submits arda
        output to iReceptor or VDJServer.

- [x] **Multi-node sharding.** `arda cluster split-fasta` round-robins a huge FASTA/FASTQ into N
      shards (one pass); `arda cluster merge` concatenates per-shard AIRR TSVs (single
      header); `arda slurm` renders/submits a `submit.sh` chaining split →
      `sbatch --array` annotate → merge via an `afterok` dependency
      (`arda.cluster`). Split/merge/script are unit-tested; the cluster run is
      pending a live SLURM test.

- [x] **Nucleotide junction re-mapping → recombination-scenario counts.** Shipped 2026-09-24 as
      `arda scenarios` (`src/arda/scenarios.py`, `project/design-scenarios.md`, `docs/scenarios.rst`).
      EM over the scenario set, writing the same long `locus/kind/key/value` table `d_prior.tsv`
      uses — a drop-in, so `arda.dpost` can stop marginalising OLGA's numbers.
      Never: the counts are EXPECTED counts summed over scenarios, not one MAP reading — 4,346
      tuples reproduce one real human TRB junction exactly, and counting one of them biases every
      distribution toward less trimming and shorter inserts.
      Never: an insertion costs its own SEQUENCE (`0.25^len`), not just its length. Without it EM
      is degenerate — measured, `insVD` mass walked to 10–11 nt on 503 real TRB junctions with the
      log-likelihood rising the whole way. With it, `insVD` peaks at 4 nt and `dlen` for
      `TRBD1*01` at 4–5, reproducing `dpost`'s independently measured median of 5.
      Open: **adopting** an estimate as the shipped `d_prior.tsv` is a measurement and a release
      decision, deliberately not a side effect of this.

  - [ ] **What is still open: the `cdr3nt` path through `arda markup` / `arda.cdr3fix`.**
        `arda scenarios` reads nucleotide junctions and `annotate.dmap.map_d_junction` maps a D
        into one, so the *model* half of this entry is done. What is not is bare-record REPAIR on
        nucleotides: `cdr3fix` is 715 lines of amino-acid alignment (`templated_aa`, the
        anchor-distance rule, the `Failed*` ladder) and `arda markup` exposes only that. A VDJdb
        record with `cdr3nt` still has to be translated to be repaired, which throws away exactly
        the resolution — germline boundaries exact rather than codon-quantised — that made the
        nucleotide path worth having. Separable, and nothing in the pipeline is blocked on it.

- [x] **`arda.hmm` — a hidden semi-Markov model of V→N1→D→N2→J.** Shipped 2026-09-24
      (`src/arda/hmm.py`). As this entry predicted, it was not a second project: the E-step of
      `arda scenarios` **is** the forward-backward pass, so `scenarios.lattice` is the one
      implementation and `hmm` is it read as inference — `log_likelihood` (marginal, not Viterbi),
      `posterior_d`, `best_scenario`. Semi-Markov because deletions and insertion lengths have
      tabulated non-geometric durations; conditioned on the mmseqs V/J call, which keeps the state
      space to `(delV, insVD, D, delDl, delDr, insDJ, delJ)` and the cost per clonotype.
      `model_for(prior=...)` scores against a table `arda scenarios` wrote, so a model fitted on
      one cohort can score another.
      Never: it does NOT replace `arda.dpost` — that answers the amino-acid question, where there
      are no nucleotides to run this on — and it **gates nothing**. Nothing in the annotation path
      calls it, because of the two measured negatives below.
      Open: `_allowed_d` (TRBD2×TRBJ1 = 0) is still enforced in three places rather than as one
      transition zero; folding it in is a refactor with no measured payoff, so it waits for one.

      **Two measured negatives, so nobody re-runs them.** (i) Re-ranking nt D candidates by
      `λ·S + log P(insVD) + log P(dlen) + log P(insDJ) + log P(D|J)` changes *nothing*:
      gene accuracy 98.9→97.8 % (IGH), 94.2→94.5 % (huTRB), flat elsewhere, identical call
      rate and `d_start` error. With 10–18 matched nt, `λ·S` is 11–20 nats and the prior
      moves ±3. (ii) Replacing the E-value gate with a present/absent Bayes factor buys
      IGH ~+3 pp recall at matched FP (93.6 % vs 90.7 % @ ~2 % FP) and nothing for TRD —
      but needs a *per-locus* threshold (BF>6 for IGH, BF>10 for TRD at the same FP),
      reintroducing the four knobs the E-value removed, and only exists for 4 of the 15
      (organism, D-locus) pairs arda ships. So: the HMM's value is posteriors, uncertainty,
      one home for the constraint, and the EM E-step — **not** D-call accuracy.

      Prerequisite for IGH: a per-allele-per-position SHM model. Without one the HMM will
      explain mutated germline as N-region and do *worse* than the current exact-match
      anchors, which at least fail safely (SHM truncates the match, widening the interior,
      never clipping the D). Also fix `_map_d`'s aa path while there: it searches the three
      translated D frames as independent database entries, tripling `n`, when the prior over
      `insVD` induces a prior over frame (as `dpost` already knows).

- [x] **`cdr3fix` repairs to canonical, always.** `cdr3_repaired` is accepted only when it
      opens with Cys104 and closes with Phe/Trp118 (`_canonicalise`); otherwise the submission
      comes back untouched. `good` implies canonical by rule. Trimming flanking framework
      (`_MAX_TRIM = 3`) is budgeted separately from inventing germline (`_MAX_FIX = 2`) —
      removing residues the germline never explained is a much smaller risk than adding ones
      never observed. Each trimmed residue costs `_TRIM`: a *free* trim could tie the untrimmed
      alignment, win the tie-break, and eat the conserved Phe of clean short IGK junctions.
      Result on the 250-row fixture: VDJdb's repair reproduced 100/100 (was 98), zero novel
      rewrites, and `good`/`vCanonical`/`jCanonical` agree with VDJdb on every row. OLGA false
      repairs unchanged at 0.35 %. `family*01` was dropped from the allele ladder: dead for all
      five shipped organisms. `FailedReplace` turned out to be reachable after all (three
      substitutions beside a long J anchor, `max_replace >= 3`), and is now tested.

- [x] **Productivity: FR4 is scanned for stops.** `productive` and `stop_codon` covered the
      V-side regions and the junction, and the junction ends AT [FW]118 -- FR4's first residue --
      so residues 2..n of the J were looked at by neither. Fixed; `tests/synthetic/
      test_fr4_stop.py` exercises it. Never: measured before changing. On the committed IgBLAST
      fixtures every evaluable row carries an FR4 (mean 10.1 aa human, 8.6 mouse) and **none of
      4,217** contains a stop, so the fix left all 4,843 annotated rows byte-identical -- a real
      hole, unexercised by real data, which is why it needed a test written on purpose.
  - [ ] **The rest of "full AIRR productivity" is aimed at the wrong thing, measured 2026-09-24.**
        The entry used to ask for the start codon and a whole-VDJ stop scan. The stop scan is
        done (above). A start codon is not assessable on a read fragment and IgBLAST does not
        require one either. And arda-vs-IgBLAST productivity was measured on the committed
        fixtures: of 3,601 rows where both tools were evaluable AND called the same
        `junction_aa`, they disagree on **8** (0.22 %) -- 3 human, 5 mouse. In every one arda
        says `T`, IgBLAST says `F`, and **both** report `stop_codon = F`. So not one of them
        would be changed by any productivity rule; all 8 are IgBLAST reporting
        `vj_in_frame = F` where arda reports `T`.
  - [x] **`vj_in_frame`: the 8 disagreeing rows were re-tested against the C gene, and arda is
        right on all of them** (2026-09-24). The first pass asked the wrong question -- it looked
        for the called J's `templated_aa` in each frame, which tests whether the junction's 3'
        residues *resemble germline*, not whether the rearrangement is in frame. The author
        supplied the correct test: `vj_in_frame` means the constant exon spliced onto that J
        translates, so **splice C on and count stop codons**.

        Done, on every disputed read. `JX277388.1` (TRBJ2-5*01) + `TRBC2*01`: 0 stops over 375
        nt, yielding the real mouse TRBC2 protein (`...FGPGTRLLVL|EDLRNVTPPKVSLFEPSKAEIANK...`).
        `JX277345.1` (TRBJ1-3*01) + `TRBC2*01`: 0 stops. `MH918759.1` (TRAJ31*01) needed no
        reconstruction -- the entry carries 66 nt of TRAC and translates clean.

        What these reads actually have is an **indel inside the J**, after the anchor codon
        (`JX277388.1` carries 2 nt `TRBJ2-5*01` does not, `JX277345.1` carries 1). So the J's 5'
        germline residues and its own FR4 cannot both be in the germline frame. arda lands on
        the FR4/C side, which is the side that decides function; a J-frame lookup keyed on the
        J's 5' alignment start lands on the other one. That is the whole disagreement.

        `X60894.1` (TRAJ9*01) is **unanswerable, not wrong**: the entry is 45 nt and stops 3 nt
        into FR4, so it carries one FR4 residue. Completing the J from germline and splicing
        `TRAC*01` gives 0 stops over 261 nt, but one residue is not evidence.

        **No `partial_nt` change is wanted.** `phase = (fwr4_start - v_coding_start) % 3` is
        already the correct test. The open decision is closed -- do not reopen it on
        `templated_aa` evidence.

    - [ ] **Never: do not measure this by comparing the junction's tail to `templated_aa`.**
          Tried, and it is confounded by SHM, not by anchoring: it fires on 51.8 % of IGL,
          51.4 % of IGK and 39.8 % of IGH (mean `v_identity` .955-.968) against 3.4 % of TRB and
          7.1 % of TRA (mean `v_identity` .992). It measures hypermutation reaching the J, which
          is biology -- and it is what produced the wrong 5:3 verdict above. The only sound test
          is translating into the constant exon. Documented for users in `docs/productivity.rst`.

- [ ] **Performance.** Optional per-chunk process-pool for inputs where mmseqs is
      not the bottleneck (mostly-non-receptor bulk RNA-seq); mmseqs index reuse.

- [x] **RNA-seq mode (`arda rnaseq`).** Staged pipeline for bulk RNA-seq (~1-5%
      receptor reads). `map`: recall-first, paired-FASTQ (`--r1/--r2`), streams and
      writes only mapped reads (keyed by read id = read-id → junction map) + optional
      candidate FASTA + a JSON report (`arda.rnaseq.map`, reuses `annotate.mapper`
      with a `mapped_only` fast-path). `correct`: CDR3 error correction, a port of
      vdjtools `Corrector` (parent:child count ratio, ≤2 mismatch, 20×) over `seqtree`
      neighbour search (`arda.rnaseq.correct`; core dep `seqtree`).
      `arda igblast`: all-loci gold-standard AIRR (`refbuild.gold`). Benchmarked in the
      `arda-benchmark` repo vs assembly-based extractors (speed) and IgBLAST (accuracy).
  - [x] **Seamless pipeline integration.** A one-shot mode command (map + assemble + correct
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
        index (new `src/_vjprefilter/` nanobind ext) that rejects non-receptor reads and
        emits V/J allele hints to prune alignment. Gated on measured need.
