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

### Next up — ranked, 2026-09-25

Everything below this block is the full backlog, ordered by subsystem rather than by priority.
This is the short list, and each entry says what would make it *done* rather than what it is.

Items 1 and 2 shipped on 2026-09-25 and are kept here with what they actually bought; item 3 is
what is left of the junction-recall gap. Evidence and method:
`~/vcs/projects/2026-arda-benchmark/results/round27` and `results/round28`.

1. ✅ **The Cys104 gate is loosened, and its cost was 92 reads — not the 512 round 27 claimed.**
   `transfer.py`'s `v_anchor_ok` now admits a junction whose exact prefix fails if **6 of its
   first 8 bases** match some called V's own `germline_nt`. ⛔ **Do not requote round 27's "512
   reads, 44.6 % of the gap."** Those are reads arda refuses *and MiXCR gets right*, but arda's
   own ungated junction is wrong on 420 of them, so they were never reachable from the gate; the
   gate refuses 1,452 junctions of which **92 are correct and 1,360 are genuine 5' over-extensions**,
   exactly as `tests/unit/test_junction_and_j_evidence_gates.py` has said since round 18. Measured
   on both sides and on a **held-out locus** (the TRB amplicon `SRR5233641`, never looked at while
   the rule family was written): **+144 correct junctions and ONE extra over-extension across
   92,466 truth junctions**. Junction recall TRA .9473 → **.9481**, TRB .9898 → **.9922** against
   MiXCR's .9944 — the held-out gap halved — with every other metric unchanged to four places.
   ⚠ **On TRB the exact gate was a net loss**: 116 correct junctions discarded to catch 21 wrong
   ones, for 0.00013 of precision. There is nothing further here: arda's *ungated* ceiling on the
   TRA amplicon is .9493, still short of MiXCR's .9708 (item 3 and the note below it).

2. ✅ **`build_germline_dbs` honours `locus.v_shared`, and the shipped TRD exposure measured
   zero.** The IgBLAST V database is now built from the same allele set the scaffolds are, deduped
   by seq id (`makeblastdb` dies on `Duplicate seq_ids ... LCL|TRAV14/DV4*01` otherwise). On the
   chimera-enabled TRA locus that is **7 of 483 scaffolds → 49 of 483** with every V call correct.
   ⚠ **`Locus("TRD", …, v_shared=("TRAV", "/DV"))` was never actually exposed** — the live TRD
   reference keeps **73/110** complete markup either way, because IMGT files the dual-use
   `TRAV*/DV*` genes under **both** stems and the TRDV germline file already carried all 15 of
   them. The fix makes the database match the scaffold allele set by construction rather than by
   IMGT filing accident; it moves no shipped locus.

3. **TRDV × TRAJ is now the whole remaining named junction gap, and `--allow-chimeras` is the
   wrong gate for it.**
   680 truth reads (1.4534 % of this library) carry a TRDV and **all 680 pair it with a TRAJ, not
   one with a TRDJ**; arda gets .100 of them, MiXCR .968, IgBLAST calls them at `v_score ≥ 70`.
   Recovering the class is **.9473 → .9604**. The flag exists because the pairing is a domain
   judgement (`refbuild/loci.py:113-141`) and that has not changed — but the cost of the current
   default is now priced, and the ceiling that made turning it on nearly a no-op (+7 scaffolds) is
   item 2 plus one open IgBLAST question: with the V db fixed, **434 of 483 still get no `j_call`**,
   and the split is purely J-driven — exactly **7 of 69 TRAJ alleles work, with all 7 TRDV alleles**
   (`TRAJ13*01/02`, `TRAJ16*02`, `TRAJ24*01/02/03`, `TRAJ39*01`). Ruled out: germline length (47 of
   the 62 failing are at least as long as the shortest working one), scaffold geometry (`n_pad`
   0/1/2 exactly as for pure TRA), and anchor availability (all 7 TRDV carry a functional
   `anchor_nt`, all 73 TRAJ anchors present). **Ask before moving the default either way.**

   ⛔ **What is left after item 3 is not addressable from the V side.** The 1,360 over-extensions
   the gate correctly catches are **461 distinct junctions**, one TRAV25\*02 sequence accounting
   for 593 reads (42.5 % of every wrong junction on the library) and 648 of them over-extended by
   exactly 9 nt — the scaffold's `V + 9 nt N-pad + J` bridge. Re-finding the true start by sliding
   the V germline along the junction is measured and dead: at the **true** offset the called V's
   `germline_nt` matches **0 bases on 1,040 of 1,369**, because the V really was chewed back past
   Cys104. Every slide rule tried is a loss on both libraries (TRA +10 right/+289 wrong, TRB
   +0/+26). Recovering that class needs the CDR3 start located without V-germline evidence, which
   is a mechanism arda does not have today.

4. ⛔ **ANSWERED on IGH, and the answer is the opposite of the TCR result** (round 31,
   2026-09-25). On a real human IGH multiplex V-primer amplicon (ngsik `BCR_Multiplex`, 100,000
   pairs, 251 nt, 3 reps, medians, idle machine) Stage-3 assembly is not rescuing a handful of
   clonotypes — **it is producing all of them**. Neither mate spans V into J on its own: R1
   stops at the V's 3' end and R2 never reaches the J, so the per-fragment AIRR carries a
   `junction` on **170 of 98,282 rows (0.2 %)** and the assembled contigs carry one on
   **24,655 of 24,655 (100 %)**. Measured cost and yield:

   | arm | wall (s) | CPU (s) | clonotypes | AIRR rows |
   |---|---:|---:|---:|---:|
   | `arda amplicon` | 255.64 | 1,635.91 | **13,040** | 196,502 |
   | `--no-assemble` | **169.07** | **1,212.22** | 108 | 196,502 |

   Assembly costs **1.51x wall and 1.35x CPU** and buys **12,932 of 13,040 clonotypes
   (99.2 %)** — against 0.025 % on a TRA amplicon and 0.013 % on TRB. ⛔ **`--no-assemble` must never become an amplicon preset.**
   On TCR it costs a rounding error; on IGH it costs the library. The original entry, which is
   what it is answering:

   **Amplicon Stage-3 assembly costs 29 % of wall to rescue 12 reads.** `arda amplicon` against
   `--no-assemble` on a 100 k TRA amplicon, 3 reps, medians: **13.46 → 9.38 s wall (1.43×)**,
   30.05 → 22.93 s CPU, for **5 clonotypes of 19,841 (0.025 %)** and +3 reads. The reason is
   structural — the mode exists for reads that span V into J, so **89.4 % of its mapped reads
   already carry a complete junction** and assembly has nothing to build; on bulk that figure is
   10.8 % and the same stage rescues 2,584 reads from 1,931 complete contigs, so `rnaseq` must keep
   it on. Round 28 **replicated it on a second library** — the TRB amplicon `SRR5233641`:
   **14.21 → 11.18 s wall (1.27×)**, 30.14 → 24.49 s CPU, for 3 clonotypes of 22,589 (0.013 %) and
   **16 contigs built across 100,000 reads**. ⚠ **Still not shippable**: both are TCR amplicons
   from one patient, so that is evidence the TRA result was not a fluke and no evidence about IGH.
   IGH RepSeq has long hypermutated CDR3s and is where amplicon assembly would pay, and it is on
   aldan3. Done means that A/B runs first; until then `--no-assemble` is a documented option, not
   a preset change.

5. ✅ **`--adaptive` is re-priced and the two stale read-path claims are gone (2026-09-25).**
   Help text, the comment at `_ADAPTIVE_TRIGGER`, `rnaseq/map.py`'s dnaio docstring and its
   `--chunk-size` sweep all now carry the round-27 numbers instead of the fixture's. Nothing
   about the default changed -- `--adaptive` stays off -- and the reason is now in the help
   text a user actually reads: the junction movement is **one-directional** (89 of 93 rows go
   to EMPTY) and **91 of 93 sit above the trigger**, so it is not a tuning problem. What
   follows is the original entry, kept for the measurement:

   Two
   documentation defects of the class `CLAUDE.md` doc-invariant 1 is about — a number that was true
   and stopped being true. (a) `--adaptive` is **1.84× on bulk map wall and 3.04× on CPU**
   (16.34 → 8.90 s, 198.20 → 65.13 s, 1,033 → 771 MB on 660 k pairs) with the read set preserved
   exactly, but it moves `junction` on **93 of 35,795 rows and the movement is one-directional —
   89 to empty, 4 the other way**, a net −85 against 3,856 reads carrying one. 91 of the 93 sit at
   90–150 bits, so `_ADAPTIVE_TRIGGER = 90` does not catch them. The default stays off; the help
   text should quote this rather than "3 of 453 reads" from a fixture. (b) `rnaseq/map.py:479`
   claims reading is "65 % of a bulk run" — measured now at **0.73 s of map's 16.19 s = 4.5 %**,
   because the dnaio port that comment motivated made its own premise false. `map.py:333`'s chunk
   sweep is stale for the same reason: 1 chunk vs 4 is **16.48 vs 16.34 s**, byte-identical output.

6. **IG accuracy: the gap is READ COVERAGE, the one instrument left is a wider V tie list.**
   ⛔ **Round 29's "IGH `v_gene` .9004, seven times the junction gap" is WITHDRAWN**, and so is
   the somatic-hypermutation suspicion it raised. Measured on two bulk libraries arda did not
   select (benchmark round 30): stratified by IgBLAST's own `v_identity` the deficit is
   **non-monotonic** — `>= 99 %` (unmutated) .9425, **`97–99 %` .7518 (worst)**, `92–95 %` .9697
   (best). Stratified by how much V germline the read covers it is monotonic and steep:
   **`< 60 nt` .1170, `>= 200 nt` .9896** — and .9896 *is* the TRA amplicon's .9867. 5.9 % of
   reads produce ~56 % of every V miss. The top confusion `IGHV1-69 → IGHV1-18` is **4,745 of
   ~9,080 misses**, and on those reads IgBLAST aligns germline **242–296 = 55 nt of the V's 3'
   end**: a 5'RACE read runs C → J → V, so the bases separating one IGHV from another are not in
   the read, and IgBLAST — which lists 1.64 V genes per read itself — is making a different
   tie-break on the same missing evidence. `skills/arda/references/genotype.md` already measured
   the threshold on TCR: TRAV ~175 nt from the 3' end, TRBV ~150.
   ✅ **Corrected: arda's IG V-gene accuracy is .93–.98 on bulk and .9896 when the read carries
   >= 200 nt of V.** Replicated on `SRR5233639`/`SRR5233640`, IGH/IGK/IGL, the two samples within
   0.5 points.
   ⛔ **A V-call evidence gate was priced and REJECTED — do not re-propose it without new
   evidence.** Dropping a call that starts at germline `>= 150` and spans `< 60` nt is a 7.6 : 1
   win on IGH 5'RACE (precision .92585 → .97693 for −0.69 pt recall) and a clear loss on bulk
   IGH / IGK / IGL (−3.24 / −6.22 / −2.39 points of recall for essentially no precision). A
   shipped constant would be wrong more often than right, which is why this repo ships none.
   ✅ **The widening instrument is SWEPT and the answer is per-LOCUS (round 32, 2026-09-25).**
   arda already ships the widening — `arda resolve-ties`, off by default — so the sweep needed no
   new mechanism, just a scorer that prices the cost. Eleven arms, five geometries, three loci,
   two receptors, against an `arda igblast` truth. **Exact `v_gene`-SET agreement**:

   | arm | locus | truth genes/read | exact before | exact after | Δ |
   |---|---|---:|---:|---:|---:|
   | `SRR5233639` / `SRR5233640` bulk | **IGK** | 1.60 / 1.62 | .5707 / .5558 | **.9373 / .9308** | **+36.7 / +37.5** |
   | `SRR5233639` / `SRR5233640` bulk | **IGL** | 1.09 | .8586 / .8577 | **.9137 / .9115** | **+5.5 / +5.4** |
   | `SRR5233641` amplicon | **TRB** | 1.08 | .9299 | **.9694** | **+4.0** |
   | six IGH arms (bulk, 5'RACE x3, multiplex) | IGH | 1.12–1.64 | .8255–.9755 | .6860–.9608 | **−0.6 to −14.9** |

   ⚠ **Intersection recall rises on ALL ELEVEN arms**, IGH included. Reading recall without the
   exact-set number would have shipped a 10-point IGH regression as a 5-point IGH win; that is
   what `scripts/vcall_widening_cost.py` exists for.
   ✅ **The mechanism is predictive, not empirical.** The deciding quantity is the ambiguity
   DEFICIT — IgBLAST's genes/read minus arda's, before widening. **arda's IGK call is 0.42
   genes/read too NARROW** (1.18 against 1.60), and closing that gap is the whole 37-point win.
   **arda's IGH call is already as wide as IgBLAST's** (deficit 0.00–0.07 on every arm, 1.589
   against 1.64 on 5'RACE), so there is nothing to recover and each arm pays the overshoot — up to
   +0.40 genes/read against a 0.06 deficit on bulk IGH. A user can run that check without a truth
   file: compare `genes/read` before and after.
   ✅ **Shipped: `arda resolve-ties --loci IGK,IGL`.** Verified on one mixed bulk library, one
   command: IGH byte-identical to untouched (exact .8255, 108 clonotypes), IGK and IGL taking the
   full win. ⚠ The cost is `solo` — the share of IGK reads naming exactly one gene falls
   .8195 → .3883 — which is the honest statement of what an IGK read supports and why this stays
   an opt-in command, not a default. Clonotype cost: **+2 on IGK, +3 of 23,559 on TRB, ZERO on
   IGL**. Evidence: `results/round32`.
   ⛔ **Do not re-propose a bit-score MARGIN for this.** The instrument that was wanted is
   measured; what it recovers is an IGK-specific narrowness, and a margin would be a second
   mechanism for the same job.
   ✅ **Shipped 2026-09-25: the coverage floor is documented**, in `docs/usage.rst`
   ("What an IG V call means when the read is short") — both span tables, why position beats
   length, why SHM does not order it, and why no threshold ships. It is the one IG finding a user
   can act on today, and `v_sequence_start`/`v_germline_start` are already AIRR columns arda
   writes, so they can filter on the span themselves.
   ✅ **The second-patient check is finished (SPX8151, `SRR5233637`/`SRR5233638`).** It is a second
   PATIENT, not a second geometry, and it replicates the only thing it was asked to: `v_gene` recall
   is **monotonic in V germline span on all six arms** (two samples x IGH/IGK/IGL) —

   | locus | SRR5233637 | SRR5233638 | < 60 nt | 60–90 | 90–120 |
   |---|---:|---:|---:|---:|---:|
   | IGH | .9103 | .9058 | .6887 / .6763 | .8238 / .8538 | **.9649 / .9546** |
   | IGK | .9770 | .9782 | .9137 / .9013 | .9684 / .9804 | **.9890 / .9890** |
   | IGL | .8512 | .8793 | .4286 / .4490 | .5568 / .6449 | **.9588 / .9822** |

   1,317–1,558 truth reads per locus per sample at `v_score >= 70`. The two samples agree within
   0.5 points on IGH and IGK; IGL differs by 2.8 points on 605 reads each.
   ✅ **And it names where the tie list CANNOT be the fix.** `truth genes/read` — how many distinct
   V genes IgBLAST itself listed — is **1.49–1.71 on IGK** (the judge could not separate them
   either, so a wider list is exactly right) but **1.02–1.09 on IGL**, where IgBLAST resolves a
   single gene and arda still misses 55 % of the `< 60 nt` bin. Those are tie-BREAK errors on
   evidence that exists, not missing evidence, and they are a separate instrument from item 6's.

   ⚠ Also open and **the author's call, 476 reads**: `IGHV3-52` and `IGHV3-71` have **0 scaffolds**
   in arda because `load_functional_alleles` excludes IMGT ORF/pseudogenes by design while
   IgBLAST's DB includes them — **82 of its 121 human IGHV genes**. `CLAUDE.md` ("Reference
   vocabulary — three checks before a gene joins or leaves") is now the rule for deciding this,
   and `SOURCES.md` carries the evidence already gathered: **IGHV3-71 is a pseudogene
   (`ENSG00000254056`) that GTEx nonetheless shows transcribed wherever B cells are** (spleen
   3.45 TPM, ileum 1.30, EBV-lymphocytes 0.378, ~0 in ~35 of 54 tissues), and it sits at
   **0.9172** 3' identity to `IGHV3-49` — **separable, so those reads are being mis-assigned
   today, not merely tied**. Contrast `IGHV3-23`/`IGHV3-23D` and `IGHV3-30`/`IGHV3-30-3` at
   **1.0000**, where no read can ever separate them and a longer tie list is all that is on offer.
   ⚠ The *scaffold* set and the *call vocabulary* are separate decisions — scaffold ids are
   positional, so adding one renumbers the locus and invalidates every precompiled index. Ask
   before changing either.

7. ✅ **`arda.dpost` consumes what `arda scenarios` writes — `arda markup --d-prior PATH`
   (2026-09-25).** `load_d_prior(organism, path)` and `posterior_d(..., prior_path=)` take a
   table; the CLI flag implies `--d-posterior`, because a prior nothing reads is the failure
   mode this repo keeps hitting. ⚠ The shipped table is ALLOWED to be missing (that is what
   `None` means for 11 of 13 pairs); a path the caller typed is a request, so it **raises**
   rather than scoring on an empty prior. This separates *using* an estimate from *adopting*
   one, which is what item 8 below is the decision about — and item 8 is now approachable
   without installing anything. Original entry:

   **`arda.dpost` cannot consume what `arda scenarios` writes.** `docs/scenarios.rst` calls the
   output a drop-in for `d_prior.tsv`, and it is — by *format*. But `dpost.load_d_prior` is
   `@lru_cache`d on the organism and reads one fixed path (`dpost.py:108-110`), so the only way to
   use a fitted table is to overwrite a file inside the installed database. `arda.hmm.model_for`
   already takes `prior=`; `dpost` does not. Thread a path through `load_d_prior` /
   `posterior_d` and expose it as `arda markup --d-prior PATH`. Small, and it separates *using* an
   estimate from *adopting* one — which is the decision the entry below is about.

8. **11 of the 13 shipped (organism, D-locus) pairs have no `d_prior.tsv` at all.** Not a
   regression: OLGA has no model for them, which is the whole reason the table is derived rather
   than measured. Verified coverage —

   | organism | D loci with germlines | loci with a prior |
   |---|---|---|
   | human | IGH, TRB, TRD | IGH, TRB, TRD |
   | mouse | IGH, TRB, TRD | TRB |
   | rabbit | IGH, TRB, TRD | — |
   | rat | IGH | — |
   | rhesus_monkey | IGH, TRB, TRD | — |

   `load_d_prior` returns `{}` for a missing organism and `posterior_d` then returns `None`, so
   `arda markup` silently has no D posterior for three of five organisms. `arda scenarios` fits
   exactly this table from real junctions, so the blocker is **a cohort per (organism, locus)**,
   not code.
   ✅ **The A/B this entry asks for is DONE on human TRB, the only locus with both tables
   (2026-09-25).** 45,536 records fitted in 5 EM iterations (log-likelihood −1,057,744.6 →
   −976,796.3), judged against arda's own **nucleotide** D call at E ≤ 0.05 — independent of both
   priors — on 5,570 distinct clonotypes:

   | prior | agreement | TRBJ2 only | confident | confident agreement |
   |---|---:|---:|---:|---:|
   | shipped (OLGA) | .9339 | .9043 | .4736 | .9996 |
   | fitted (`arda scenarios`) | **.9363** | **.9076** | **.4876** | .9996 |

   Read the **TRBJ2** column: on TRBJ1 both the posterior and the nucleotide caller enforce the
   same TRBD2 × TRBJ1 prohibition, so agreement there is guaranteed rather than earned.
   ⚠ **The fit is in-sample** — same library — so +0.33 points is an upper bound. What it
   establishes is that a fitted table is **usable and not worse**, and that the path now works end
   to end. ⛔ **It also found a real defect**: `arda scenarios` writes a `#` provenance line above
   its header and `load_d_prior` skipped line 1 by POSITION, so the one file `docs/scenarios.rst`
   calls a drop-in raised on read. Fixed.
   **Still open, and it is DATA**: a cohort for the 11 pairs with no table at all. `aldan3`'s ngsik
   registry carries **3,616 M. musculus library rows**, which is where mouse IGH and mouse TRD
   would come from.

9. **`--error-rate`'s single default is wrong for variant preservation.** At the default `1e-3`,
   `rnaseq correct` erases both published MIGEC spike-in variants; `1e-5` recovers both exactly,
   and `1e-4` kept both while removing 72 % of real PCR errors on an independent cloud. Not a
   defect — no abundance method separates signal-to-noise ~1, which is why UMI consensus exists —
   but one constant cannot serve both regimes. Wanted: a per-library calibration *rule*, or an
   estimate off the data, not a re-tuned constant. ⚠ Whatever it becomes, it is not a QC threshold
   and must not turn into one.

10. ✅ **The SHM model is fitted and shipped as `arda shm-model` — on CONTEXT, not on
   per-allele position (2026-09-25).** `src/arda/shmmodel.py`, `project/design-shm.md`,
   `docs/shm.rst`; evidence in benchmark `results/round33`.

   ⛔ **The parameterisation this entry asked for was measured and REFUSED — do not re-propose a
   per-allele-per-position table.** Four bulk IGH libraries, two donors, on arda's own round-30
   annotation with no new alignment work:

   | object | within one donor | between donors |
   |---|---:|---:|
   | per-allele per-position profile | r .8578 – .9897 | **r .2557 – .5601** |
   | pooled positional profile | r .9708 / .9599 | r .5994 – .6224 |
   | **5-mer context** | r .9718 / .9785 | **r .7424 – .7854** |

   Inside one donor a positional table repeats at r ≈ .99 because the same expanded clones carry
   the same mutations at the same positions — `IGHV3-23*05` transfers at .9301 within SPX6730 and
   at .2557–.5601 between donors, with no change of method. ⚠ **And context is the only one of the
   two that reaches the positions the consumer needs**: the junction model wants
   `P(observed nt | V germline)` for the V tail *inside* the junction, which is exactly what
   `arda.shm` scopes out as unidentifiable. A position there has no clean data ever; a 5-mer takes
   its rate from every allele that carries it in framework sequence.

   ✅ **Two parameter families, both earned.** AID's WRCY/RGYW motifs carry **4.66× / 5.00× /
   4.77× / 4.69×** the rate of every other covered position across the four libraries, on an
   overall rate that itself moves 1.6× between the donors. And context does not exhaust the
   CDR:FWR contrast — fitting contexts on one donor and predicting the other's per-region counts
   leaves FWR1 **0.630 / 0.626** and CDR2 **1.389 / 1.291**, i.e. 4.65× raw becomes a 2.2×
   residual that reproduces between people who share no clones. The shipped estimator, run on the
   two donors independently, gives FWR1 **×0.683 against ×0.675**.
   ⚠ **The scale is the SAMPLE's, not the model's** (.031 against .050, same shape), so it is
   written as provenance and never applied.

   ✅ **S2 shipped: the junction model reads it (`arda scenarios --shm-model`).**
   `lattice(..., shm=)`, `accumulate`, `estimate` and the CLI flag. A templated stretch is scored
   by `Π(1 − μ)` over matches and `μ/3` over mismatches. ✅ **The exact-match bound turns out to
   be the `μ = 0` case of that emission**, not a separate rule — with every rate 0 no templated
   length past the common prefix survives, which is what `_common_prefix` computes, and a test
   compares term weights between the two paths. `shm=None` is the default and byte-identical.

   Measured on 26,619 real IGH junctions with the model fitted on a **different donor's** library
   (round 34):

   | parameter | exact | SHM | change |
   |---|---:|---:|---:|
   | `delV` mean | 4.148 | **2.708** | **−1.440** |
   | `delV` P(0) | .1426 | **.2388** | +.0962 |
   | `insVD` mean | 12.065 | **10.972** | **−1.093** |
   | `insDJ` mean | 12.342 | 12.144 | −0.198 |
   | `delJ` mean | 12.247 | 12.253 | **+0.006** |

   **The exact bound was charging 1.44 nt of V germline per IGH rearrangement to deletion** and
   1.09 nt to V-side insertion; `delV`'s mode moves from 1 to 0. ✅ The specificity is the
   result's own control — the model is applied to the V side only, so `delJ` at +0.006 against
   `delV` at −1.440 says the change reached what it should and nothing else. D posterior against
   arda's independent nucleotide caller: .7885 → .7898 on 10,849 junctions, **+15 genes net**,
   concentrated in the most mutated bin (+0.0040) and −0.0019 on unmutated. Cost +4.8 % wall.
   ⛔ **Never compare the two arms' log-likelihoods** — different models, the SHM arm carries an
   emission term per templated base.

   **Still open, and unrelated to SHM**: `_map_d`'s amino-acid path searches the three translated
   D frames as independent database entries, tripling `n`, when the prior over `insVD` already
   induces a prior over frame.

11. **Personalized germline — the consumer side shipped, the inference needs a confidence model.**
   `arda resolve-ties --genotype` applies an allele set (`v_call_genotyped`, `v_call` untouched,
   no re-alignment and no reference rebuild) and `arda genotype` infers one. `stats.py`'s
   `allele_candidate` is untouched and stays *"a shortlist to look at, never a call"* — the
   genotype is a separate, narrower object about the **reference**, never the repertoire.

   The call is a **likelihood ratio between diploid genotypes**. The junction is de facto a UMI,
   so a clonotype is one independent draw from the donor's two chromosomes and disagreement WITHIN
   a junction is error rather than allele -- which is where the miscall rate is measured from.
   Read and clonotype coverage are both reported. ⚠ The first version was TIgGER's frequency rule
   and had to be replaced: with no error model it called `TRBV11-2` off **754 of 757** clonotypes
   and `TRBV20-1` off **1 of 2,544**, reported both as `explained = 1.0000, ok`, and returned 43 of
   53 genes with none heterozygous.

   ⚠ **Open: read length is the binding constraint, not the rule.** Separating a TRBV gene's
   alleles needs a median of 150 nt from the 3' end; a 151 nt TRB amplicon covers a median of
   **72 nt** of V (56 after anchor-clipping, none reaching 150), so 33.1 % of clonotypes can be
   assigned and **17 of 53** genes are called -- 6 of them genuinely inferred. On a TRA amplicon at
   the same read length (`SRR5233635`, 21,710 clonotypes, error rate 5.53e-04) it is **20 of 44**,
   and one of them is the first heterozygous call: `TRAV36/DV7` = `*01`/`*04`, 78 clonotypes
   against 143, log10 BF **211**. Applying that genotype narrows **281 of 47,743** rows and
   contradicts 75 -- you cannot restrict what you could not genotype, and 14 of the 20 called genes
   have one catalogued allele. The whole inference is 1.53 s / 442 MB on 49,748 reads.
   The honest next step is to run this on a library that HAS the resolution (full-length, 5'RACE,
   or `arda cells` contigs at N50 536 nt) and measure against a known genotype; nothing about the
   rule can be judged on a library that cannot separate the alleles in the first place.

   Also open, and deliberately not in this cut: **novel-allele discovery** (needs the per-position
   SHM model of item 10 for IGH; TIgGER's y-intercept regression needs mutated reads, which TCR
   does not supply), **J-gene genotyping** (J targets are 38–69 nt and the framework-scoped span
   sits right on `TieResolver.MIN_SPAN`; measure before adding), and any **per-donor reference
   rebuild** — scaffold ids are positional, `build-db` needs IgBLAST, and the mmseqs freshness
   contract records no allele-set identity.


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
