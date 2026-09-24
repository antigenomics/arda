# Speed flags and the machinery behind them

Loaded on demand from `SKILL.md`. Each flag's measured regime, why the two configurations do not
compose, and the internals a change to the hot path has to respect.

## The regime rule

| regime | configuration |
|---|---|
| amplicon / RepSeq | `--two-pass --fast-segments --v-only-on-segment` |
| bulk RNA-seq | `--prefilter` |

They do **not** compose, and each is a loss in the other's regime. `--two-pass` alone is a loss
(0.762x on bulk, 0.87x on an IGH amplicon). The mode NAME picks the preset; `--exact` turns every
lever off; `arda map` still exposes all five individually for A/B work. The predictor is
`fast_fraction` in the run report -- the fraction of reads hitting both a V and a J segment --
not the library's name.

### `--prefilter`: the bulk lever

`arda map --prefilter` drops reads that share no exact 16-mer with the reference **before**
they reach MMseqs2. On bulk RNA-seq the receptor fraction is 0.02–3 %, so `mmseqs search` spends
essentially all of its time proving reads are *not* receptor reads; the fitted cost model says so
directly — `wall ≈ reads/46,353 + hits/350`, dominated by the **read count**, not the answer.

| sample | receptor % | before | after | speedup | vs TRUST4 |
|---|---|---|---|---|---|
| SRR10611239 | 0.024 | 82.58 s | 7.79 s | **10.6×** | **0.332** |
| SRR6926533 | 0.123 | 39.11 s | 5.26 s | **7.4×** | 0.603 |
| SRR8363894 | 0.772 | 46.86 s | 13.19 s | **3.6×** | 1.149 |

**Off by default**, because it costs ~0.5 % of real reads. The shape of that loss is what makes it
usable: every lost read scored **75–79 bits** against the `--min-score 75` cutoff (nothing
confident is ever dropped), `junction_aa` moved on **zero** shared reads, and the loss is
**entirely IG with zero across all four TR loci** — one substitution destroys k consecutive exact
windows, and SHM supplies substitutions. Check `prefilter_stats` (`seen`/`passed`) in the run
report to see whether it earned its keep on a given library.

**Measured against an independent IgBLAST truth on 344,554 real IGH-amplicon mates**, the loss is
now on the axis that causes it rather than a locus label: recall **.99209** vs the default's .99582,
and **1,283 of the 1,285 lost reads are in the `<90 %` V-identity bin** — 0 above 95 %, 2 in
90–95 %. On the cell lines the two reads it drops from Raji sit at 89.66 % and 77.78 % identity
(mean identity lost 83.72 vs 93.01 kept), and the 77.78 % read is the most mutated in the sample.
So it is safe for germline-proximal work and progressively costly with SHM load.

⚠ **Bulk only, and that is now quantitative.** Amplicon runs 46–90 % receptor, so almost everything
passes and the scan is pure overhead: on a 90 %-receptor IGH amplicon it is worth **exactly 1.00×**
(319.15 s vs 319.74 s) while still losing those 1,285 reads. Its value is confined to cold bulk
(**3.17×** on a 0.15 %-receptor library), where IG content is negligible anyway. ⚠ **The index is built from `Reference.target_fasta`, the same FASTA MMseqs2 searches.**
Never build it from a hand-listed segment set: against a `V+pad+J`-only reference the loss is
16.29 %, **69.27 % of it J→C reads**. Deriving it from the search target makes that unreachable.

⚠ **Do not reach for an MMseqs2 flag instead.** A 15-setting sweep found **no lossless candidate
above 1.05×**. `--min-ungapped-score 30` is free and gives *zero* speedup — which is the proof
that the cost is the k-mer stage, not the ungapped extension. And MMseqs2 can only prefilter reads
already in a DB, so the FASTA write and `createdb` (19.6 % of a 4 M-read run) are paid regardless.

⚠ **`--adaptive` and `--two-pass` are both WORSE on top of `--prefilter`** (measured: 11.13/16.01/
8.19 s and 14.69/22.84/9.17 s against 10.20/14.66/6.42). Once the scan term is gone, `--adaptive`'s
re-search is pure overhead, and prefiltered survivors being 84 % hitting does *not* put them in
`--two-pass`'s amplicon regime.

### The segment reference and the rescue guarantee

`build-index` also writes `segments.fasta` + `segments.markup.tsv` (`arda.refbuild.segments`):
every V allele, every J allele and every **C** allele as its own target — **924 human targets
(775 V + 124 J + 25 C) against 15,414 V×J scaffolds**, 16.7× fewer. It is *derived* from
`alleles.fasta` + `markup.tsv` in under a second, so it is neither committed nor shipped in the
release tarball, exactly like the mmseqs indexes.

The two-pass search uses it: hit the segment reference, take each read's best V and best J, look
the pair up in `combinations.tsv` — that names exactly **one** V×J scaffold, so the second
alignment is one target per read instead of ~277. A read running J into C names its J+C scaffold
the same way, through `Reference.jc_combinations()` — `(j_call, c_call) → scaffold id`.

⚠ **The constant region is its own target, and that was not always so.** Through 2.7.2 the 345 J+C
scaffolds were copied through verbatim, and they are a J×C product (IGH 14 J × 11 C, IGL 9 × 7,
TRB 16 × 2) in which every scaffold of a locus ends in the *same* constant sequence. A read
reaching C was therefore aligned against all of them to learn one `c_call`: measured on a TRA
amplicon, **27.7 % of the targets drew 76.4 % of the alignments**, 4,977 per target against 603 for
a V. Collapsing it (345 → 25) is **1.89× on the segment search and 1.33× on the whole `map`**, and
it *improves* the calls — a J call decided by a whole-scaffold bit score whose constant half is
arbitrary is a worse J call (`v_call` disagreements vs the one-pass 147 → 85, `j_call` 401 → 296).
A C target is kept for **every** locus, not only the informative ones: only IGH's 11 alleles
separate anything reportable (7 classes = the isotype), but a C target is also the sole segment
target a read lying wholly inside the constant region can hit, and dropping the other 14 makes 14
of 453 fixture reads vanish.

⚠ **`--two-pass` pays off when reads SPAN V INTO J — not by library type.** It needs a read to hit
both a V and a J segment, and the predictor is **`fast_fraction` in the report**, nothing else:

| library | receptor | fast path | vs one-pass |
|---|---|---|---|
| TCR amplicon (TRA) | 48 % | 85 % | **3.51×** |
| 5′-RACE human TRB | 100 % | **95.6 %** | **2.96×** |
| 5′-RACE mouse TRA | 100 % | 89 % | **2.64×** |
| 5′-RACE human **IGH** | 100 % | **16.3 %** | **1.03× slower** |
| IGH RepSeq amplicon | 90 % | **5.2 %** | **0.89× slower** |
| bulk RNA-seq | 2.74 % | 5 % | **0.762×, 31 % slower** |

Never: This was documented as "an amplicon optimisation, do not reach for it on bulk" and that framing
is **wrong in both directions**. The two 100 %-receptor rows are the *same dataset, same tool*: IGH
reads there cover V and stop short of the short IGHJ target, so they hit a V and no J and the fast
path collapses, while TRB reads span the junction. Read as amplicon-only, the flag gets left off
exactly where it is worth 3×. Bulk is separately a *scan-term* problem, which is `--prefilter`'s
job; this lever only touches the align term.

**Never: Every `fast path` number above is MMseqs2-specific.** `fast_fraction` is a property of
**(reads × segment mapper)**, not of the reads. On the same 100,000 IGH RepSeq pairs it is **0.052**
with the mmseqs segment search and **0.5018** with `--fast-segments` — 9.7×, from changing nothing
but the mapper, with `v_only` rescues falling 169,004 → 85,933. MMseqs2 misses the short IGHJ on
95 % of these reads and `_segmap`'s ungapped extension finds it on half. Re-read this table with
`--fast-segments` on before concluding a library is a bad fit for the fast path.

Off by default because there is no library type that predicts it — run a sample and read
`fast_fraction`.

### `--fast-segments`: the segment pass without a homology search

`arda map --two-pass --fast-segments` replaces the segment `mmseqs search` with
`arda._segmap`, a C++ seed → vote-by-diagonal → ungapped-extension mapper over the same
`segments.fasta`. It only **nominates**: every candidate is still aligned against the full
`V+pad+J` scaffold and scored by MMseqs2, so this is not a second aligner.

Why it can be exact: the segment pass asks only for the best V allele, the best J allele and their
coordinates (`_SEGMENT_FORMAT`) over a fixed 236 kb germline reference — **no cigar, no backtrace,
no gaps**. Germline V/J carry no indels relative to a read except sequencing error and IG SHM, so
ungapped extension is sufficient, and matches are near-identical to germline rather than remote
homologs.

| | segment step, 100 k amplicon reads, 8 threads | agreement with `mmseqs search` |
|---|---|---|
| `mmseqs search` | 2,770 ms | — |
| `_segmap` | **74 ms** | V allele .9997, J allele .9998, C allele 1.0000 |

End to end it is **1.49×** on `--two-pass` (50 k pairs, 9.10 s → 6.12 s) with `locus` identical on
every read, and on the 13-dataset benchmark panel **1.22× the default at identical recall, with
FP 62 → 59 and zero control FP** — better on every accuracy axis measured there.

**Where it is worth the most: IGH amplicon, 1.87× at 41 % less memory** (100 k pairs, 90 % receptor,
319.74 s → 170.77 s, RSS 4,016 → 2,382 MB), because it raises `fast_fraction` 0.052 → 0.5018 there.
It also **reverses** `--two-pass`'s bulk penalty: 0.73× → **1.24×** on a 0.15 %-receptor library,
where `no_segment_hit` is 99.87 % and those reads skip the full search entirely. `--two-pass` has
always been a structure-preserving prefilter that keeps the V/J call; it just was not worth paying
an `mmseqs search` for on bulk.

⚠ **Cost, measured against an independent IgBLAST truth on 344,554 real mates** (two IGH RepSeq
amplicons): recall **.99525** vs the default's .99582 — **197 reads**, of which **200 are in the
`<90 %` V-identity bin and 0 are above 90 %**. Those losses are *not* a seeding failure: `_segmap`
seeds more reads than mmseqs, not fewer (`no_segment_hit` 395,601 vs 395,682). It promotes ~100 more
reads onto the fast path, so a hypermutated V-only read gets seated on an implied V×J scaffold via a
weak J hit and the alignment then falls under `--min-score 75`. Same class as the round-5 narrowing
refutation — **do not "fix" it by lowering k.** Still off by default for that reason.

Requires the extension — check `arda.segmap.available()`. Without it the flag is a silent no-op, so
assert it in any job that claims to measure it.

**Never: `_segmap` CANNOT rank `V×J` scaffolds — only segments.** Pointed at the 15,414-scaffold
reference it indexes and maps fine and calls garbage (`v_gene` agreement **.3430**): a
junction-spanning read sits on a scaffold at **two** diagonals, V at one offset and J shifted by the
N-pad plus the non-templated junction, so one ungapped extension scores `max(V, J)` and never their
sum. The true scaffold therefore does not outscore a wrong one sharing its better-covered half.
That is the same fact that makes the segment reference work, read backwards — and it is why the
rescue stays MMseqs2's job. Do **not** retry it with chaining or a banded aligner "to fix the
diagonal": chaining two anchors across a non-templated insert is precisely the general-aligner cost
the whole approach exists to avoid.

Two constants make it equivalent, and both were wrong first:

- **`K = 12`**, the `-k` arda already passes MMseqs2 — *not* 16, which is `prefilter`'s value and is
  calibrated for **rejection**. Seed length sets sensitivity to mismatches: at k=16 the mapper
  seeded 53,048 reads against mmseqs' 53,121, and a read with no segment hit is assumed
  non-receptor and **never rescued**, so those 73 were lost outright.
- **`MIN_SCORE = 40`.** mmseqs applies `-e 1e-3`; this scheme has no e-value, so with no floor
  43,010 reads pick up a constant-region hit against mmseqs' 473 — half of them scoring exactly 38,
  a bare seed plus a couple of flanking matches.

### `--indel-rescue`: what an ungapped extension structurally cannot score

`arda map --two-pass --fast-segments --indel-rescue`. An ungapped extension follows ONE
diagonal, so a read carrying an indel relative to germline scores only up to the indel — measured
in the unit test: a 120 nt read scores **240** clean and **120**, exactly half, with one base
deleted. Its scaffold is then chosen on truncated evidence.

Measured on **341,294 real IGH mates** (IgBLAST gapped-alignment strings; a `-` in either row *is*
an indel, so nothing is inferred):

| V identity | `>=98` | `95-98` | `90-95` | `<90` |
|---|---|---|---|---|
| reads carrying a V indel | 0.74 % | 1.63 % | 3.56 % | **8.00 %** |

3.18 % pooled, tracking SHM load because AID makes indels and not only substitutions.

The signature is in the seed votes before any extension runs — two well-supported diagonals on the
**same** target, offset by the indel length — and the votes are already sorted by
`(target, diagonal)`, so one pass reads it.

**Rerouted, never dropped.** A flagged read is demoted from `implied` to `rescue` and realigned
gapped, so a false positive costs a little speed and *cannot* cost a read. That asymmetry is why
`MIN_DIAG_SEEDS` is deliberately low. Counted as `indel_rescued` in the report.

⚠ **Recall cannot measure this flag.** Recall is identical with and without it (174,066 of 174,226
amplicon fragments either way) because these reads are found regardless. What moves is the CALL,
which is load-bearing: the clonotype key is `(locus, v_call, j_call, junction)` at allele level.
Adjudicated against IgBLAST at gene level on exactly the reads whose call moved:

| sample | v_call moved | correct without | correct with |
|---|---|---|---|
| IGH_repertoire (hypermutated) | 586 | 53.24 % | **84.13 %** |
| IGH_naive (low SHM) | 100 | **77.00 %** | 63.00 % |
| **pooled** | **686** | 56.71 % | **81.05 %** |

**+167 reads, +24.3 points** — and the sign flip is the mechanism confirming itself. Where real
indels are common the flag fixes truncated calls; where they are rare its false positives (repeats
read as two diagonals) dominate. So **turn it on for hypermutated IG and leave it off elsewhere**:
on 13 bulk RNA-seq datasets it demotes **zero** reads and the output is byte-identical, which is
correct — bulk TR carries ~0 indels.

**Reads are never dropped by the fast path.** `arda.annotate.shortlist.shortlist()` partitions
every read into `implied` (took the fast path) or `rescue` (goes back to the full reference), and
asserts the partition is total. Anything that does not resolve — V only, J only, a V×J pair the
reference does not contain, a failed second alignment — is realigned, not discarded.

Traps, each of which produced *correct output* while quietly breaking something:

- **`JC|` targets are named by scaffold id, not by allele.** Feed the raw target name to
  `shortlist()` and every J→C read returns `no_such_combination`. Resolve through
  `Reference.segment_j_call()`. Measured cost of getting this wrong: the fast path collapsed
  from 85.3 % to 0.1 %. (Only reachable on a pre-2.8.0 reference; `C|` targets are named by allele.)
- **Never `top_hit` the segment pass.** One best hit per read destroys the V+J pairing the whole
  scheme depends on — `implied` goes to 0.
- **`_SEGMENT_SIDE` is ONE mapping shared by `_segment_rows` and `_segment_best_hits`.** Spell the
  target-kind rule out twice — once in polars, once in Python — and they drift: when `C|` was added
  to the loop alone, the polars reduction discarded every C row, `best_c` stayed empty, no
  constant-only read was rescued, and 15 J→C reads vanished **without `no_segment_hit` moving**,
  because the rows were gone before any counter saw them.
- **Nominate the J+C contest from the J, never from a C hit.** A J→C read spans the J/C boundary,
  so its constant overlap is often below the search threshold on its own even though the
  concatenated scaffold cleared it. Gating on C evidence re-invented `TRBV12-3*02`, destroyed a
  `c_call`, and fabricated a `junction_aa` on a read the one-pass calls V-less.
- **`segments.fasta` is generated, and a deploy does not regenerate it.** A new mapper reads stale
  `JC|` targets through its back-compat path and *succeeds* while measuring the old reference.
  Assert the target composition, not the file's existence.

**TRD is TRAV/DV + TRDJ; the J (and C) decides the locus, not the V.** arda's reference already
encodes it — of 1,050 scaffolds built from a TRAV/DV segment, 1,005 are locus TRA (with a TRAJ)
and 45 are locus TRD (with a TRDJ). There is deliberately no "α/δ is ambiguous" rule:
`combinations.tsv` is the arbiter, a pair it contains is real biology and one it does not is a
genuine chimera.

### `annotate.project`: the junction from coordinates, no alignment

```python
from arda.annotate.project import project_junction, UNVALIDATED_LOCI

proj, refusal = project_junction(strand_seq, qlen, v_row=v_row, j_row=j_row,
                                 v_call=..., j_call=..., anchors=ref.anchors,
                                 split_checked=...)   # REQUIRED keyword, deliberately no default
if proj:
    proj.junction   # AIRR: Cys104 .. [FW]118 INCLUSIVE
    proj.cdr3       # IMGT: anchors excluded, two residues shorter
```

AIRR `junction` runs Cys104 -> [FW]118. The segment pass already returns `(target, tstart, qstart)`
per read per side and `cdr3_anchors.tsv` already records `anchor_nt`, so the position is arithmetic:

```
offset = (anchor_nt + 1) - tstart
pos    = qstart + offset                    # forward
pos    = (qlen - qstart + 1) + offset       # reverse complement
```

**Never: Three coordinate systems disagree about their origin here.** `tstart` is **1-based** on the
forward target (`segmap.cpp`), `anchor_nt` is **0-based** in the germline (`cdr3fix.Anchor`), and a
minus-strand `qstart` is in **forward** coordinates while the sequence it indexes is the reverse
complement. Each is an off-by-one that still yields a plausible-looking junction — right length,
starts with a codon, ends with a codon. `strand_seq` must be the strand the hits were MEASURED on;
pass `revcomp(read)` for a minus-strand hit rather than reflecting coordinates afterwards.

Measured vs IgBLAST at `v_score >= 70` on 254,867 reads: IGH .99977 (naive) / .99634 (91.77 % median
V identity), TRA .99947, TRB .99949, IGK and IGL exact. >= .993 in every V-identity stratum.

**It refuses rather than degrades** — `no_anchor`, `unvalidated_locus`, `indel_split`,
`strand_mismatch`, `order`, `off_read`, `bad_codon`. A well-formed junction that is wrong is worse
than no junction: the reference-geometry bug shipped junctions that started `C`, ended `[FW]`,
passed `--complete-only`, and were short by exactly the allele's truncation.

**Never: `UNVALIDATED_LOCI = {"TRD"}`.** Not because TRD is known bad — because it has **zero** coverage.
Across two TR amplicons the segment pass never handed a TRD read both anchors, so all 767 TRD
truth junctions went to the aligner and TRD never appeared. Absent reads like fine in every
aggregate. The locus is taken from the **J** anchor, never the V: TRAV/DV rearranges to either TRAJ
or TRDJ and the J decides.

**Never: A `[FW]GXG` motif check is NOT equivalent to reading `anchor_nt`.** `TRAJ35*01`'s anchor codon
decodes **Cys (TGC)** — it is `status = ok` and a functional IMGT `F` gene, with a real Cys six
codons past the FGXG. A motif check deletes the gene silently (33/33 amplicon reads lost).

⚠ Yield, not accuracy, sets the scope: 87.0 % / 77.3 % of hit amplicon reads carry both anchors
against **7.2 % of bulk** reads, which mostly do not span a junction at all.
