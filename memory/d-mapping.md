# D-segment mapping: interior, gate, genomic order, alphabets

Everything here was measured. Don't re-derive it; the ROADMAP records the dead ends too.

## 1. The interior must come from the anchors, not the projection

A scaffold is `V + 9 nt N-pad + J`. A real read has a 20–40 nt N-D-N region there, and
mmseqs — unable to align anything to a run of N — parks those bases against the flanking
V and J. Both projected boundaries then march inward. On the human fixtures the projected
`v_sequence_end` sits +6 nt too far and `j_sequence_start` −13 nt too early (IGH, medians),
leaving an 11 nt "interior" where the truth is 37 nt: **a window too small to hold an IGH D
at all.**

`transfer._anchored_vj_bounds` fixes this from the per-allele germlines in
`cdr3_anchors.tsv` — the V templates the junction's 5′ end, the J its 3′ end, so the bounds
are the longest common prefix / suffix. vs IgBLAST: `v_sequence_end` within 2 nt for 85 % of
IGH (projection: 43 %), `j_sequence_start` for 99 % of TRB (projection: 49 %).

**This — not "paralogous germlines + SHM" — was the cause of the old 46–69 % IGH D
concordance.** ROADMAP said the latter for months. It was wrong.

SHM truncates the exact match early, which *widens* the interior. That is the safe
direction: a D is never clipped, only surrounded by a little more sequence.

## 2. One E-value, not four score floors

Four hand-tuned per-locus floors (9/7/7/6 for huIGH/huTRB/huTRD/moTRB) became one knob:

    E = K·m·n·e^(−λS)   →   S_min = ceil( ln(K·m·n/E) / λ )

with `λ = ln((1−p)/p)` for `p` = the chance two residues match. **p = 1/4 recovers λ = ln 3**
for nucleotides. For amino acids p is *not* 1/20: N-region inserts and D germlines are both
G/S/Y-rich, and the measured p over real middles × real D frames is **0.0613 → λ = 2.7285**.

The absolute E is not calibrated (inserts are Markov, not uniform; K is borrowed). Treat it
as an m·n-corrected score and the threshold as an operating point. nt: `E ≤ 0.2`. aa: `E ≤
0.05` plus a floor of 4, because a 22–38 residue aa database leaves the E-value badly
under-calibrated at small n (at 0.2 the IGH false-call rate reaches 13 %).

`d_support` ships the E-value so a consumer can re-threshold.

## 3. Genomic order forbids TRBD2 × TRBJ1

The TRB locus runs `TRBD1 – TRBJ1 cluster – TRBC1 – TRBD2 – TRBJ2 cluster – TRBC2`, and
V(D)J joining *deletes* the intervening DNA. TRBD2 therefore cannot reach any TRBJ1. This is
genomic order, **not** a usage preference, and it holds in human, mouse, rat and rhesus.
IGH and TRD place every D 5′ of every J: nothing is forbidden there (verified — 0 structural
zeros, and IGH's 29 sub-1e-3 cells are low-usage genes, not a block).

Unenforced, TRBD2 (16 nt) simply outscores TRBD1 (12 nt) on noise: **17 % of real human TRB
J1-cluster D calls (18/104) and 12 % of mouse (7/59) were impossible**, at E-values (median
0.096) sitting in the chance band — against 0.014 for the producible TRBJ2 × TRBD2.

Two traps found while fixing it:

* **OLGA's human TRB model does not encode the constraint, and hides that behind an allele.**
  `TRBD2*01` correctly falls to ~1e-5 on J1 rows, but `TRBD2*02` absorbs 21–27 % of every J1
  row. At *gene* level — which is what `dpost` consumes — the unmasked model claims
  `P(TRBD2 | TRBJ1) ≈ 0.23`. The mouse model does learn it (~0.003). Mask, don't trust.
* **`dpost` backs off to the marginal when the J allele is outside the model, and the human
  model has no `TRBJ1-6*01` row at all** — so the mask must be re-applied after the backoff
  (`_mask_forbidden`), or TRBD2 walks back in through the door.

**Mask-only, never mask-and-rescore.** `_D_MAX_EVALUE` is an operating point calibrated
against the *full* locus D set, so shrinking `n` to the masked set silently loosens the gate:
measured, it admitted 58 new calls at median E = 0.098, squarely in the chance band. Holding
`n` fixed, a J1 record can lose an impossible call but never gain a weak one. Result: 18
human impossible calls gone, 3 recovering as TRBD1 at E ≤ 0.054, the 86 good D1 calls
untouched.

Where the accuracy actually landed: **all of it in the nt caller.** The aa posterior changed
*zero* calls on real data — the unmasked prior already favoured TRBD1 3.3:1 and the aa
evidence agreed. What masking buys `dpost` is honest confidence (posterior 1.0, entropy 0,
instead of 0.92/0.40 for a certainty) and an unrepresentable rather than merely unlikely
wrong answer.

## 4. aa input: three frames, one allele

A trimmed D has no knowable reading frame, so the aa reference carries all three translations
of every allele (`reference._load_d_germlines_aa`). Two frames of one allele can tie on one
span — `_best_d` de-duplicates, because that is one allele, not an ambiguity.

At E ≤ 0.05, floor 4, against a composition-preserving shuffled null (the "D trimmed away
entirely" null has n = 0..42 per locus and is useless):

| locus | call rate | gene acc | false call | ambiguous |
|---|---|---|---|---|
| human IGH | 69 % | 99 % | 2.0 % | 5 % |
| human TRB | 8 % | 100 % | 0.3 % | 0 % |
| human TRD | 12 % | 100 % | 0.0 % | 0 % |
| mouse TRB | 11 % | 100 % | 0.7 % | 58 % |

Out of model, on real fixtures: human IGH calls a D on 36 % of records where nt manages 68 %,
agreeing with it on 98 %. SHM and codon degeneracy cost half the recall, almost no precision.

Mouse TRBD1/TRBD2 translate to near-identical poly-glycine — hence 58 % ambiguity. That is
the honest answer, and it is why `_best_d` reports lists rather than picking a winner.

`d_germline_*` and `d_cigar` are **withheld** on aa: those offsets index a reading frame, not
the D germline.

## 5. D-D is limited by surviving length, not orientation

Order-respecting injection (TRBD1 5′ of TRBD2, the only producible orientation, and it can
only join a J2-cluster J): sensitivity 14.7 % human / 12.9 % mouse TRB, false `d2_call` on a
true single-D junction 0.1 % / 0.0 %. Stratified: **~0 % when either D survives under 6 nt,
57–72 % once both survive 7+.** Trimming usually leaves less. IGH reaches 61 %.

Orientation was never the limit: only 3 of 223 human calls came back as the wrong pair/order.

## 6. Dead ends (measured; see ROADMAP)

* Re-ranking nt D candidates by the length prior (`λS + log P(insVD) + log P(dlen) +
  log P(insDJ) + log P(D|J)`) changes **nothing**: IGH gene accuracy 98.9 → 97.8 %, huTRB
  94.2 → 94.5 %, identical call rate and `d_start` error. With 10–18 matched nt, `λS` is
  11–20 nats and the prior moves ±3. In protein the ratio inverts — which is exactly why
  `dpost` needs the prior and `_map_d` does not.
* A present/absent Bayes-factor gate buys IGH ~+3 pp recall at matched FP (93.6 % vs 90.7 %
  @ ~2 % FP) and nothing for TRD — but needs a **per-locus** threshold (BF > 6 for IGH,
  BF > 10 for TRD at the same FP), reintroducing the knobs the E-value removed, and the
  priors exist for only 4 of the 15 (organism, D-locus) pairs arda ships.
