# A model of somatic hypermutation — what it is parameterised on, and why

`ROADMAP.md` item 10. `arda.shm` reports **where** a read differs from its germline; this is the
model of **where a difference is expected**. Two entries wait on it: the V → N1 → D → N2 → J
model (`arda.scenarios` / `arda.hmm`) cannot be pointed at IGH without one, and novel-allele
discovery needs one because TIgGER's regression runs on mutated reads.

Status: **S1 shipped** (`src/arda/shmmodel.py`, `arda shm-model`). S2 and S3 below are not
started and each needs a measurement before it is written.

## The decision this document exists for

The roadmap asked for **per-allele-per-position**. Benchmark `results/round33` measured whether
the data can carry that, on four bulk IGH libraries from two donors, using arda's own round-30
annotation and no new alignment work.

| object | within one donor | between donors |
|---|---:|---:|
| per-allele per-position profile | r .8578 – .9897 | **r .2557 – .5601** |
| pooled positional profile | r .9708 / .9599 | r .5994 – .6224 |
| **5-mer context** | r .9718 / .9785 | **r .7424 – .7854** |

⛔ **A per-allele-per-position table is a portrait of one donor's clonal history.** Inside a donor
it repeats at r ≈ .99 because the same expanded clones carry the same mutations at the same
positions. `IGHV3-23*05` transfers at .9301 within SPX6730 and at .2557–.5601 between donors, with
no change of method. Do not re-propose one as a shipped artifact.

⚠ **And the consumer needs rates where a positional model has none by construction.** The junction
model wants `P(observed nt | V germline)` for the V 3' tail **inside** the junction — exactly the
positions `arda.shm` scopes out, because a difference there is hypermutation, chew-back or
N-addition and sequence cannot say which. A germline position past Cys104 has no clean data ever.
A 5-mer does: its rate is estimated from every allele that carries it in framework sequence. The
model is estimated where the evidence is and applied where the consumer needs it, and the transfer
number above is what licenses that.

## The two parameter families, and why each is earned

**Context**, because it transfers (above) and because AID's targeting is a context property:
positions matching WRCY on the sense strand or RGYW on the other carry **4.66× / 5.00× / 4.77× /
4.69×** the rate of every other covered position across the four libraries — on an overall rate
that itself moves 1.6× between the donors. Three of the four hottest 5-mers in SPX6730a carry
`AGCT`, the canonical WRCY.

**A region multiplier**, because context alone does not exhaust the CDR:FWR contrast. Fitting
contexts on one donor and predicting the other's per-region counts leaves:

| region | fit SPX6730 → test SPX8151 | fit SPX8151 → test SPX6730 |
|---|---:|---:|
| FWR1 | 0.630 | 0.626 |
| CDR2 | 1.389 | 1.291 |

Raw CDR2-over-FWR1 is 4.65×; the residual after context is 2.2×, and it reproduces between people
who share no clones. ✅ The shipped estimator, run independently on the two donors, gives FWR1
**0.683 against 0.675** — agreement to 0.008 on a parameter fitted from different individuals.

**No third family.** Nothing else is measured to be necessary.

## What is deliberately not a parameter

| not this | why |
|---|---|
| the overall rate | it is the **sample's**, not the model's: .031 against .050 between two donors with the same shape. `scale` is written as provenance and a consumer rescales to the sample in front of it |
| a per-allele-per-position table | measured above; it does not transfer |
| the substitution destination (`P(base | mutated)`) | the entries carry it (`G279C`) and collecting it is nearly free, but nothing reads it yet. A consumer that needs `P(x | g)` rather than `P(mutated)` adds a fourth column then — uniform 1/3 until one exists |
| a shipped table in `database/` | same stance as `d_prior.tsv`: fitting a model and adopting one are separate decisions, and the fit is per sample |
| a mutability model for J | J targets are 38–69 nt and the framework-scoped span is short; measure before adding |
| TR loci | SHM is an IG property. A TR fit measures allele mismatch in the templated framework, which is real but is not this |

## S1 — `arda shm-model` (shipped)

```bash
arda shm-model -i mapped.airr.tsv -o shm.tsv --locus IGH
```

Reads what `map` / `amplicon` / `rnaseq` already write under the default `--shm framework`.
Output is one `#` provenance line then `kind` / `key` / `value` rows: one `scale`, one `context`
per 5-mer clearing `--min-context`, one `region` per estimated region. 0.76 s on 4,661
observations.

Never: **the denominator is capped at `v_anchor_nt`.** `v_mutations` is framework-scoped, so
counting every covered position would divide a scoped numerator by an unscoped denominator and
report a rate too low by however far the read runs past Cys104. This is the single easiest thing
to get wrong here.

Never: **a tie list is dropped, not resolved.** `IGHV3-23*01,IGHV3-23D*01` does not say which
germline the read came from, so its differences cannot be attributed to one of them.

Never: **`--weight unique` is the default and the choice is written into the header.** A bulk read
carries no junction, so there is no clonotype key to collapse on; a distinct
`(call, span, mutation set)` tuple is the closest available proxy for one molecule, and counting
reads lets one expanded clone vote thousands of times.

Never: **a region with no observed substitution is OMITTED, never given a multiplier of 0.** Zero
does not read as "thin evidence here", it reads as "a substitution in FWR1 is impossible", and
`rate()` would then return 0 for every position in that region — correct-looking output, silently
wrong. An omitted region takes the default multiplier of 1.

Never: **the loader skips comments, blanks and the header by WHAT THEY ARE, not by position.**
`load_d_prior` dropped line 1 by position and broke the day a generator wrote a provenance line
above the header. `tests/unit/test_shmmodel.py` pins it with a file that has a comment above the
header, a second comment lower down, a repeated header mid-file and a malformed row.

## S2 — the junction model reads it (not started)

`scenarios.lattice` bounds the templated V length with `_common_prefix(junction, v_nt)`, an
**exact** match. Under SHM a substitution in the templated tail truncates that bound and the
remainder is explained as insertion. The roadmap's own note says the failure is *safe* — it
widens the interior and never clips the D — which is why this is an upgrade, not a bug fix.

The change is to score a templated stretch by `Π (1 − μ)` over matches and `μ · 1/3` over
mismatches, with `μ` from `ShmModel.rate(germline, pos)` and no region term (the tail belongs to
no V region). ⚠ **Done means measured, not written**: the existing D-call accuracy and junction
recall on IGH must not move against round 28's numbers, and the two negatives the roadmap already
records must not be re-run.

## S3 — novel-allele discovery (not started)

TIgGER's y-intercept regression, which needs S1's rates as the expected-mutation baseline: a
position whose mutation frequency stays high in reads with *few* other mutations is a
polymorphism, not a hotspot. ⚠ This is also the reason a per-position table is not just
untransferable but actively misleading — the donor's own alleles show as 100 %-mutated positions,
which is item 11's payload and this model's contaminant.
