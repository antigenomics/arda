# Recombination scenarios from a nucleotide junction

2026-09-24. Implements `ROADMAP.md`'s *"Nucleotide junction re-mapping → recombination-scenario
counts"*. `project/design-singlecell.md` and `design-qc.md` are the house style this follows.

## What closes the loop

`arda.dpost` places and identifies a D from an **amino-acid** junction by marginalising a
generative model. That model ships as `database/vdj/<org>/d_prior.tsv`, and every number in it is
**borrowed from OLGA / vdjrearm**. The roadmap item is not "add a feature": it is *stop borrowing*.

The table is already the right shape — long, `locus / kind / key / value`, the same idiom
`stats.py` uses — and `load_d_prior` reads exactly five counted kinds:

| kind | key | what it is |
|---|---|---|
| `insVD` | `<n>` | P(number of non-templated nt between V and D) |
| `insDJ` | `<n>` | P(number between D and J) |
| `dlen` | `<d_allele>:<n>` | P(surviving D length \| allele) |
| `d_marginal` | `<d_allele>` | P(D) |
| `d_given_j` | `<d_allele>\|<j_allele>` | P(D \| J) |

(`beta` is a fitted temperature, not a count, and is left alone.)

So the deliverable is: read `(junction_nt, v_call, j_call)` records, count recombination
scenarios, and write a table of that exact shape. It is then a **drop-in** for the shipped file.

## Why the counts cannot be read off one reading of the junction

`annotate.dmap.map_d_junction` already gives *one* reading — `v_sequence_end` from the longest
common prefix against the V germline, `j_sequence_start` from the longest common suffix, and a
D called under an E-value gate. That is a MAP-ish point estimate and it is **not** what an E-step
wants.

Never: **a scenario is not identifiable from sequence.** Exonuclease chew-back and N/P addition
mean several `(delV, insVD, delDl, delDr, insDJ, delJ)` explain one junction exactly — a first
non-templated base that happens to match germline is indistinguishable from one less nucleotide of
trimming. Counting the MAP reading assigns weight 1 to one member of that set and 0 to the rest,
which biases every trimming and insertion distribution toward *less* trimming and *shorter*
inserts. The same ambiguity is why `CLAUDE.md` forbids "fixing" a V/J boundary inside a junction.
So the tallies are **expected counts under the current model**, summed over scenarios.

Having a weight requires a model; a model plus expected counts plus renormalisation *is* EM. That
is why this stage produces an estimator rather than a counter, and it is the E-step
`arda.hmm` (ROADMAP's next item) will reuse — the forward-backward pass over the same state space.

## S0 — the scenario set

`src/arda/scenarios.py`, pure: no I/O, no new dependency.

```python
@dataclass(frozen=True, slots=True)
class Scenario:
    del_v: int; ins_vd: int; d_call: str
    del_dl: int; del_dr: int; ins_dj: int; del_j: int
```

Junction space throughout, as in `cdr3fix` and `dmap`: Cys104 → [FW]118 inclusive. `Vg` and `Jg`
are the **untrimmed** junction-region germlines already shipped per allele in
`cdr3_anchors.tsv` (`germline_nt`), so nothing here reads a FASTA.

A scenario is admissible when the junction is reproduced exactly:

```
junction  ==  Vg[:len(Vg)-del_v]  +  N1  +  Dg[del_dl:len(Dg)-del_dr]  +  N2  +  Jg[del_j:]
                                    └ ins_vd ┘                            └ ins_dj ┘
```

Enumerated by **D placement first**, which is what prunes it. For each D allele and each start
`p` in the junction, extend the longest exact match at every germline offset; each surviving
`(p, m, del_dl, del_dr)` then admits `v_len` in `0..min(v_max, p)` and `j_len` in
`0..min(j_max, L-p-m)`, with `ins_vd = p - v_len` and `ins_dj = L - p - m - j_len`. `v_max` /
`j_max` are the longest common prefix / suffix — the same two quantities `map_d_junction` already
computes, reused rather than recomputed.

Never: **`m == 0` is admissible and common**, not a failure. The shipped table's own most massive
entry is `IGHD1-1*01:0 = 0.8759` — the D survived nothing. A VDJ locus always rearranged a D;
"no D found" means it was trimmed away, and dropping those records would truncate `dlen` at 1 and
inflate every insertion distribution by the length of the D that was actually there.

Never: **at a VJ locus there is no D and no `insDJ`.** The middle is one `insVD` run, and the
locus contributes to `delV` / `delJ` / `insVD` only. TRA/TRG/IGK/IGL are not skipped — they carry
most of the trimming evidence.

## S1 — expected counts, and the EM loop

One pass accumulates, per record, `P(scenario | junction)` under the current parameters and adds
it to every sufficient statistic the scenario touches. Renormalise; repeat.

Initialisation is the **shipped prior** where it exists, so iteration 1 is already a sensible
E-step rather than a uniform one. `delV` / `delJ` have no shipped distribution (nothing consumed
them), so they start uniform over their feasible range and EM moves them.

Never: **weight by the model, never by 1/n.** Uniform weighting over an ambiguity set is itself a
strong and wrong prior — it says a 12-nt insertion is as likely as a 2-nt one.

Never: **a record contributes 1.0 in total, not one count per scenario.** A junction with 400
admissible scenarios must not outvote one with 3.

Never: **an insertion must cost its own sequence, not just its length** — `P(seq | len) = 0.25^len`.
Without it the junction's 5' end is explained either as templated V or as an insertion that
happens to match V germline, and the second is *free*, so EM walks into the corner where
everything is insertion. Measured before the term existed: three iterations on 503 real human TRB
junctions moved `insVD` mass onto **10–11 nt** with the log-likelihood rising monotonically the
whole way — EM doing exactly what it was asked. With the term, the same data gives `insVD` peaking
at **4 nt** and `dlen` for `TRBD1*01` peaking at **4–5 surviving nt**, which independently
reproduces the *"median surviving D is 5 nt for human TRB"* figure in `dpost.py`'s own docstring.

## S2 — `arda scenarios`

```bash
arda scenarios -i sample.clones.tsv -o d_prior.tsv --organism human
```

Reads `junction` / `v_call` / `j_call` (the clonotype table's own names, so a run's output feeds
straight in) and writes the five shipped kinds plus `delV` / `delJ`, keyed `<allele>:<n>`.

Never: **abundance-weighted, not row-weighted, by default.** A clonotype table row is a clonotype,
not an observation of the recombination process. `duplicate_count` is the number of observations.
⚠ but a clonal expansion is one recombination event seen many times, so `--weight rows` is there
and the choice is stated in the output header rather than defaulted silently.

## Deliberately not here

| not this stage | why |
|---|---|
| `P(V)`, `P(D,J)` full IGoR parameterisation | nothing in arda consumes them; `d_prior.tsv`'s five kinds are what `dpost` reads |
| `delDl` / `delDr` as aggregate kinds | unidentifiable when `m == 0`, which is the modal case; kept per-record where a placement is concrete |
| SHM-aware scenarios | a mutated IGH junction breaks the exact-match premise; this stage is for TR and unmutated IG |
| replacing the shipped `d_prior.tsv` | generating one is this stage; **adopting** it is a measurement and a release decision |
| a per-record MAP scenario artifact | `enumerate_scenarios` is public and returns the whole set for one junction, which is what inspection wants; a MAP column would be the one reading this module exists to argue against |
| a fitted insertion composition model | uniform `0.25/base` is what breaks the degeneracy (see below); a first-order Markov is a refinement and one more table initialised from nothing |
| `arda.hmm` | ROADMAP's next item, and the E-step here is its forward-backward pass |
