# Changelog

What changed for you, per release. Anything not listed is internal.
Full release notes: <https://github.com/antigenomics/arda/releases>.

## 2.34.1

**Faster, with identical answers.** `_extend` — the function that walks a germline run, and the
hottest one in `markup_batch` at 159,240 calls per 20,000 junctions — used to walk the same span
three times. It walks it once now: **16,287 → 18,338 junction keys/s**, so VDJdb's whole 189,596-key
corpus marks up in 10.3 s in one process. The audit digest over every one of those keys is unchanged.

**Two invariants are now tests rather than accidents.** Every residue `v_end` / `j_start` credits is
germline of the allele it names — measured at 185,636 of 185,636 V keys and 187,272 of 187,272 J keys,
zero disagreements — and `cdr3_repaired` is a fixed point of itself. Separately, a functional V opens
with Cys104 and a functional J closes with Phe/Trp118 on 3,129 of the 3,130 functional entries across
all five shipped organisms; the one exception, `TRAJ35*01`, is real IMGT germline and is named in the
test rather than repaired away.

## 2.34.0

**New — a blank V or J call is proposed from the junction instead of refused.** A submission is
allowed to leave one side out, and a blank is not a reason to refuse the record: the locus comes
from the side that *is* named, and the junction is evidence about the missing one. The new
`Cdr3Markup.proposed` (and the `proposed` column) says which side was never curated — a different
fact from `allele`, which means the submission named a different allele of the same gene.

Over VDJdb's 192,726 distinct curation keys, 3,130 leave a side blank (644 no V, 2,947 no J); 2,532
now get a segment, 2,504 of them `good` with both boundaries placed. The 598 that stay refused name
neither side, so there is no locus to propose within.

An **unresolvable** call is still refused: `TRBVnope*01` keeps its `FailedBadSegment`, because a
submission that names something wrong has a defect a curator must see, where one that names nothing
has a gap the junction can fill.

## 2.33.0

**Breaking — `arda.dpost` moved to vdjtools.** A posterior over the D gene is a
recombination-model question, so it now lives with the model: `vdjtools.model.posterior_d`, plus a
new batch form `posterior_d_batch` (one row out per row in, in input order). The module was ported,
not rewritten, so its answers are unchanged. `arda markup --d-posterior` and `--d-prior` are gone
with it. arda still ships the prior table and still fits one (`arda scenarios`).

**Breaking — `cdr3fix` repairs less and reports more.** The rule is now *when it is not certain,
flag it, do not change it*. Four edits are admissible: re-call the V/J when another allele explains
two more residues contiguously from its own anchor; trim framework past an anchor; restore germline
the submission was cut inside of; and substitute the **first or last residue only**, to the anchor.
A residue disagreeing with germline *inside* the templated run is reported with the new `mismatch`
flag and left exactly as submitted — a curation error and an allele IMGT does not record cannot be
told apart from a single junction.

What you may notice: `good` no longer trips on an internal disagreement (it means "neither side has
an anchor the rules cannot restore"); `max_replace` bounds the anchor only; and a junction that lost
just its Cys104 or Phe/Trp118 now gets it back.

Against the VDJdb 2026-06-03 release over all 187,488 curation keys, the repaired junction agrees on
**99.84 %** (was 98.36 %), and arda ships **677** non-canonical junctions where that release ships
715.

**New — `v_alts` / `j_alts`.** The alleles a junction cannot separate travel with the answer, chosen
one first, instead of being resolved by name order. Hand them to `vdjtools.model.infer_nt_batch`,
which scores a list of alleles per row.

**New — `map_d_junction(v_end=, j_start=)`** takes the interval to search instead of re-deriving it.

**Deprecated — `arda.hmm`**, with a warning. For the most likely nucleotide reading of a junction,
use `vdjtools.model.infer_nt_batch`.

## 2.31.0

The germline V/J boundary in **nucleotides**: `Cdr3Markup.v_end_nt` / `j_start_nt`, and
`cdr3fix.boundary_nt`. An exonuclease does not stop on a codon boundary, so a residue-level boundary
rounds; against external nucleotide truth `v.end` goes 71.8 % → **92.9 %** exact.

## 2.30.0 – 2.30.1

`v_end` / `j_start` stop where germline evidence stops rather than where the aligned extent ends:
`j.start` exact on 8,164 of 8,334 records against 8,041 before. Of the two, 2.30.1 is the one that
reached PyPI.

## Earlier

2.5.0 – 2.29.0 are described in their git tags (`git tag --sort=-v:refname`) and in the GitHub
release notes.
