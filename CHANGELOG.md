# Changelog

What changed for you, per release. Anything not listed is internal.
Full release notes: <https://github.com/antigenomics/arda/releases>.

## 2.37.0

**One CPU allocation across mapping phases.** Automatic thread selection respects
process affinity, Linux cgroup CPU quotas and SLURM allocations. Streaming mapping reserves one CPU for
its reader and reads synchronously with a one-CPU budget. MMseqs preparation and
reduction no longer default to every host core; child OpenMP limits prevent nested
teams. The CLI initializes a one-thread Polars pool; embedded mapping rejects an
already-initialized table pool that exceeds its selected CPU budget. A bounded reader stops on
consumer failure as well as end of input. Python callers can pass `threads` to
`annotate_sequences` while retaining the existing batch interface.

**Opt-in TRDV × TRAJ reference construction.** Missing J markup can be restored
from exact scaffold geometry and existing germline anchors. The opt-in builder
uses anchor-consistent J frames and refuses to package incomplete TRDV × TRAJ
coverage. This changes neither the default pairing policy nor the shipped
reference; real-read recall validation is tracked in issue #149.

**Documentation and examples.** Navigation separates tutorial, how-to, reference
and explanation. The performance guide describes CPU ownership and resource
measurements; embedding examples correctly describe internal search processes
and temporary files. The cost-model benchmark uses the current `arda map`
command. Junction example regeneration no longer passes retired D-posterior flags;
markup TSV uses the same unquoted dialect as the other writers. Polars stays below its unvalidated major version 2.

## 2.36.2

Requires seqtree 1.1.1 or newer within major version 1. The productivity roadmap now
points to the implemented FR4/frame/stop contract instead of promising a separate rule set.
TRDV × TRAJ reference completeness remains an explicit follow-up; this release does not
change the existing opt-in pairing policy or claim a junction-recall improvement.

**Reference indexes are checked by content.** Source discovery recognizes ARDA's own
source module rather than another project's `database/` and `pyproject.toml`.
Explicit `ARDA_HOME` remains authoritative. Legacy timestamp-only MMseqs indexes
are preserved but rebuilt into immutable generations. Reuse validates source FASTA,
MMseqs version and every sequence/header/index file, including split payloads.
Atomic generation publication keeps concurrent readers on coherent files.

`arda.reference_index_identity(organism, seqtype="aa")` resolves and validates the
actual searched index and returns SHA256 provenance for its source and payloads.

## 2.36.1

**Shared TRA/TRD V alleles retain both locus memberships.** `load_anchors` previously
kept one `Anchor.locus` per allele, so the later TRD row hid the TRA membership of
15 human and 42 mouse alleles, including `TRAV38-2/DV8*01`. `Anchor.loci` now retains
all declared memberships; candidate selection uses them, while `Anchor.locus` remains
available for existing callers.

**A shared V uses its named J to determine the chain**, including accepted aliases
such as `TRAV38-2` and `DV8`. Junction markup, scenario construction and D calling
resolve aliases before choosing the locus, so a shared V with `TRDJ` remains TRD.
With no J named, the canonical shared V continues to default to TRA.
Comma- and `+`-joined call lists consistently use their leading candidate in all
three consumers; later shared candidates do not change the leading V's locus.

Regression checks cover membership and reversed reference-row order across all five
bundled organisms, plus human and mouse aliases through all three consumers.
No reference sequences or coordinates changed.

## 2.36.0

**A junction naming NEITHER V nor J now gets a locus proposed too.** 2.34.0 proposed the missing
side from the junction, but only within the locus the *other* side named — so a record naming
neither came back refused, with no locus, no boundary and nothing for a nucleotide stage to work
on. There is no new rule here: every locus the organism ships now competes, scored by the same
`anchor_depth` from each germline's own anchor, and a locus only wins by explaining residues at
**both** ends. Both sides must name something; half an explanation is not a locus.

522 of VDJdb's curated `chunks` records name neither side (461 distinct keys). All 461 now get a
locus and **459 come back `good`**, against none before. The proposed locus agrees with the
`cdr3.alpha` / `cdr3.beta` column the record was filed under on **457 of 461 (99.13 %)** — 318 TRB,
139 TRA — and all four disagreements are `CACD…DKLIF`, TRDV2's own anchor and TRDJ1's own ending,
in a schema with no delta column. `proposed` reads `V,J` on every one, so a curator can see exactly
which record had nothing to go on.

Each end must still explain its **own** anchor residue, or no locus is named: without that floor a
junction agreeing with nothing (`QQQQQQQQQQQQ`) was handed one. Over the 461 real keys the winning
locus clears it on every one.

**`scripts/audit_cdr3fix.py` now keeps blank-side keys**, which it had been dropping. It required
both a V and a J to be named, so the A/B instrument was blind to exactly the rows 2.34.0 and this
release change. The corpus goes 189,596 → **192,726** keys and the digest rebases to
`0df8541ed1060fe1`; on the 189,596 keys the old filter kept, the digest is **unchanged** at
`2b75491b380310aa`, so nothing that already worked moved.

## 2.35.0

**`arda.hmm` is no longer deprecated — it is the B-cell entry point.** 2.33.0 deprecated it on the
grounds that nothing consumes it and `vdjtools.model.infer_nt_batch` answers the same question
faster. That surveyed the T-cell path and missed the `shm=` parameter.

Without a somatic-hypermutation model the templated V length is bounded by an **exact common
prefix**, so a single substitution in the V tail forces the whole rest of it to be re-read as N
region — measured on `IGHV3-30*18` / `IGHJ4*02`, `del_v` goes from `0` with a model to
`>= len(v_nt) - 3` without one. For a hypermutated IGH junction that is the normal case, and the
replacement recommended in 2.33.0 has no SHM term at all, so it prices a mutated V tail as insertion
in exactly the same way. The deprecation warning is gone and the module stays.

Unchanged: `arda.hmm` is still not on the annotation path, and the two measured negatives behind that
still stand — re-ranking nucleotide D candidates by a scenario likelihood changes nothing, and
replacing the E-value gate with a Bayes factor would need a per-locus shipped threshold. Both are
statements about germline TCR junctions, where an alignment already settles it.

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

**Deprecated — `arda.hmm`**, with a warning. *(Reversed in 2.35.0 — see that entry.)*

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
