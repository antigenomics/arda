# Roadmap

What arda does today, and what is planned. The maintainer's working backlog — measured dead ends,
per-round evidence, open questions — is kept out of the published tree.

## What works

**Reference build.** An offline IMGT + IgBLAST build of every in-frame V·J germline scaffold with
FR1–4 / CDR1–3 markup, for five organisms (human and mouse full IG+TR; rat, rabbit and rhesus IG),
with a per-organism locus-coverage manifest.

**Annotation.** MMseqs2 (or arda's own `segmap`) maps a query to that reference and a C++ hot path
transfers the markup through the alignment, emitting spec-valid AIRR Rearrangement rows: per-segment
V/D/J/C CIGARs, sequence and germline alignments, reverse-complement handling, out-of-frame junction
translation, and D-segment mapping including D-D fusions on IGH and TRD.

**Junction repair** (`arda markup`, `arda.cdr3fix`). For a bare `(junction_aa, V, J, species)`
record: confirm or re-call the V and J, repair the junction against the germline it names, and place
the V/J boundary in residues and in nucleotides. A side the record leaves blank is proposed from the
junction rather than refused — and where it names neither side, the locus is proposed too. The
repair policy is deliberately conservative — see the 2.33.0 entry in `CHANGELOG.md`.

**Scale.** Streaming, bounded-memory FASTQ I/O; multi-file samples (read groups); sharding across a
cluster via SLURM, Nextflow and Snakemake.

**Single cell** (`arda cells`): reference-free per-cell contig assembly, chain pairing, doublet
flagging, UMI counts.

**QC.** `arda stats` per run and `arda qc` across a cohort, with a self-contained HTML dashboard.

**A generative model of the rearrangement** (`arda scenarios`): EM over the recombination scenario
set, which fits the D prior arda ships. `arda.hmm` is the same recursion read as inference — score a
nucleotide junction, or ask which D made it. Pass `shm=` for **B cells**: without a model the
templated V length stops at the first mismatch, so a single hypermutation in the V tail makes the
whole tail read as N region.

## Where the boundaries are

arda is **the germline reference and the alignment against it**. Recombination probabilities — Pgen,
the most likely nucleotide reading of an amino-acid junction, the D posterior — belong to
[vdjtools](https://github.com/antigenomics/vdjtools), which reads arda for germline and markup. If
you want "which D, and where" from an amino-acid junction, call
`vdjtools.model.annotate_junctions`; it runs arda for the markup and the D alignment.

## Planned

- **Promotion of the broader TRDV × TRAJ reference.** The shipped default remains unchanged:
  the recorded real-amplicon comparison contains 680 truth reads of which it recovers 10 %.
  Version 2.37.0 repairs the opt-in construction and validates junction preservation;
  promoting a different default vocabulary is a separate reference-release decision.
- **A live multi-node SLURM run**; the sharding is unit-tested without a cluster.

## Release follow-up — 17 October 2026

[RELEASE_AUDIT.md](RELEASE_AUDIT.md) records the implemented structural-annotation
requirements and 60 focused checks. [Issue 149](https://github.com/antigenomics/arda/issues/149)
records the repaired opt-in TRDV × TRAJ construction and matched real-read validation
for 2.37.0: all 483 additional scaffolds pass coordinate checks, while all junctions
and J calls remain unchanged in the fresh 100,000-read comparison.
The productivity contract is implemented and documented in `docs/productivity.rst`: FR4 frame,
stops across the annotated span, and empty calls when the junction cannot be evaluated.
A start-codon requirement is not appropriate for partial reads; matching the junction tail
to germline is not a frame test.
They are not prerequisites missing from the repaired reference-index API.
