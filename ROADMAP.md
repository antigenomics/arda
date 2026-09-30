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
the V/J boundary in residues and in nucleotides. The repair policy is deliberately conservative —
see the 2.33.0 entry in `CHANGELOG.md`.

**Scale.** Streaming, bounded-memory FASTQ I/O; multi-file samples (read groups); sharding across a
cluster via SLURM, Nextflow and Snakemake.

**Single cell** (`arda cells`): reference-free per-cell contig assembly, chain pairing, doublet
flagging, UMI counts.

**QC.** `arda stats` per run and `arda qc` across a cohort, with a self-contained HTML dashboard.

**A generative model of the rearrangement** (`arda scenarios`): EM over the recombination scenario
set, which fits the D prior arda ships.

## Where the boundaries are

arda is **the germline reference and the alignment against it**. Recombination probabilities — Pgen,
the most likely nucleotide reading of an amino-acid junction, the D posterior — belong to
[vdjtools](https://github.com/antigenomics/vdjtools), which reads arda for germline and markup. If
you want "which D, and where" from an amino-acid junction, call
`vdjtools.model.annotate_junctions`; it runs arda for the markup and the D alignment.

## Planned

- **Retire `arda.hmm`** (deprecated in 2.33.0). Nothing consumes it and vdjtools answers the same
  question faster.
- **Junction recall on TRDV × TRAJ rearrangements**, the largest remaining named gap against MiXCR
  on a real amplicon: 680 truth reads of which arda currently recovers 10 %.
- **Full AIRR productivity rules**, beyond the stop-codon and frame checks shipped today.
- **A live multi-node SLURM run**; the sharding is unit-tested without a cluster.
