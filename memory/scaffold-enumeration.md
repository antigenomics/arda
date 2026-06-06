# Scaffold enumeration: V×J only, dedup, N-spacer

## Decision
The reference DB enumerates **V×J** scaffolds for every locus — NOT the full
V×D×J the original spec mentioned for IGH/TRB/TRD.

## Why
The FR/CDR region *coordinates* we transfer at runtime are fully determined by:
- **V**: FR1, CDR1, FR2, CDR2, FR3, and the CDR3 start (conserved Cys104).
- **J**: the CDR3 end and FR4 (conserved `[FW]GXG`).

The D segment lies **inside** CDR3, whose sequence is somatic/query-specific at
runtime. Enumerating D therefore adds no markup information but multiplies the DB
~50× (human IGH: 5,306 V×J vs 254,688 V×D×J). Confirmed with the user.

For VDJ loci we still insert a short **frame-neutral N spacer** (`DEFAULT_D_SPACER_NT
= 9`, a multiple of 3) where D would sit, so IgBLAST annotates a plausible CDR3 +
FR4. See `refbuild/combinations.py`.

## Dedup
Byte-identical assembled scaffolds collapse to one DB entry; `combinations.tsv`
records every (V,J) allele pair mapping to each scaffold. Dedup yield is low
(alleles differ by SNPs) — e.g. IGH 5,306 → ~5,124 — but it's free and lets the
AIRR output report ambiguous allele calls correctly (comma-joined `v_call`).

## Scale (human, measured)
IGH ~5.1k, TRB ~2.6k, TRA ~8.6k unique scaffolds; ~25–30k total across 7 loci per
species. Watch committed-size if this grows; consider gzip if it becomes heavy.
