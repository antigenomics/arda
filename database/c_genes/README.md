# Constant-region (C) genes — CH1 exons

The IMGT V-QUEST reference directory that `arda refbuild` downloads ships **only V/D/J**. These files
supply the missing constant region.

## Why it matters

A V-J scaffold ends at the J segment's 3′ terminus, so a read running from the J across the splice
junction into the constant region has **nowhere legal to align**. That class is ~22 % of all real
receptor fragments in bulk RNA-seq, and a `V + pad + J` reference recovers only about half of it.

## How they are used — additively, not by appending C to every scaffold

**A J→C read has no V.** That is the definition of the class, so C never has to sit downstream of a V.

`refbuild` therefore builds a small, separate set of `J + C[:150]` scaffolds and adds them *beside* the
V-J scaffolds (see `refbuild/constant.py`). For human that is **345 extra scaffolds, +2.0 %**. A read
spanning V→J→C still maps to its V-J scaffold and soft-clips the C tail; a V-less J→C read finally has
somewhere legal to end.

The alternative — appending `C[:150]` to *every* V-J scaffold — is accuracy-identical and far more
expensive: it multiplies each locus's scaffolds by its number of distinct C stubs (11 for IGH), taking
human from 17,244 scaffolds / 6 Mnt to 73,360 / 37 Mnt. Measured, that raises MMseqs2's per-query cost
1.71× and peak RSS by 367 MB, for no accuracy gain. Do not reintroduce it.

`markup.tsv` records `vj_end` for every scaffold: its full length for a V-J scaffold, the J length for a
`J + C` scaffold. An alignment with `tstart >= vj_end` lies wholly inside the constant region — real
receptor mRNA, but carrying no V(D)J and therefore no clonotype.

## Provenance

Each record is the **CH1 exon** of one constant-region gene, lifted from the reference genome assembly
by exon coordinates.

- CH1 is the exon that splices directly onto J. It begins **mid-codon** — the codon straddles the J–C
  splice — so `J + CH1` reconstructs the mRNA and translates contiguously. That property is the
  bundle's acceptance test: `IGHJ4 + IGHG1` → `WGQGTLVTVSS|ASTKGPSVFP`, `IGHJ4 + IGHM` →
  `WGQGTLVTVSS|GSASAPTLFP`, `IGKJ1 + IGKC` → `WTFGQGTKVEIK|RTVAAPSVFI`. `tests/unit/test_rnaseq.py`
  asserts all three; a bundle that fails them is rejected rather than written.
- Alleles are named `*01`. The genomic source carries no IMGT allele table, so the suffix is nominal.
- Minus-strand genes are stored with **descending** coordinates. A naive slice silently loses the first
  base — `GAACTGTGG…` where the truth is `CGAACTGTGG…` for IGKC — and every downstream translation comes
  out as garbage. Extract by coordinate direction; never hand-slice.

Report the isotype **class**, never the subclass: IGHG1–4 are ~95 % identical over CH1, so the
top-scoring *gene* ties on ~27 % of real reads while the top *class* is unique on essentially all of
them (`refbuild.constant.isotype_class`).

Gene names follow IMGT's locus-prefixed convention (`TRGC1`, not `TCRG-C1`). A name whose first three
characters are not a locus is dropped **silently**, taking that locus's whole constant region with it —
`tests/unit/test_rnaseq.py` asserts the mapping is total for every shipped bundle.

| organism | CH1 exons | loci with J+C scaffolds | J+C scaffolds |
|---|---|---|---|
| human | 25 | 7 | 345 |
| mouse | 21 | 7 | 210 |
| rabbit | 15 | 3 | 137 |
| rhesus_monkey | 8 | 3 | 153 |
| rat | 3 | 0 | 0 |

**rat has no J+C scaffolds, and that is correct.** Its three CH1 exons are all TR (`TRAC1`, `TRBC1/2`),
and the IMGT V-QUEST reference directory ships **no TR directory at all** for rat — so there is no J to
splice them onto. `refbuild` logs the skip per locus rather than failing the species build.
