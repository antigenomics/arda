# Personalized germline — `arda genotype`, `arda resolve-ties --genotype`

A reference catalogues every allele **anyone** carries. No donor carries all of them, and none
carries more than two per gene — so a call naming a third is wrong before any sequence is looked
at, and a read that "could not choose" between two alleles the donor does not have was never
ambiguous. TIgGER measured the size of this on full-length BCR: restricting V calls to an inferred
genotype took ambiguous assignments from 11.2 % to 1.5 %.

Two halves with very different confidence, and they are separate commands on purpose.

## Applying a genotype — solid, and useful with zero inference

```bash
arda resolve-ties -i sample.airr.tsv -o sample.genotyped.tsv --genotype genotype.tsv
```

The simplest valid genotype file is a one-column TSV of allele names (`allele` header, then
`TRBV19*01`, …). An OGRDB set, a library MiXCR inferred, a list typed by hand — all first-class.

- Adds **`v_call_genotyped`**; `v_call` stays byte-identical, so a reader can re-restrict at a
  different genotype without re-running anything.
- Never: it **re-assigns, it never re-aligns and never rebuilds a reference.** Given the span a
  read already aligned over, the restriction is a set intersection. Three arda-specific reasons a
  per-donor reference is a trap: scaffold ids are positional (`f"{locus}_{idx}"`), so changing the
  allele set renumbers every scaffold in the locus; `build-db` needs IgBLAST plus IMGT network
  access a `pip install` does not have; and the mmseqs freshness contract is an mtime with **no
  allele-set identity recorded**, so a donor-specific FASTA beside the shared one silently
  invalidates the index for every concurrent process. TIgGER's `reassignAlleles` reaches the same
  conclusion for its own reasons.
- Never: **an allele survives unless its OWN gene was genotyped and did not name it.** Tie lists
  routinely span genes, so a gene-blind test empties every read whose tie list merely brushes a
  genotyped gene — measured on the 453-read fixture with a one-gene genotype, 79 rows came back
  empty and every one was a read of some other gene.
- Never: **an unknown allele raises.** 884 human V alleles are functional in IMGT and **801** reach
  a scaffold, so naming one of the other 83 is an easy mistake that would otherwise restrict
  nothing, silently.

`v_call_genotyped` has three value classes and they are deliberately distinguishable: **narrower**
than `v_call` (the point), **equal** to it (the tie machinery declined, or no candidate's gene was
genotyped — a refusal to answer is not a contradiction), **empty** (every candidate's gene was
genotyped and none is carried — this read contradicts the genotype).

## Inferring one — a likelihood ratio, not a coverage rule

```bash
arda genotype -i sample.airr.tsv -o sample.genotype.tsv --loci TRB
```

**The junction is de facto a UMI.** Junctional diversity makes a nucleotide junction essentially
unique to one rearrangement, so grouping on `(locus, gene, J gene, junction)` groups by molecule.
Two things follow and the inference rests on both:

- **The clonotype, not the read, is the unit of observation** — one junction is one independent
  draw from the donor's two chromosomes however many reads carry it. Read depth is reported
  *beside* the clonotype count (`clonotypes`, `reads`), never substituted for it.
- **Disagreement within a junction is error, not allele** — every read of one rearrangement carries
  the same V allele by construction, so a read naming another is SHM or sequencing error. That is
  where the error rate comes from, measured on the library instead of picked.
  ⚠ Measured over **all** reads, not the germline-exact subset used for assignment: those were
  *selected* for carrying no mismatch, so they never disagree and the estimate collapses silently
  onto its floor (0 discordant of 77,345 reads on the amplicon below).

Every single allele and every pair is then scored by its multinomial likelihood under that rate,
and the winner is called only if it beats the runner-up by `--min-log10-bf` (default 1.0).
Homozygous and heterozygous are the same formula at genotype size 1 and 2.

⚠ **The first version was TIgGER's frequency rule and had to be replaced.** "Fewest alleles
explaining 7/8 of the calls" has no error model, so it cannot tell overwhelming evidence from none:
on a real TRB amplicon it called `TRBV11-2` off **754 of 757** clonotypes and `TRBV20-1` off **1 of
2,544**, reported both as `explained = 1.0000, ok`, and returned 43 of 53 genes with none
heterozygous.

Ambiguous clonotypes still constrain the answer — discarding them is wrong, not conservative. If a
donor is `*01/*02` and only `*02` is separable at the read length, every unambiguous observation is
`*02`, the likelihood says homozygous, and the restriction then empties every `*01` read on every
sample with that genotype. `--fraction-to-explain` is the guard against exactly that.

`note` values: `ok` · `single_allele` (the reference has one allele for this gene) · `low_support`
(< `--min-clonotypes` assignable) · `undetermined` (best genotype does not beat the runner-up) ·
`unexplained` (fits the counts, leaves too many ambiguous clonotypes unexplained). A gene that
could not be called gets a row with an empty `allele` — **omitted with a reason, never absent.**

## ⚠ Read length is the binding constraint, not the rule

Separating a gene's alleles, measured from the 3' end of the V germline over all 801 committed
human V germlines:

| | multi-allele genes | median nt needed | separable within 100 nt of the 3' end |
|---|---|---|---|
| TRBV | 44 / 56 | 150 | 16 / 44 |
| TRAV | 30 / 45 | 175 | 8 / 30 |
| IGHV | 57 / 75 | 230 | 9 / 57 |

Measured end to end on `SRR5233641` (human TRB amplicon, 151 nt paired, 99,839 mapped reads →
33,440 clonotypes / 77,345 reads, error rate 5.25 × 10⁻⁴): reads cover a median of **72 nt** of V
germline — **56 nt** after clipping at Cys104 — and **not one** reaches 150. Only **33.1 %** of
clonotypes can be assigned an allele at all. Result: **17 of 53 genes called**, 11 of them
`single_allele` and **6 genuinely inferred** (log₁₀ BF 223 for `TRBV11-2` over 754 clonotypes, 251
for `TRBV5-6` over 846, 11.7 for `TRBV5-8` over 39); 36 refused, including `TRBV10-3` with **1,037
clonotypes** that still cannot separate its alleles.

The refusals are the output, not a failure of it. Judging the rule needs a library with the
resolution — full-length, 5'RACE, or `arda cells` contigs at N50 536 nt.

**What it looks like when a gene IS separable.** TRA amplicon `SRR5233635`, 100,000 reads, same
151 nt (21,710 clonotypes / 45,007 reads, error rate 5.53 × 10⁻⁴ — two independent libraries, same
order): **20 of 44 genes called and one heterozygous** — `TRAV36/DV7` = **\*01 / \*04**, 78
clonotypes against 143, **log₁₀ BF 211** over the gene's 290. Five more inferred homozygous at
log₁₀ BF 53–221, 14 `single_allele`, **24 refused** `low_support`. The whole inference is **1.53 s
/ 442 MB** on 49,748 mapped reads — about a tenth of what `map` cost.

⚠ **Applying it narrowed 281 of 47,743 rows** (75 contradicted, 47,387 unchanged). You cannot
restrict what you could not genotype, and 14 of the 20 called genes have one allele to begin with.
The 11.2 % → 1.5 % TIgGER reports is a full-length-library number.

Never: **`narrowed` counts rows that lost an ALLELE, not rows whose STRING changed.** The first run
of that restriction reported 20,587 narrowed; 20,306 of them were `TRAV20*02,TRAV20*01` becoming
`TRAV20*01,TRAV20*02` — the tie list comes back sorted by name while `v_call` carries the aligner's
order. Survivors now keep `v_call`'s order, so a no-op restriction is byte-identical.

## What this is not

- **A claim about the reference, never about the repertoire.** Clonal composition, diversity and
  overlap are `vdjtools`'.
- **Not** `arda stats`' `allele_candidate`, which stays a shortlist of recurrent high-quality V
  mutations to look at and is deliberately never a call. The genotyper does not read it.
- **No novel alleles.** Restricted to what the reference already catalogues. TIgGER's y-intercept
  regression needs mutated reads (which TCR does not supply) and IMGT numbering (which arda has
  nowhere — `refbuild/imgt.py` actively ungaps).

## Where the germlines come from

`arda.germline.segment_germlines(organism, "v"|"j")` derives per-allele germlines from the
**committed reference** (`markup.tsv` + `alleles.fasta`), not from the IMGT source tree — which is
absent on every `pip install`, and is why `arda resolve-ties` used to raise there.

- Never: **`markup.tsv`'s `v_call` is not always one allele.** Scaffolds are deduplicated by
  assembled sequence, so alleles whose scaffolds come out byte-identical collapse into one
  comma-joined row — 23 such groups hiding 49 alleles in human V. Keyed on the group string they
  are invisible to every consumer.
- Never: **take the MODE of `v_sequence_end` across scaffolds, not first- or last-wins.** It is
  IgBLAST's markup of one assembled scaffold, not a property of the allele, and it wobbles — mouse
  disagrees on 4 of 897 V alleles, always by 1–2 nt with a landslide majority (`TRAV16*02`: 288 nt
  on 58 scaffolds, 290 nt on 1). Ties break by longer, then by sequence: a total order.
- 49 human V alleles in 23 groups have byte-identical germlines, so **no read length separates
  them** — two of the groups are TR. A genotyper has to carry such a group whole rather than pick a
  member.
