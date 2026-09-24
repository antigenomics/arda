# Read groups, sample sheets and cluster adapters

Loaded on demand from `SKILL.md`. Everything here is about *scheduling* arda, not about what it
computes: the answer is byte-identical however the work is divided.

## Samples split across files (read groups)

A sample delivered as one FASTQ per Illumina lane, or per in-house chunk, is **one repertoire**
and gets **one** clonotype table. `--r1` / `--r2` / `--id` are repeatable and matched BY POSITION;
repeat an id to merge. `--samples sheet.tsv|csv` takes nf-core's `sample, fastq_1, fastq_2`
columns, repeated `sample` merging in row order.

```python
from arda.samples import load, read_sheet   # -> [Sample(id, pairs=((r1, r2), ...))]
samples = load(r1=[...], r2=[...], ids=["A", "A", "B"])   # two samples, three read groups
samples = read_sheet("sheet.tsv")
```

```bash
arda rnaseq -d out/ --r1 a1_1.fq --r2 a1_2.fq --id A --r1 a2_1.fq --r2 a2_2.fq --id A
arda rnaseq --samples sheet.tsv -d out/
```

**A multi-file sample is an ALREADY-SHARDED sample.** `pipeline.run` takes read groups, maps each,
concatenates in DECLARED order, and calls the same `finish` a one-file run calls. Nothing below
Stage 1 changes — not the clonotype key, not the contig layout, not the AIRR columns. Verified
byte-identical on `tests/data/rnaseq_real` cut into four.

**Never infer the grouping from filenames.** `RNA-SAMPLE_ID:12:00XX919:3_1.fastq.gz` has no field a
rule can call the sample; guessing wrong splits one repertoire into four with no error. More than
one `--r1` without `--id` is refused.

**Never sort read groups.** `A_L010` sorts before `A_L002` and the clonotype fold is not
permutation-invariant. Declared order = command-line order = sheet row order.

**Never split one sample into two to parallelise it.** Contigs that tile across the split are never
built: 57 reads as one sample against 41 + 17 = 58 as two on the test fixture. (With `--no-assemble`
they do partition exactly, 35 + 16 = 51.)

## Read groups are the unit of parallel work, samples are not

Several samples in one CLI call run SEQUENTIALLY, each with every core — mmseqs threads internally,
so N samples at cores/N is slower than N in a row. Scale out at read-group granularity instead:
`arda cluster plan` writes `readgroups.tsv` + `samples.tsv` (one `arda map` per row, then one
`arda cluster reduce` per sample) and prints both commands FILLED IN;
`arda cluster submit-samples` renders them as two SLURM arrays, map over read groups then reduce
over samples gated with `afterok`. 5 samples × 4 lanes = 20 jobs, not 5. Snakemake
(`integrations/snakemake/arda/`) and Nextflow schedule the same way.

**Never hand-assemble the flags for a split `map` + `reduce`.** `arda rnaseq` defaults
`--ec-mode rnaseq`, which reads a column Stage 1 writes only under `--junction-quality`; the mode
body wires that itself and two separate jobs cannot. Ask `arda.cluster.regime_flags(regime)` — it
returns both halves and is what `plan`, both submit renderers and both workflow integrations use.
`arda cluster reduce` also takes `--ec-mode` now; without it a sharded run denoises with `fast`
while the mode command used `rnaseq`, and nothing says so.

## Cluster: two adapters, and picking the wrong one fails silently

| input | command | shard unit |
|---|---|---|
| amplicon / single-end FASTA | `arda cluster submit-fasta` (`cluster split-fasta` + `cluster merge`) | one record |
| bulk paired RNA-seq, one pair | `arda cluster submit` (`arda cluster split` + `cluster reduce`) | one read **pair** |
| any sheet, one or many samples | `arda cluster submit-samples` / `cluster plan` (no split — the files are the shards) | one read **group** |

**Never point `arda cluster split-fasta` / `submit-fasta` at paired RNA-seq.** They write FASTA — dropping the
quality strings `--reconstruct` needs — and round-robin *records*, which puts a fragment's two
mates in different shards. There is no error; the numbers just come out wrong.

**Never shard Stage 2 or Stage 3.** `correct` counts distinct fragments and collapses error
variants globally; `assemble` grows contigs across reads. Per shard, a clone split across N
shards is counted N times and the long-CDR3 contigs Stage 3 exists for are never built, because
the reads that tile them never meet. `arda cluster submit` distributes only `map` and runs the
rest once, through the same function the mode commands use — so a sharded run is
**byte-identical** to a single-node one (verified on real data, all three artifacts).
