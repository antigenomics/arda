# Annotation Reference

The runtime annotation API: in-memory and file-streaming entry points, parameters,
the full AIRR output schema, and performance characteristics.

**Contents:** entry points · parameters · AIRR output fields · D-segment mapping ·
performance.

## Entry points

```python
import arda
arda.annotate_sequences(sequences, seqtype="nt", organism="human", map_d=True)
```

`sequences` is an iterable of raw strings or `(id, sequence)` tuples. Returns a
`list[dict]` of AIRR records. This is the convenience wrapper; for full control
use the mapper:

```python
from arda.annotate.mapper import annotate_records, annotate_file, build_index

annotate_records(
    records,                 # list[(id, seq)]
    organism="human",
    seqtype="nt",            # "nt" | "aa"
    threads=0,               # 0 = all cores
    sensitivity=None,        # None -> tuned default (7.0)
    strand="both",           # nt only: "both" (search + re-orient) | "forward"
    map_d=True,              # map D segments (VDJ loci, nt only)
) -> list[dict]

annotate_file(
    input, output,           # FASTA/FASTQ (gz ok) -> AIRR TSV
    organism="human", seqtype="nt",
    threads=0, sensitivity=None, strand="both",
    chunk_size=50_000,       # streaming chunk -> flat memory for huge FASTQ
    map_d=True,
) -> Path
```

`annotate_file` runs a background reader thread that prefetches the next chunk
while the current one is annotated (mmseqs releases the GIL), so memory stays flat
for arbitrarily large inputs and read parsing overlaps compute.

## Parameters

- **seqtype** — `"nt"` is the more complete path (D mapping, productivity, frame
  bridging). `"aa"` returns region `*_aa` directly; coordinates are in aa space.
- **strand** — nt only. `"both"` (default) searches both strands and re-orients
  reverse-complement hits (`rev_comp="T"`); `"forward"` searches the plus strand
  only (use for germline/sense input to avoid spurious revcomp hits).
- **sensitivity** — MMseqs2 search sensitivity; default 7.0 (tuned for short
  germline-similar queries). There is **no coverage filter**, so partial reads
  still map.
- **threads** — `0` uses all cores.
- **map_d** — map D segment(s) into the V..J interior for VDJ loci; nt only.

## AIRR output fields

Column order (`arda.annotate.transfer.AIRR_COLUMNS`):

```
sequence_id, sequence, locus, v_call, d_call, d2_call, j_call, rev_comp, productive,
v_sequence_start, v_sequence_end,
d_sequence_start, d_sequence_end, d2_sequence_start, d2_sequence_end,
j_sequence_start, np1, np2, np3, junction, junction_aa,
<for each region in fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4>:
  {region}_start, {region}_end, {region}, {region}_aa
```

- All coordinates are **1-based closed**, in query space.
- `{region}` is the nucleotide (or aa, for aa input) slice; `{region}_aa` is the
  amino-acid translation (V-side regions read in the V frame; FR4 in the J frame).
- `v_sequence_end` = CDR3 start − 3 nt; `j_sequence_start` = CDR3 end + 1 (FR4 start).
- `productive` = "T" only when in-frame and free of stop codons / N-bridge.
- Round-trip invariant: `query[{r}_start-1 : {r}_end] == record[{r}]` for every
  covered region.

## D-segment mapping (VDJ loci, nt only)

D germlines are short and trimmed, so they are mapped by a gapless C++ local
alignment of every locus D allele against the V..J junction interior (not via the
scaffold DB). For IGH/TRD a second non-overlapping D is sought; the two are ordered
5'→3' as `d_call`/`d2_call` with `np1`/`np2`/`np3` between V, the D(s), and J.

## Performance

- IgBLAST work is offline (DB build); annotation is MMseqs2 + C++ projection.
- ~4–8× faster than IgBLAST at annotation, ~97–99% region concordance.
- A single search annotates a mixed bulk file across all loci at once (one
  combined reference DB per organism).
