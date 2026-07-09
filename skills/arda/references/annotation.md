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
- **strand** — nt only. `"both"` (default) searches both strands; a reverse-strand hit
  is annotated on the coding strand and flagged `rev_comp="T"` (its `sequence` stays as
  submitted, see below). `"forward"` searches the plus strand only (use for
  germline/sense input to avoid spurious revcomp hits).
- **sensitivity** — MMseqs2 search sensitivity; default 7.0 (tuned for short
  germline-similar queries). There is **no coverage filter**, so partial reads
  still map.
- **threads** — `0` uses all cores.
- **map_d** — map D segment(s) into the V..J interior for VDJ loci; nt only.

## AIRR output fields

The output is a **spec-valid AIRR Rearrangement** file: every record passes
`airr.schema.RearrangementSchema.validate_row` (all 14 required fields present and
typed). Column order (`arda.annotate.transfer.AIRR_COLUMNS`):

```
sequence_id, sequence, locus, v_call, d_call, d2_call, j_call, c_call, c_class,
mmseqs2_score, mmseqs2_evalue, mmseqs2_identity,
mmseqs2_{qstart,qend,qlen,tstart,tend,tlen,t_vend,t_jstart,t_vjend},
rev_comp, productive, stop_codon, vj_in_frame, v_identity,
sequence_alignment, germline_alignment,
v_cigar, j_cigar, c_cigar,
v_germline_start, v_germline_end, j_germline_start, j_germline_end,
v_sequence_start, v_sequence_end,
d_sequence_start, d_sequence_end, d2_sequence_start, d2_sequence_end,
d_germline_start, d_germline_end, d_cigar, d2_germline_start, d2_germline_end, d2_cigar,
j_sequence_start, np1, np2, np3, junction, junction_aa,
<for each region in fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4>:
  {region}_start, {region}_end, {region}, {region}_aa
```

- `sequence` holds the read **as submitted**; the annotation (coords, CIGARs, alignment
  strings) is computed on the coding strand. A reverse-strand hit sets `rev_comp="T"`,
  which per AIRR means all output data are on the **reverse complement** of `sequence` —
  `sequence` itself is not re-oriented (earlier builds stored the reverse complement here).
- All coordinates are **1-based closed**. `*_sequence_start/end` are in query space;
  `*_germline_start/end` are in the germline allele; `{region}_*` are in query space.
- `sequence_alignment`/`germline_alignment` are the aligned query / germline strings
  (the scaffold's non-templated stretch reads as `N`, per AIRR).
- `v_cigar`/`j_cigar`/`c_cigar`/`d_cigar` follow the AIRR CIGAR spec: leading `S`
  (query 5′ offset) then `N` (germline 5′ offset), an `M`/`I`/`D` body, a trailing `S`.
- `{region}` is the nucleotide (or aa, for aa input) slice; `{region}_aa` is the
  amino-acid translation (V-side regions read in the V frame; FR4 in the J frame).
- `productive` = "T" only when in-frame and free of stop codons / N-bridge;
  `stop_codon` and `vj_in_frame` surface the two facts it collapses.
- The `mmseqs2_*` columns carry the scaffold hit's alignment score and geometry;
  `t_vend`/`t_jstart`/`t_vjend` are the scaffold's V-end / J-start / V-J-end, which
  tell V-J hits from `J + C` constant-region hits (`tstart ≥ t_vjend` ⇒ wholly in C).
- Round-trip invariant: `query[{r}_start-1 : {r}_end] == record[{r}]` for every
  covered region.

## D-segment mapping (VDJ loci, nt only)

D germlines are short and trimmed, so they are mapped by a gapless C++ local
alignment of every locus D allele against the V..J junction interior (not via the
scaffold DB). For IGH/TRB/TRD a second non-overlapping D is sought; the two are
ordered 5'→3' as `d_call`/`d2_call` with `np1`/`np2`/`np3` between V, the D(s), and
J. `d_call`/`d2_call` are comma-joined lists when alleles tie (7 pairs of human IGH
D germlines are byte-identical across different genes).

## Constant region & isotype

The reference includes `J + C` scaffolds — the CH1 exon of each constant gene
spliced onto each J — so a read spanning the J→C splice (no V, hence no junction)
still maps and yields a `c_call` (CH1 exon allele) and a `c_class` isotype. **Report
the class, never the subclass:** IGHG1–4 are ~95% identical over CH1, so the top
gene ties often while the top class is unique; `c_class` collapses to `IGHG`/`IGHM`/
`IGHA`/…, or the locus constant (`IGHC`) when calls straddle classes. In paired bulk
RNA-seq the isotype of a CDR3-bearing read is recovered from its constant-region
mate (`arda rnaseq map`).

## Contig annotation

`arda.annotate.contig` gives a caller-supplied assembled contig valid
`v_cigar`/`j_cigar`/`c_cigar` two ways, which produce the **same record** (both feed a
synthetic `hit` into `transfer_hit`, so output is field-for-field comparable):

- `reannotate_contigs(records, ...)` — treat each contig as one long query and re-align
  it through `annotate_records` (a fresh mmseqs pass, then `segment_cigars`).
- `merge_contig(contig, reference, ...)` / `merge_contigs(...)` — stitch the member
  reads' existing scaffold alignments into the contig's via the C++
  `_markup.merge_alignment`, skipping the alignment pass; wins at ~10⁵ contigs/sample
  (scRNA-seq).

arda does not yet assemble contigs de-novo — these paths only re-annotate contigs the
caller already assembled.

## Performance

- IgBLAST work is offline (DB build); annotation is MMseqs2 + C++ projection.
- ~4–8× faster than IgBLAST at annotation, ~97–99% region concordance.
- A single search annotates a mixed bulk file across all loci at once (one
  combined reference DB per organism).
