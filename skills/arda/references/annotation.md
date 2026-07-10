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
    map_d=True,              # map D segments (VDJ loci; works on aa input too)
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

## D-segment mapping (VDJ loci)

D germlines are short and trimmed, so they are mapped by a gapless C++ local
alignment of every locus D allele against the V..J junction interior (not via the
scaffold DB). For IGH/TRB/TRD a second non-overlapping D is sought; the two are
ordered 5'→3' as `d_call`/`d2_call` with `np1`/`np2`/`np3` between V, the D(s), and
J. `d_call`/`d2_call` are comma-joined lists when alleles tie (7 pairs of human IGH
D germlines are byte-identical across different genes).

Three things constrain the call:

- **The interior comes from the anchors, not the projection.** A scaffold has a 9 nt
  N-pad where a read has 20–40 nt of N-D-N, so mmseqs parks those bases against the
  flanking V/J and the projected interior collapses (real IGH: 37 nt of truth, 11 nt
  projected). `transfer._anchored_vj_bounds` recovers the bounds from the per-allele
  germlines in `cdr3_anchors.tsv` — longest common prefix/suffix of the junction.
- **The gate is a Karlin–Altschul E-value**, `_D_MAX_EVALUE = 0.2` (`d_support` in the
  output), replacing four hand-tuned per-locus score floors. λ = ln((1−p)/p) with p the
  chance two residues match: 1/4 → ln 3 for nt; a *measured* 0.0613 → 2.7285 for aa.
- **Genomic order forbids some (D, J) pairs.** TRBD2 lies 3′ of the whole TRBJ1 cluster
  and V(D)J joining deletes the intervening DNA, so TRBD2 × TRBJ1 cannot be produced.
  Unenforced, TRBD2 (16 nt) outscored TRBD1 (12 nt) on noise and took 17 % of real human
  TRB J1-cluster D calls. `transfer._allowed_d` masks the candidate set while holding the
  E-value's `n` at the full locus size, so a J1 record can lose an impossible call but
  never gain a weak one. IGH and TRD put every D 5′ of every J: nothing is masked there.

**aa input** gets the same call against each D germline's three translated frames (a
trimmed D has no knowable frame). Informative for IGH — a D on ~36 % of real records,
98 % gene agreement with the nt call — and mostly silent for the TR loci. `d_germline_*`
and `d_cigar` are withheld there: the offsets index a reading frame, not the germline.

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

De-novo assembly lives in `arda.rnaseq.assemble` (Stage 3, anchored greedy
overlap-extension); it calls `reannotate_contigs` and carries the contig's D call onto
every member read. The two paths above only *annotate* a contig, whoever assembled it.

## Bare records: junction markup, repair, and D without a read

`arda.cdr3fix` marks up a `(junction_aa, v_call, j_call, species)` record — a VDJdb row —
against `cdr3_anchors.tsv`: which residues each germline templates, where the submitted
junction disagrees, and how far. Everything is **junction space** (Cys104..Phe/Trp118,
both anchors included), which is what VDJdb's `cdr3` column holds and is *not* arda's
`cdr3` field. Repair applies only anchor-adjacent edits (`_MAX_REPLACE`) and reports the
rest; on 102,990 VDJdb records it reproduces VDJdb's own repair on 96.4 % of those it
flags. CLI: `arda markup`.

`arda.annotate.dmap.map_d_junction` maps D (and D-D) into a bare nucleotide junction with
no mmseqs pass — the anchors give the interior directly. `arda.dpost.posterior_d` infers
the D gene and its position from the junction *length*: the length pins
`insVD + |D surviving| + insDJ`, so the D is placed to a median 1–3 nt even when the
protein shows none of it. Priors ship in `d_prior.tsv` for human IGH/TRB/TRD and mouse
TRB only; other pairs return `None` rather than guessing.

## Performance

- IgBLAST work is offline (DB build); annotation is MMseqs2 + C++ projection.
- ~4–8× faster than IgBLAST at annotation, ~97–99% region concordance.
- A single search annotates a mixed bulk file across all loci at once (one
  combined reference DB per organism).
