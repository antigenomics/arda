# Regions, Junction & Bare Germline Segments

How arda places region coordinates, the junction/CDR3 convention, and how to get
per-allele germline FR/CDR subsequences from isolated V or J segments.

**Contents:** region projection · bare germline V/J annotation · junction/CDR3 ·
amino-acid frames.

## Region projection (no coverage filter)

`transfer_hit` projects the reference scaffold's region coordinates onto the query
through the MMseqs2 alignment (`qaln/taln/qstart/tstart`). Regions whose projected
start falls outside the query coverage (`qs < 0`) are simply omitted. Consequences:

- Mutations don't move coordinates; an insertion of L nt inside CDR3 shifts FR4
  downstream by L and preserves the `[FW]GXG` motif.
- Truncated, partial, and reverse-strand queries annotate correctly (revcomp hits
  are re-oriented; `rev_comp="T"`).

## Bare germline V or J segments

Because there is no coverage filter, a lone germline segment maps to its scaffold
and returns only the regions it covers — no synthetic rearrangement needed:

```python
from arda.annotate.mapper import annotate_records

recs = annotate_records(
    [("TRBV9*01", v_germline_nt), ("TRBJ2-7*01", j_germline_nt)],
    organism="human", seqtype="nt", strand="forward", map_d=False,
)
# V record -> fwr1, cdr1, fwr2, cdr2, fwr3  (+ v_sequence_end = CDR3 start)
# J record -> fwr4                          (+ j_sequence_start = CDR3 end / FR4 start)
```

Notes:
- V-side regions (FR1–FR3, CDR1–CDR2) are V-gene-determined; J-side (FR4, and the
  J contribution to CDR3) is J-gene-determined.
- Use `strand="forward"` for germline (sense) input.
- For a J, the CDR3-contributing residues are those 5' of `fwr4_start` (in the J
  reading frame, i.e. starting at `(fwr4_start - 1) % 3`).
- Sanity-check that the returned `locus` matches the expected locus; treat
  TRA/TRD as equivalent for the dual `TRAV.../DV` genes.

This is exactly how mirpy bakes per-allele FR/CDR subsequences into its gene
library; see `tests/synthetic/test_germline_segments.py`.

## Junction & CDR3

- `junction` (nt) spans the Cys104 codon through the [FW]118 codon that opens FR4;
  `junction_aa` starts with C and ends with F/W for a canonical rearrangement.
- `cdr3` is **J-anchored**: from the V-anchored CDR3 start up to just before the
  [FW] that opens FR4 (so somatic length is query-specific).
- Out-of-frame junctions (V and J in different frames): 1–2 N bases are inserted
  after the V germline end to restore the J frame; the codon holding an inserted N
  is rendered `_`. FR4 still reads.

## Amino-acid frames

- For nt input, V-side region `*_aa` are translated in the V coding frame (derived
  from the alignment phase — works even without FR1); `fwr4_aa` reads in its own
  (J) frame regardless of productivity.
- For aa input, regions are already amino acids; no frame bridging is applied.
- Helpers in `arda.refbuild.translate`: `translate(nt, frame=0)`,
  `detect_coding_frame(nt)`, `reverse_complement(nt)`, `back_translate(aa)`,
  `aa_coords_from_nt(nt_start, nt_end, coding_start)`.
