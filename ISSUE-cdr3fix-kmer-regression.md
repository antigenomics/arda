# `cdr3fix` declines repairs the k-mer fixer it replaced applied

Measured against VDJdb's curation corpus on 2026-09-30: **190,902 distinct `(species, cdr3, v, j)`
keys**, `arda-mapper` 2.31.0, `markup_cdr3` one call per key.

`arda.cdr3fix` replaced VDJdb's `Cdr3Fixer` (vendored copies of the original in
`reference-legacy-cdr3fixer/` beside this file). On one class of input it detects the defect, names it
exactly, and does not apply the repair, so a malformed junction reaches the shipped database where the
retired fixer had corrected it.

## 1. Framework past an anchor is reported, not repaired - and asymmetrically

```
YFCASSQSPGGVAFFGQG   TRBV14 / TRBJ1-1   ->  CASSQSPGGVAFFGQG
   v_fix=FixTrim             j_fix=FailedNoAlignment
   V ins@0 extra 'YF'  d=0
   J ins@13 extra 'GQG' d=0 (reported, not repaired)

YFCAVVGTGLGYTFGSG    TRBV9  / TRBJ1-2   ->  YFCAVVGTGLGYTF
   v_fix=FailedNoAlignment   j_fix=FixTrim
   V ins@0 extra 'YF'  d=0 (reported, not repaired)
   J ins@14 extra 'GSG' d=0
```

The same defect at both ends of two junctions, and the decision flips per end. Expected
`CASSQSPGGVAFF` and `CAVVGTGLGYTF`; the 2026-06-03 VDJdb release ships both, produced by the k-mer
fixer.

## 2. A junction cut inside the germline is not completed

```
CAISGEFGSGA          TRBV10-1 / TRBJ2-6  ->  CAISGEFGSGA
   v_fix=NoFixNeeded         j_fix=FailedReplace
   J del@10 missing 'NVLTF' d=0 (reported, not repaired)
```

`TRBJ2-6*01` templates `SGANVLTF`; the junction ends `SGA`, the germline's first three residues, so
five templated residues were cut. Expected `CAISGEFGSGANVLTF`, which the release ships as `FixAdd`.
`_MAX_FIX = 2` is the bound that refuses it. That bound is sound where residues are *invented* against
an ambiguous boundary, and this is not that case: the germline plainly continues and three residues of
it are still present.

## 3. `NoFixNeeded` is returned without an alignment

This is the largest count and the one with no counterpart in the legacy behaviour.

```
CASSQQQQQQQQQF       TRBV14 / TRBJ1-1    ->  v_fix=NoFixNeeded  j_fix=NoFixNeeded
CFALRGNNAGNMLTF      TRAV20 / TRAJ39     ->  v_fix=NoFixNeeded  v_end=1
```

`TRBJ1-1*01` templates `NTEAFF` and the first junction shares only its terminal Phe. `TRAV20*01`
templates `CAVQ` and `CF` shares only the Cys. Legacy returns `NoFixNeeded` **only** on a k-mer hit
with `start_in_segment == 0 and start_in_cdr3 == 0`, i.e. at least `min_hit_size = 2` residues
agreeing at the anchor (`Cdr3Fixer.fix`, lines 137-146 of the vendored copy). One residue is the
anchor coincidence, not an alignment.

**2,165 keys** are passed this way - 1,396 at the V end, 769 at the J end, every one at depth 1. They
concentrate on particular alleles rather than scattering, which suggests the calls or the anchors
rather than 2,165 independent curation errors:

| V call | templated | junction starts | n (sample) |
|---|---|---|---|
| `TRAV12-1*01` | `CVVN` | `CAAD`, `CAAG`, `CAAF` | many |
| `TRAV10*01` | `CVVS` | `CAAC`, `CAAF` | many |

Full table: `cdr3fix-nofixneeded-without-alignment.tsv` (2,165 rows, columns `species`, `submitted`,
`arda_output`, `v_call`, `j_call`, `end`, `call`, `templated`, `germline_depth`, `arda_fix`).

## Counts

382 of 190,902 keys end up with a junction that does not run Cys104 to the anchor its own segment
encodes - 267 wrong at the J end, 115 at the V end. 263 human, 100 mouse, 19 rhesus; 279 TRA, 103 TRB.

| class | n | dominant `v_fix` / `j_fix` |
|---|--:|---|
| J: anchor absent or mis-read, no germline support | 215 | `NoFixNeeded` / `FailedReplace` (122) |
| V: no Cys104 near the start | 99 | `FailedNoAlignment` / `NoFixNeeded` (72) |
| J: framework after the anchor | 35 | `NoFixNeeded` / `FailedReplace` (17) |
| J: truncated inside the germline | 17 | `NoFixNeeded` / `FailedReplace` (15) |
| V: framework before Cys104 | 13 | `FailedNoAlignment` / `FixTrim` (7) |
| V: one framework residue before Cys104 | 3 | `FailedNoAlignment` / `NoFixNeeded` (3) |

Full table: `cdr3fix-noncanonical-junctions.tsv` (382 rows, columns `species`, `submitted`,
`arda_output`, `v_call`, `j_call`, `v_templated`, `j_templated`, `v_fix`, `j_fix`, `bad_v`, `bad_j`,
`class`).

Worked examples per class, V and J calls as the curator wrote them:

```
J framework after the anchor
  CAASEGNNRLAFW          TRAV13-1       TRAJ39*01    j_fix=FailedNoAlignment
  CAASIGFGNVLHCG         TRAV12-2       TRAJ35       j_fix=FailedReplace
  CAASMGFGNVLHCG         TRAV12-2       TRAJ35       j_fix=FailedReplace
  CAATIGFGNVLHCG         TRAV12-2       TRAJ35       j_fix=FailedReplace

J truncated inside the germline
  CASSPEQGYQET           TRBV7-9        TRBJ2-5      j_fix=FailedReplace
  CASSSGLLSNTG           TRBV9          TRBJ2-2      j_fix=FailedReplace
  CASSSGQLTNTE           TRBV9          TRBJ1-1      j_fix=FailedReplace
  CASSSGQVSNTG           TRBV9          TRBJ2-2      j_fix=FailedReplace

V framework before Cys104
  FRAPCSCKDDHKLMF        TRAV8-1        TRAJ16       v_fix=FailedNoAlignment
  GLNCGGSQGNLTF          TRAV26-2       TRAJ42       v_fix=FailedNoAlignment
  HPCTKEEEMRNTF          TRAV26-2*01    TRAJ48*01    v_fix=FailedNoAlignment
  LRCETSRGKLIF           TRAV1-1*01     TRAJ23*01    v_fix=FailedNoAlignment
  LCCGGDNQGGKLIF         TRAV3*01       TRAJ23*01    v_fix=FailedNoAlignment
  LCREHFGNEKLTF          TRAV12-2*01    TRAJ48*01    v_fix=FailedNoAlignment
  WCQQPPEGPPDTEAFF       TRBV14         TRBJ1-1      v_fix=FailedNoAlignment

V no Cys104 near the start
  ASMYKGGGNEKFTF         TRAV26-2*01    TRAJ45*01    v_fix=FailedNoAlignment
  FAHTLNEGGFRIMF         TRAV9-2        TRAJ9        v_fix=FailedNoAlignment
  FALRGNNAGNMLTF         TRAV20         TRAJ39       v_fix=FailedNoAlignment
  FAPVGGWRTLQF           TRAV17*01      TRAJ24*01    v_fix=FailedNoAlignment
```

Two notes on the last class, both checked rather than assumed. Prepending the anchor is **not** the
answer: `ASMYKGGGNEKFTF` becomes `CASMYKGGGNEKFTF`, and `TRAV26-2*01` templates `CILRD`, so `CASM`
still agrees on one residue - and arda then returns `NoFixNeeded` for it, which is defect 3 again.
`CFALRGNNAGNMLTF` is accepted the same way and is not a sequence anyone would write. These records
need the segment re-called from the junction, which is what `Cdr3Fixer.guess_id` did.

`LCCGGDNQGGKLIF` and `CSSCALLWSASTQYF` carry **two** candidate Cys and no germline depth to choose
between them. Those should stay unrepaired and be tagged abnormal; they are in the table for
completeness, not as repair targets.

## Why the reference is the cause

`Anchor.templated_aa` runs Cys104 through [FW]118 inclusive and stops. A junction carrying framework
past an anchor therefore has nothing to align that framework *to*, so the boundary can only be inferred
from the junction side and the aligner declines. Every case above is a consequence of that one gap, and
the legacy fixer did not have it because it sliced the segment with flanks:

```python
reference_point = int(segment['reference_point'])
sequence = segment['sequence'][:reference_point + 4] if is_j_segment \
    else segment['sequence'][reference_point - 3:]
```

One codon of FR3 before Cys104, and the `[FW]GXG` motif after 118. With the flanks present the repair
is positional - `start_in_segment` and `start_in_cdr3` are the whole answer - and needs no alignment
scoring, no per-end floors and no tie rules.

**arda already ships what this needs.** `database/vdj/<organism>/alleles.fasta` holds full allele
sequences and `cdr3_anchors.tsv` carries `anchor_nt` per allele, which is the reference point. So:

* V scanned sequence = `translate(allele_nt[anchor_nt - 3:])`
* J scanned sequence = `translate(allele_nt[:anchor_nt + 4], from_end=True)`

## Starting point in this branch

`src/arda/cdr3kmer.py` on `feature/kmer-cdr3fix` has the scanner and the translation ported verbatim
from the legacy code, with the four-outcome table documented. **Not finished**: it does not yet build
the flank-extended segment sequences from `alleles.fasta` + `anchor_nt`, and has no `fix()` entry
point, no batch API and no tests. It is a starting point, not a working replacement.
