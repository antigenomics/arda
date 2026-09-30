# Bare records: a CDR3 amino acid, a V call, a J call, no read

Loaded on demand from `SKILL.md`. Markup, repair and D inference for rows that have nothing to
align -- VDJdb-style records.

## Bare records — a CDR3 amino acid, a V call, a J call, no read

VDJdb-style rows have nothing to align. The V and J germlines still template a known run of
residues into each end of the junction, and arda ships those per allele.

```python
from arda.cdr3fix import markup_cdr3, markup_records   # markup_records: a whole polars frame
from arda.annotate.dmap import map_d_junction          # D (+ tandem D-D) on a bare nt junction
from vdjtools.model import posterior_d_batch           # D gene + position from junction LENGTH

mk = markup_cdr3("CAIRDDKII", "TRAV12-3*01", "TRAJ30*01", "HomoSapiens")
mk.cdr3_repaired             # 'CAIRDDKIIF'  -- the Phe118 anchor restored
mk.v_end, mk.j_start         # residues templated by V / index of the first J residue
[str(e) for e in mk.errors]  # ["J del@8 missing 'F' d=0"]
mk.v_flags, mk.j_flags       # ('ok',), ('add',)   -- the per-side verdict
mk.good                      # neither side `impossible` AND both anchors present
mk.to_cdr3fix()              # VDJdb's `cdr3fix` JSON object, key for key
```

CLI: `arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report - [--d-posterior]`.

## 2.32.0 — the engine is a PORT of VDJdb's `Cdr3Fixer`, not an alignment of its own

The 2.16.0 semi-global Needleman-Wunsch is gone. One **gapless local alignment** per side
(`_markup.d_local_align`, the same C++ the D caller uses) places the germline's templated run
anywhere in the junction, and legacy's positional table decides the outcome from the two offsets.

| `start_in_segment` | `start_in_cdr3` | outcome |
|---|---|---|
| 0 | 0 | `NoFixNeeded` |
| 0 | > 0 | `FixTrim` — framework outside the junction |
| > 0 | 0 | `FixAdd` — germline the submission cut |
| > 0 | ≤ `max_replace` | `FixReplace` — substitute at the anchor |
| > 0 | > `max_replace` | `FailedReplace` |

⚠ **What legacy scanned is `templated_aa`, and the flank slicing is why that is confusing.**
`_load_segments_data` slices `sequence[reference_point - 3:]` / `[:reference_point + 4]` — but
VDJdb's `reference_point` is offset differently from arda's `anchor_nt`, and the sequence those
slices yield IS the templated run. The release proves it by trimming both ends *completely*:
`YFCASSQSPGGVAFFGQG` → `CASSQSPGGVAFF`, `vEnd=5`, `jStart=10`. A segment carrying a flanking
residue cannot do that — it hands back `FCASSQSPGGVAFFG`.

⛔ **Searching every offset is the point; do not re-anchor at index 0.** `CAMYLCASSLFGSPLHF`
against `TRBV9` (`CASSV`) has a spurious `CA` at offset 0 and the real `CASS` at offset 5. The NW
engine scored the first and called the record clean; the release trims five residues (`vEnd=4`).
⚠ But **the anchored placement wins a TIE** — `CASSQQQQQQQQQF` on TRBJ1-1 scores 1 either way, and
taking the shifted one appended an F to a junction that already ended in one.

⛔ **`templated` is the LEADING EXACT run**, which is what `fix_both` passed through as `vEnd` and
`len - jStart`. Counting a mismatched residue as templated moved the boundary off the release on
2–4 % of records. Restored residues never count: `CAISGEFGSGA` → `CAISGEFGSGANVLTF` reports
`jStart = 13`, the surviving `SGA`.

⛔ **Two bounds 2.16.0 added on top of `max_replace` are GONE, both verified against the release.**
`_MAX_TRIM = 3` refused `CAMYLCASSLFGSPLHF`, which legacy trims; `_MAX_FIX = 2` was a second knob
on one decision, so the outcome depended on which bound bit first. `max_replace` is the whole
budget, and `_canonicalise` is what protects the anchor.

### Measured against the authoritative 2026-06-03 release

189,596 distinct `(species, junction, V, J)` curation keys; **184,765 join (97.45 %)** — coverage
before any rate, because an inner join gives each side its own denominator.

| | 2.31.0 (NW) | 2.32.0 (port) |
|---|---|---|
| repaired junction agrees | 180,374 (97.6235 %) | **181,892 (98.4451 %)** |
| `vEnd` agrees | 98.2013 % | **98.3571 %** |
| `jStart` agrees | 97.8316 % | **99.7242 %** |
| the release's own repairs reproduced | 584 / 779 | **623 / 779** |
| repairs the release did NOT make | 4,196 | **2,717** |
| `good` beside an unrepaired disagreement | 19,436 | **0** |
| throughput, one process | 12,900 keys/s | **20,400 keys/s** |

## The verdict is a SET of flags per side, not one label

`v_flags` / `j_flags` carry any combination of seven names; the TSV writes them comma-joined.
Counts over the corpus above.

| flag | meaning | V | J |
|---|---|--:|--:|
| `ok` | the junction agrees with this side's germline | 174,426 | 172,909 |
| `allele` | markup used an allele the submission did not name | 6,334 | 1,545 |
| `sub` | residue(s) substituted to germline | 1,466 | 1,598 |
| `add` | residue(s) restored from germline | 6 | 304 |
| `trim` | residue(s) removed — framework outside the junction | 47 | 197 |
| `shallow` | agrees on the conserved anchor alone — call uncorroborated | 1,425 | 579 |
| `impossible` | no segment, no alignment, **or** a declined disagreement | 7,109 | 13,279 |

**`good` == neither side `impossible`, plus both anchors present.** Read it off the flags.

⛔ **`v_fix`/`j_fix` are VDJdb's worst-wins names, kept only so the `cdr3fix` JSON stays key-for-key
comparable.** One label cannot say both "I trimmed a flank" and "I found a substitution I will not
touch", and that is how a declined repair came to read as a clean record: the fix type was computed
from the *applied* edits alone, so `CAISGEFGSGA` reported `V sub@2 I>S d=2` and still returned
`NoFixNeeded` and `good`. **19,436 keys did that.**
⚠ And the name must stay the APPLIED fix — `markup_cdr3` gates the repair on `v_fix in _GOOD`, so
renaming a side that really did repair something throws the repair away.

⚠ **`shallow` keeps `good`; it is a caveat, not a failure.** `CGGSARSGELFF` against TRBV9 (`CASSV`)
agrees on the Cys alone and that is the RIGHT answer — the V is trimmed back to Cys104 and `GGS` is
N region, pinned since 2.16.0. `CASSQQQQQQQQQF` on TRBJ1-1 agrees on the Phe alone and looks like a
mis-call. **No property of a junction separates the two**, so both are flagged and both keep their
boundary. ⛔ Never turn the depth floor into a refusal — it gates only repairs that INVENT residues
(`CALRPA` on TRAJ17, which the release also refuses).

## ⛔ A contradicted call changes the ALLELE, not the sequence

`guess_allele` (legacy's `guess_id`, rebuilt on the same scan) runs before any substitution is kept:
if a different allele of the locus explains the submitted residues as they stand, the call is
re-assigned and the junction is left alone, with the `allele` flag recording it.

This exists because the alternative corrupts data. Of the anchor-adjacent substitutions 2.31.0 wrote
into an **already-canonical** junction, **74.7 % were records whose own 3' end matched a different J
allele better than the called one — against 1.1 % of untouched records**, a 68× enrichment.
`CASSLRGAATDTQYF` is a clean TRBJ2-3 junction called TRBJ2-1, which 2.31.0 rewrote to
`CASSLRGAATDTQFF` — a string no germline supports. It is now re-called to TRBJ2-3, unchanged.

⚠ **Legacy's own V guesser never worked**: `guess_id` puts `return ""` INSIDE the five-prime loop,
so it tried one prefix length and gave up — 3 non-empty V guesses in 4,000 sequences against 3,797
for J, whose branch has the statement correctly in a `for ... else`.

`posterior_d` infers the D gene *and where it sits* from the junction's nucleotide length,
which pins `insVD + |D surviving| + insDJ`. Shipped for human IGH/TRB/TRD and mouse TRB only
(the pairs with a published generative model); **every other pair returns `None` rather than
guessing** — do not substitute a human proxy.

```python
posterior_d(junction_aa, v_call, j_call, "rhesus_monkey", prior_path="fitted.tsv")
```

`prior_path` (CLI: `arda markup --d-prior PATH`, which implies `--d-posterior`) scores against a
table `arda scenarios` fitted instead of the shipped `database/vdj/<org>/d_prior.tsv`. That is the
only way to reach the **11 of 13 shipped (organism, D-locus) pairs that have no prior at all** --
OLGA has no model for them -- without overwriting a file inside the installed database, which
would silently change every later run on the machine. Using an estimate is not adopting one.
⚠ A `prior_path` that does not exist **raises**: the shipped table is allowed to be missing, a
path the caller typed is not.
