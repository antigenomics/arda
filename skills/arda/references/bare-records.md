# Bare records: a CDR3 amino acid, a V call, a J call, no read

Loaded on demand from `SKILL.md`. Markup, repair and D inference for rows that have nothing to
align -- VDJdb-style records.

## Bare records — a CDR3 amino acid, a V call, a J call, no read

VDJdb-style rows have nothing to align. The V and J germlines still template a known run of
residues into each end of the junction, and arda ships those per allele.

```python
from arda.cdr3fix import markup_cdr3, markup_records   # markup_records: a whole polars frame
from arda.annotate.dmap import map_d_junction          # D (+ tandem D-D) on a bare nt junction
from arda.dpost import posterior_d                     # D gene + position from junction LENGTH

mk = markup_cdr3("CAIRDDKII", "TRAV12-3*01", "TRAJ30*01", "HomoSapiens")
mk.cdr3_repaired             # 'CAIRDDKIIF'  -- the Phe118 anchor restored
mk.v_end, mk.j_start         # residues templated by V / index of the first J residue
[str(e) for e in mk.errors]  # ["J del@8 missing 'F' d=0"]
mk.good                      # both sides repaired AND both anchors present
mk.to_cdr3fix()              # VDJdb's `cdr3fix` JSON object, key for key
```

CLI: `arda markup -i vdjdb.txt -o marked.tsv --vdjdb --report - [--d-posterior]`.

> **The single biggest correctness trap.** These coordinates are **junction space**: Cys104
> through Phe/Trp118, **both anchors included**. That is what VDJdb's `cdr3` column holds. It
> is **not** arda's `cdr3` field, which excludes both — `junction_aa` is two residues longer
> than `cdr3_aa`. Conflating them silently corrupts every coordinate, and downstream corrupts
> Pgen, clustering and matching.

Repair is deliberately conservative and its two decisions are separate:

- Every germline disagreement is **reported** (side, kind, position, extent, distance from the
  anchor). Only edits *adjacent* to a conserved anchor are **applied**; deeper ones are left
  alone, because there a mismatch is as likely to be the real V/N boundary as a typo.
  `Cdr3Error.applied` is true only when the edit reached `cdr3_repaired`.
- **A repair always lands on a canonical junction.** If the result would not open with Cys104
  and close with Phe/Trp118, it is refused and the submission returned untouched. So `good`
  implies canonical. An allele with no derivable anchor gives `FailedBadSegment` — flagged,
  never guessed.

`posterior_d` infers the D gene *and where it sits* from the junction's nucleotide length,
which pins `insVD + |D surviving| + insDJ`. Shipped for human IGH/TRB/TRD and mouse TRB only
(the pairs with a published generative model); **every other pair returns `None` rather than
guessing** — do not substitute a human proxy.
