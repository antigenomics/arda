# Memo: implementing single cell in arda

2026-08-14, **superseded in its S3 half on 2026-08-15**. Read `design-singlecell.md` for the why;
this was the what-to-type. S0-S3 are merged; only **S4** is still a work list.

## Where things stand

| stage | state | what landed |
|---|---|---|
| S0 identifier parser | **done** | `src/arda/cell.py` — `CellKey`, `parse`, `sniff`, `make_parser`, `modal_length`. `tests/unit/test_cell.py`, 18 cases |
| S1 `cell_id` column | **done** | `--cell-from` / `--cell-regex` on `map`, `amplicon`, `rnaseq`. Hook in `map.py` `flush()`. `tests/synthetic/test_cell_id_end_to_end.py` |
| S2 `--read-map` | **done** | now `sequence_id, locus, junction, clone_row` (`correct.py:1417`) |
| S3 per-cell table | **done, differently** | `arda cells` — see below. `src/arda/singlecell.py`, `src/arda/partition.py`, `src/arda/scplot.py`, `docs/singlecell.rst`, `notebooks/singlecell_qc.py` |
| S4 `umi_count` | **blocked** | on a migec format decision, see the last section |

Verified end to end: `arda map --cell-from migec` on 29 migec-named receptor records emits
`cell_id` as the 16 nt barcode on 29/29, with `.c<k>` and `.c<k>.<m>` suffixes stripped correctly;
`clone_row` reproduces `clones.tsv`'s `duplicate_count` for every row.

### What S3 got wrong, and what shipped instead

The plan below builds the per-cell table from the **`--read-map` join over per-read junctions**.
That is not what shipped, because measuring the premise first changed the design:

- **A molecule's consensus covers one window of the transcript, not the junction.** Reads of one
  `(CB, UMI)` are co-terminal in 5' droplet chemistry, so depth does not extend a molecule. At
  `--min-reads 30` only **7,855 of 47,584** migec molecules (16.5%) carried a cell, a locus and a
  junction at all, with a mean consensus of 204 nt against a ~508 nt amplicon. A join over
  per-molecule junctions therefore throws away five sixths of the evidence before it starts.
- **But the cell's molecules TILE the transcript**, because different molecules start at different
  positions. Scored before writing any assembler: **99.90%** of the 25-mers of Cell Ranger's 943
  filtered contigs are already present in their own cell's raw molecules, and 942/943 CDR3
  nucleotide sequences appear verbatim in one of them.

So `arda cells` **assembles the cell first** — trim, 25-mer seeds, verified overlaps, union-find
layout, haplotype phasing, weighted column consensus — and annotates the contig. Everything the
memo specifies below (the columns, the slot rules, the statuses, the doublet rule on the heavy
slot) survived; only the source of the junction changed, from a per-read join to a per-cell
contig. The join columns S2 added are still the right thing for a non-assembled path.

Three things were each measured by removing them, and all three are load-bearing:
adapter trimming (CDR3 exact 0.9618 without it), depth weighting (0.967 against 0.988), and
phasing (0.9820 without it, 0.9894 with). Result against Cell Ranger: **933/943 CDR3s verbatim,
TRB 479/479, chain recall 0.9777**, in 23 s over 479 cells.

## S3 — `arda cells`, the per-cell table

### The shape of the command

A **new leaf command**, not a mode and not a flag on `correct`. It joins three files that already
exist and writes one; it maps nothing and corrects nothing, so it belongs beside `stats` and `shm`
in the "reads an artifact, writes a table" family rather than in the pipeline.

```
arda cells --airr <prefix>.airr.tsv --read-map <prefix>.read_map.tsv \
           --clones <prefix>.clones.tsv -o <prefix>.cells.tsv
           [--cell-from auto] [--cell-regex RX] [--min-molecules 1]
```

`--cell-from` is accepted here too, so a `<prefix>.airr.tsv` produced *without* `--cell-from` can
still be resolved after the fact from its `sequence_id` — otherwise the user has to re-run Stage 1
to get a column that is a pure function of a column they already have.

Wire it into `pipeline.finish` as well, so a mode run emits `<prefix>.cells.tsv` when the input
carried a cell id. Add `"cells"` to `OUTPUTS` in `src/arda/rnaseq/pipeline.py`.

### The join

```
read_map:  sequence_id -> (locus, junction, clone_row)
airr:      sequence_id -> (cell_id, v_call, j_call, c_call, productive, ...)
clones:    row index    -> the corrected clonotype
```

Group by `(sample_id, cell_id, clone_row)`. Never by `cell_id` alone: 10x barcodes come from a
fixed 737,280-entry whitelist, so the same `cell_id` in two samples is a collision by design, not
by accident. `CellKey.sample` is what supplies `sample_id`; when it is `None` (a Cell Ranger contig
id carries no sample) write `.` and say so once in the report.

Never take `v_call` / `junction` from `.airr.tsv` for the row. The per-read junction is truncated
on 42% of IGH and 29% of TRB reads at 100 bp (`correct.py:94-104`); the corrected one only exists
after `correct`. `.airr.tsv` supplies the **cell id** and nothing else that ends up in a call.

### Columns

One row per `(sample_id, cell_id, chain)`.

| column | source |
|---|---|
| `sample_id` | `CellKey.sample`, or `.` |
| `cell_id` | `CellKey.cell` |
| `locus` | `clones.tsv` |
| `slot` | `heavy` / `light`, assigned per cell — see below |
| `v_call`, `j_call`, `c_call`, `junction`, `junction_aa`, `productive` | `clones.tsv` via `clone_row` |
| `clone_row` | the join key, kept so the row is traceable |
| `molecules` | distinct `sequence_id` for this `(cell, clone_row)` — see S4 on what this counts |
| `reads` | `duplicate_count` restricted to this cell, summed |
| `chain_fraction` | this chain's `molecules` over the cell's total for its locus class |
| `rank` | 1 for the deepest chain in the slot, 2 for the next, ... |
| `status` | `paired` / `orphan` / `extra` / `doublet_candidate` |

### Slot assignment, and the trap in it

Decide the cell's **locus class first**, then assign slots inside it:

```
TR-ab: heavy = TRB,  light = TRA
TR-gd: heavy = TRD,  light = TRG
IG:    heavy = IGH,  light = IGK | IGL
```

Class is whichever set holds the most molecules in that cell. Never assign a fixed global slot to
TRD: with one table it is both "the heavy chain of a gamma-delta cell" and "a stray chain in an
alpha-beta cell", and a rule that does not branch on the class makes it both at once.

`status`:

* `paired` — the rank-1 heavy and the rank-1 light are both present and productive.
* `orphan` — one of them is missing.
* `extra` — rank >= 2 in the **light** slot. Two productive TRA chains are allelic inclusion, a
  known biological population; calling that a doublet would delete a real cell.
* `doublet_candidate` — rank >= 2 in the **heavy** slot only.

Never ship `ambient_candidate` yet. It needs a `--ambient-fraction` threshold and there is no
measurement for one. Emit `chain_fraction` and let the reader threshold it; add the call when the
PBMC_1k run gives a number.

### Files to touch

```
src/arda/cells.py          new: the join and the slot logic, pure, no I/O beyond read/write
src/arda/cli.py            new @app.command("cells"); _CELLS_HELP hoisted to module level
src/arda/rnaseq/pipeline.py  OUTPUTS["cells"], and call it from finish() when a cell id exists
docs/singlecell.rst        new page -- write it only once the command exists
tests/unit/test_cells.py   the slot and status matrix, hand-built frames, no mmseqs
tests/synthetic/test_cells_end_to_end.py   map -> correct -> cells on the committed fixture
```

### Tests that must exist

1. A TRA+TRB cell is `paired`; dropping the TRB makes it `orphan`.
2. A cell with two productive TRB is `doublet_candidate`; a cell with two productive TRA is
   `extra` on the second and still `paired` on the first.
3. A gamma-delta cell with two TRD is `doublet_candidate` — the same rule the ab class applies to
   TRB, reached through the class branch and not through a TRD special case.
4. Two samples reusing one barcode produce two rows, not one merged row.
5. A `sequence_id` whose cell does not parse is counted in a report field and does not silently
   vanish from a molecule count.
6. `arda cells` on a `<prefix>.airr.tsv` written without `--cell-from` still works via its own
   `--cell-from`.

## S4 — `umi_count`, and the decision it waits on

Never redefine `duplicate_count` or `consensus_count`. They are AIRR-spec fields, arda's bulk users
read them, and their current meanings (`correct.py:1332-1333`) are correct for bulk.

The blocker is on **migec's** side and is worth stating precisely, because it is not obvious:

* `assemble` names a molecule `<sample>[.<cell>].<umi>[.c<k>][.<m>]`.
* `.c<k>` is the overlap component -- **one molecule's fragments**.
* `.<m>` is *not* a molecule index. `assemble_group` (`migec/src/consensus.cpp:477-482`) builds one
  **flat** vector, outer loop over components and inner over the splits of each, and `.<m>` is the
  index into that flat vector.
* Both suffixes are emitted **conditionally** (`migec/src/assemble.cpp:503-505`): `.c<k>` only when
  `components > 1`, `.<m>` only when `molecules.size() > components`.

So on a plain `--contig` group the FASTQ name and `<sample>.mig.tsv` disagree about the molecule
field, and a `umi_count` keyed on the parsed name is wrong in exactly the case `--contig` exists
for. Two ways out, and **the choice belongs to migec**:

1. migec emits both suffixes unconditionally, so the name is a total function of `(component,
   molecule)`. A format change; cheap; breaks any existing parser.
2. arda defines `umi_count` as distinct `(sample, cell, umi, contig, molecule)` as parsed, accepts
   that it equals `consensus_count` on this input, and documents that a `--contig` molecule split
   across components counts once per component.

Until one is chosen, `molecules` in `cells.tsv` means **distinct `sequence_id`**, which is exact
and honest, and the column is named `molecules` rather than `umis` for that reason.

## What arda must not grow

Restated because every single-cell request eventually asks for one of these:

| ask | whose | why not arda's |
|---|---|---|
| read a barcode out of R1 | migec `checkout` | pattern matching on raw reads; arda never sees one |
| correct a barcode to a whitelist | migec `refine` | needs the background hypothesis and a measured off-list prior |
| collapse reads of a UMI | migec `assemble` | needs per-base quality across a group; arda's unit is one sequence per record |
| call cells from a knee | migec `refine` (OrdMag) | a property of the molecule-count distribution, not of the receptor |
| EmptyDrops rescue | Cell Ranger | out of scope in migec too |
| per-cell de novo V(D)J assembly | open | `annotate/contig.py` merges **existing alignments**; it does not assemble reads into a contig. Real gap, own round, not S3 |

That last row is the one to be honest about: `merge_contigs` stitches reads that are already
placed, which is not what Cell Ranger's assembler does. If per-cell assembly is wanted, it is a new
subsystem, and migec's `assemble --contig` already does the per-molecule half of it.
