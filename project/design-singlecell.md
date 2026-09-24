# Single cell in arda: what is missing, and the order to build it

2026-08-14. Supersedes the two-stage sketch at `ROADMAP.md:14-30`.

## Two premises this replaces

**"arda only works in amplicon mode."** `amplicon` and `rnaseq` are **speed presets over one
pipeline**, not structure models. `src/arda/cli.py:842` `_MODE_SPEED` is four booleans --
`two_pass`, `fast_segments`, `segment_only_v`, `prefilter` -- and both modes run the identical
`map -> assemble -> correct` through one shared body, `_mode_run` (`cli.py:849`), calling
`rnaseq.pipeline.run(**speed, **kw)`. `--exact` zeroes the preset and the pipeline is unchanged.
Nothing in this document is about a mode.

**"arda needs to be made to report read ids."** It already does. `sequence_id` is the first column
of `AIRR_COLUMNS` (`annotate/transfer.py:33`), set at `transfer.py:338` and `mapper.py:1430`;
`<prefix>.airr.tsv` and `<prefix>.assembled.airr.tsv` are per-read tables keyed by it, and
`arda correct --read-map` (`cli.py:655`, `correct.py:1417`) writes the link back out. Measured: a
migec consensus named `PBMC.AAACCTGCAGCCTGTT.CGTTTTTATC` arrives in `sequence_id` intact, because
`dnaio` truncates at whitespace and the name has none. Only `<prefix>.clones.tsv` drops it.

## What is actually missing

arda has **no `cell_id` column, no per-cell table, and no chain pairing**. `cli.py:1042` says so in
the source: "arda has no barcode or UMI concept at all, which is the real single-cell gap -- not the
assembler."

The gap is narrower than that comment implies, because **migec already does the half arda's roadmap
called "the real gap"**: barcode whitelist, barcode correction, UMI consensus and cell calling. What
arda needs is not a barcode subsystem. It is:

1. a way to **lift** a cell id out of an identifier it is already carrying, and
2. the **per-cell reasoning** that follows -- chain pairing, doublets, ambient chains, per-cell
   clonotype tables.

## The migec/arda contract

migec's `assemble` writes, per molecule (`migec/src/assemble.cpp:501-511`):

```
name  <sample>.<cell>.<umi>[.c<k>][.<m>]
tags  RX:Z:<umi>  BC:Z:<sample>  CB:Z:<cell>  MI:Z:<name>  cD:i:<reads>
```

**The name is sufficient for sample, cell and UMI**, and that is why migec puts them there. Every
per-cell and per-molecule feature in this document reads the name and nothing else.

⚠ **Corrected 2026-09-24: `dnaio` does NOT drop the comment -- arda does.** This section used to
say the tags cannot reach arda because dnaio discards the FASTQ comment. Measured on dnaio 1.2.4
(arda's floor is `>=1.2`): `SequenceRecord.comment` carries the whole tag string, `.id` is the
first whitespace-delimited token and `.name` is both. What discards it is `map.py:489`, which
yields `rec.id`. So `cD:i:<reads>` -- the read depth behind each consensus, and the only tag
carrying information the name does not -- is one field away, not a sidecar TSV away. `MI:Z:` is
the name repeated and adds nothing. See S4 for what that does and does not unblock.

Note: `migec/docs/formats.rst:141` documents a colon-delimited name that the code does not write.
Parse against the code, not the page.

## S0 -- the identifier parser

Ships **first, and inside migec's `scripts/compare_cellranger.py`**, not in arda. `arda correct
--read-map` plus this parser already gives per-cell chain calls today from two existing commands,
so the measurement decides what belongs in arda rather than the other way round.

When it moves in, it lands as `src/arda/cell.py`: a pure module, no I/O, no CLI, no new dependency.

```python
@dataclass(frozen=True, slots=True)
class CellKey:
    sample: str | None
    cell: str | None
    umi: str | None
    contig: int    # 1 unless .c<k> -- ONE molecule's fragments
    molecule: int  # 1 unless .<m>  -- DIFFERENT molecules

def parse(sequence_id: str, dialect: str) -> CellKey | None
def sniff(sequence_ids: Iterable[str], *, min_match: float = 0.95) -> str
def make_parser(dialect: str = "auto", regex: str | None = None)
```

Named dialects **and** a regex escape hatch. Named, because each dialect carries a validator a
regex cannot express. Regex, because the dialect set will never be complete (BD Rhapsody, Parse,
in-house).

| dialect | rule |
|---|---|
| `cellranger` | `^(?P<cell>[ACGTN]+)(-\d+)?_contig_\d+$` |
| `migec` | parse right to left, below |
| `prefix` | `^(?P<cell>[ACGTN]+)[_.]` |

**The migec rule must parse from the right, and that is not cosmetic.** `validate_sample_id`
(`migec/include/migec/types.hpp`) forbids `/`, `\`, `..`, a leading `-` and control characters. It
does **not** forbid a `.`, so `PBMC.rep2` is a legal sample id and a left-to-right split reads
`rep2` as the barcode. The walk: strip a trailing all-digit field to `molecule`; strip a trailing
`c<digits>` field to `contig`; require the last remaining field to be `[ACGTN]+` for `umi`; the
field before it is the cell if it too is `[ACGTN]+`; everything left is the sample.

Never: **a length check, not just an alphabet check.** A bulk sample named `TCGA`, `CAG` or any
donor code drawn from those four letters otherwise walks as a fabricated single-base-alphabet cell,
`sniff` reports 100 percent matched, and every molecule in the sample lands in one invented cell.
`sniff` learns the modal cell-field and UMI-field lengths and refuses when they vary -- a real 10x
cell field is exactly 16 nt in every record.

Never: **key on sample plus cell, never cell alone.** 10x barcodes come from a fixed 737,280-entry
whitelist, so identical `cell_id` strings across two samples are guaranteed, not unlikely. Anything
that merges shard AIRRs (`pipeline.py:259-267`) hits this.

Never: **`sniff` reservoir-samples.** migec writes in barcode order, so the first 10,000 ids share
their leading bases and are not a sample of anything.

Note: arda's own `/1` and `/2` mate suffixes are stripped by `_fragment_id` (`map.py:179`) before
any of this runs; the parser must be given the fragment id, not the raw query id.

## S1 -- `cell_id` as an AIRR column

`cell_id` is AIRR Rearrangement spec-valid, so this is a column arda should have regardless of 10x.

- New constant `CELL_ID = "cell_id"` beside `JUNCTION_QUALITY` in `src/arda/rnaseq/map.py`.
- `map()` gains `cell_from: str = ""` and `cell_regex: str | None = None`.
- The hook goes in **`flush()`**, beside the `with_junction_quality` loop (`map.py:726-731`), and
  `CELL_ID` is appended to `extra_cols` **last**.

Never: **not at `mapper.py:1430`.** That line is inside the unmapped-record branch, which
`mapped_only=True` -- what `map` always passes -- `continue`s past before reaching. Every mapped
record is built by `transfer_hit` at `mapper.py:1448`. A hook there is dead code on the only path
`map` uses.

Never: **appended last, never prepended.** The comment at `map.py:686` pins "every extra goes at
the END, in a fixed order, so a consumer reading the shipped set by position is unaffected
whichever combination is on". Prepending would move `junction_quality` whenever both are on.

CLI: `--cell-from {auto,migec,cellranger,prefix,none}` (default `none`) and `--cell-regex`, on
`map`, `amplicon` and `rnaseq`. Help text hoisted to a module-level `_CELL_FROM_HELP`.

**No `arda singlecell` mode.** The obvious preset for it is the all-False vector, which is
byte-for-byte what `--exact` already produces on either existing mode (`cli.py:849-856`) -- a
no-op mode. Keep the exit-2 stub and `tests/unit/test_cli.py:91`; change only its message to point
at `arda rnaseq --exact --cell-from migec`. Add the mode when a measured speed row differs from
both existing presets.

Acceptance: `map --cell-from migec` on a migec consensus FASTQ emits a `cell_id` column whose
values are the 16 nt barcodes, `airr_header` ends with it, and a bulk run without the flag emits a
byte-identical file to today's.

## S2 -- `locus` and the clonotype index in `--read-map`

One line at `correct.py:1418`. This is the highest-value change in the document and the smallest.

`--read-map` currently writes `sequence_id -> junction`. That does not close the per-cell join,
because **`junction` alone is not the clonotype key** -- `correct.py` keys on
`(locus, v_call, j_call, junction)`. Emit `sequence_id, locus, junction, clone_row` where
`clone_row` is the 0-based row index into `<prefix>.clones.tsv`.

Acceptance: joining `--read-map` to `clones.tsv` on `clone_row` reproduces `clones.tsv`'s own
`duplicate_count` for every row.

## S3 -- per-cell chain pairing and `<prefix>.cells.tsv`

Built from the **`--read-map` join**, never from `.airr.tsv`.

Never: **the per-read junction is not the clonotype.** `correct.py:94-104` measures it truncated on
42 percent of IGH and 29 percent of TRB reads on 100 bp input; aggregating those per cell is
exactly what that comment forbids. The corrected junction only exists after `correct`.

`<prefix>.cells.tsv`, one row per `(sample_id, cell_id, chain)`:

```
sample_id  cell_id  chain  slot  v_call  j_call  c_call  junction  junction_aa
           molecules  reads  chain_fraction  productive  status
```

`status` is one of `paired`, `orphan`, `extra`, `ambient_candidate`, `doublet_candidate`.

Slot assignment is **locus-set dependent**: decide TR-ab against TR-gd against IG per cell from the
dominant loci first, then assign heavy and light within that set. Otherwise TRD is simultaneously
in the doublet-raising slot and excluded from it, and a genuine gamma-delta cell has no rule.

`doublet_candidate` is raised on the **heavy slot only** (TRB / IGH / TRD-in-a-gd-cell). Two
productive TRA chains are allelic inclusion, a known biological population, not a doublet.

`ambient_candidate` needs `--ambient-fraction`, and there is no measurement for it yet. Ship it
unset: emit `chain_fraction` and let the reader decide, rather than shipping an unmeasured default
that looks calibrated.

## S4 -- molecule counting

Rescoped 2026-09-24. The previous version parked S4 as *"blocked on a migec format decision, not
on arda"*. It was not blocked. The blocker analysis was answering a different question from the
one AIRR's field asks.

### What AIRR actually asks for

> `umi_count` -- Number of **distinct UMIs** represented by this sequence.

Distinct UMIs. Not distinct molecules, not distinct consensus records. That distinction is the
whole of the rescope, because the thing the old analysis called a blocker lives strictly *inside*
one UMI:

```
name  <sample>.<cell>.<umi>[.c<k>][.<m>]
                             │      └─ flat index over (component x split)
                             └──────── overlap component: a CONTIG of the molecule
```

`.c<k>` and `.<m>` subdivide the reads that already share one `<umi>`. A UMI counted once is
counted once whether it arrived as one record, as three contigs, or as a split pair -- so **no
reading of the conditional suffixes can change a distinct-UMI count.** `<sample>`, `<cell>` and
`<umi>` are written unconditionally, `arda.cell.parse` already returns all three, and
`tests/unit/test_cell.py:39` already pins the suffix defaults.

The old blocker is real, but it is about *molecule identity* -- joining a record to its row in
`<sample>.mig.tsv`, where the name's conditional suffixes and the table's unconditional columns
disagree (`migec/src/assemble.cpp:511-512`). Nothing in `umi_count` needs that join.

Never: **do not implement `umi_count` as distinct `(sample, cell, umi, contig, molecule)`.** That
counts *records*, which on one-record-per-molecule input is what `consensus_count` already is --
a redundant column -- and on contig-mode input over-counts a molecule once per contig, which is
the one case where `umi_count` was supposed to say something new.

### Why it is worth a column at all, measured

`umi_count` differs from `consensus_count` exactly when **several consensus records sharing one
UMI land in the same clonotype**. Measured 2026-09-24 on real migec output rather than reasoned
about, because the first version of this section reasoned about it and got it wrong.

**It is the SPLIT case, not contig mode.** `assemble_component`
(`migec/src/consensus.cpp:469-507`) cuts one overlap component into two molecules when minor
alleles co-segregate above `linkage_threshold` (8.68) -- a saturated barcode holding two
templates. When the co-segregating positions lie outside the junction, both records carry the
*same* junction and land in the same clonotype. 60 reads on one UMI split 30/30 at three
positions in V, plus 10 reads on a second UMI:

```
S.AAAACCCC.1  cD:i:30      ->  duplicate_count 3   consensus_count 3   umi_count 2
S.AAAACCCC.2  cD:i:30
S.GGGGTTTT    cD:i:10
```

**Contig mode does NOT produce it.** `place_reads` never extends a component across a gap, so
two contigs of one molecule do not overlap and **at most one of them can carry a complete
junction**; the others are dropped by `complete_only` before they ever reach a clonotype.
Measured, 12 reads on `[0,170)` and 12 on `[230,400)` of one molecule plus a second clean UMI:
migec emits `S.AAAACCCC.c1`, `S.AAAACCCC.c2`, `S.GGGGTTTT`, arda keeps two of the three, and
`consensus_count == umi_count == 2`.

So the column's value is **on saturated barcodes**, which is the condition migec's own occupancy
reporting exists to surface -- not on fragmented molecules. On an unsaturated amplicon library it
equals `consensus_count` and carries nothing new, which is why it is opt-in rather than always on.

### What ships

One column, `umi_count`, in `<prefix>.clones.tsv`.

- `correct()` gains `cell_from: str = ""` / `cell_regex: str | None = None`, the same pair
  `map()` already takes. `pipeline.run` already resolves `auto` to a concrete dialect once
  (`pipeline.py:257-259`) and passes it to `map`; it passes the same resolved string to `correct`.
  Standalone `arda correct` gets the two options so the stage is usable on its own, as every
  other stage is.
- Beside `dup` / `cons` at `correct.py:1332`, a third list: distinct `(sample, cell, umi)` over
  each clonotype's read set, with `_strip_mate` applied first for the same reason `cons` applies
  it.
- Never: **appended LAST**, after `_flag_chimeras`, not inserted into the `pl.DataFrame` literal.
  The literal is followed by two conditional appends (`_clonotype_d`, `_flag_chimeras`); a key
  added to the literal moves `d_call` whenever `umi_count` is on. Same rule, and same reason, as
  `map.py:740-743` for `cell_id`.
- Never: **omitted, never 0**, when no dialect is active or the identifiers do not parse --
  QC trap 4. A bulk run must not read as "1 UMI per clonotype".

Acceptance: on a migec consensus FASTQ, `umi_count <= consensus_count` for every clonotype, with
strict inequality on a barcode migec split into two molecules and equality where each UMI yields
one junction-bearing record; a run without `--cell-from` writes a `clones.tsv` with no
`umi_count` column at all, byte-identical to today's.

### Deliberately not in S4

**`cD:i:<reads>` -> `consensus_count`.** AIRR defines `consensus_count` as *"the sum of the number
of reads for all UMIs that contribute to the query sequence"*, and `cD:i:` is literally that
number -- so on migec input the spec-correct `consensus_count` is the summed depth, while arda
writes distinct fragment consensuses. Now reachable (the comment survives; see the contract
section above), and deliberately separate: it **changes an existing column's values** on an
existing input type rather than adding a new one, so it wants its own flag, its own measurement
and its own release note. `umi_count` changes nothing that already ships.

**Molecule identity and the `.mig.tsv` join.** Still blocked on the same migec emitter decision,
and still migec's to make. Nothing in arda needs it today.

## Out of scope, and why

| not arda's | whose | why |
|---|---|---|
| barcode demultiplexing | migec `checkout` | pattern matching on raw reads; arda never sees a raw read |
| whitelist correction | migec `refine` | needs the background hypothesis and a measured off-list prior |
| UMI consensus | migec `assemble` | needs per-base quality across a group; arda's input is one sequence per record |
| cell calling | migec `refine` (OrdMag) | a property of the molecule count distribution, not of the receptor |
| EmptyDrops rescue | Cell Ranger | stated out of scope in migec too |
| genome alignment | minimap2, bwa | arda maps to germline scaffolds, not to a genome |

## Order

S0 (in migec's comparison script) -> S1 -> S2 -> S3 -> S4. S1 and S2 are what the
migec-against-Cell-Ranger per-cell chain axis needs; S3 is what makes the result readable without
a join written by hand; S4 is the one column that separates a molecule from its contigs.

S4 used to read "waits on a migec format decision". It does not -- see its own section. What
still waits on migec is molecule identity for a `.mig.tsv` join, which nothing here needs.
