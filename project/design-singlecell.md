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

Never: **the tags do not reach arda.** They live in the FASTQ comment and `dnaio` drops it. Only
the name arrives. That is sufficient for sample, cell and UMI, and it is why they are in the name at
all. It is **not** sufficient for `cD:i:<reads>` -- the read depth behind each consensus. If a
depth-weighted clonotype count is ever wanted, the depth must travel as a sidecar TSV keyed on the
molecule id, not by teaching arda to parse SAM tags.

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

Add a new **`umi_count`**. Never redefine `duplicate_count` or `consensus_count`: they are
AIRR-spec fields with defined meanings, and arda's bulk users read them.

Never: **`umi_count` cannot be reconstructed from the molecule name.** `.<m>` is a **flat index over
(component x split)**, not a molecule index -- `migec/src/consensus.cpp:477-482` builds one flat
vector, outer loop over overlap components and inner over the splits of each. And each suffix is
emitted **conditionally** (`assemble.cpp:503-505`): `.c<k>` only when `components > 1`, `.<m>` only
when `molecules.size() > components`. So on a plain `--contig` group the FASTQ name and
`<sample>.mig.tsv` disagree about the molecule field, and a join keyed on it refuses correct input.

Two ways out, and the choice belongs to migec, not arda: define `umi_count` as distinct
`(sample, cell, umi, contig, molecule)` and accept that it equals `consensus_count`, or make
migec's emitter write both suffixes unconditionally first. Until one is decided, S4 does not ship.

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

S0 (in migec's comparison script) -> S1 -> S2 -> S3. S4 waits on a migec format decision.
S1 and S2 are what the migec-against-Cell-Ranger per-cell chain axis needs; S3 is what makes the
result readable without a join written by hand.
