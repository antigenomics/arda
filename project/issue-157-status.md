# Issue 157 — tandem-D controls and optimized alignment

Implemented backward-compatible `map_dd=True` / `--map-dd` defaults.
`map_dd=False` / `--no-map-dd` keeps single-D mapping; `map_d=False` /
`--no-map-d` skips both. Controls propagate through annotation, assembly,
correction, single cells, pipelines, SLURM renderers and Snakemake. Nextflow
uses its existing `task.ext.args`. AIRR columns and default calls are preserved.

Profiling isolated repeated per-allele D alignment. Both greedy selections now
reuse one existing GIL-released native batch within a query. No new native
kernel, dependency, worker pool, or persistent computed-result cache was added.
One-query batches run serially inside the existing outer annotation budget.
The existing best-span-per-allele and genomic-order rules remain unchanged.

## Verification

- Full unit/synthetic/realworld suite: 1469 passed, 9 optional OLGA skips,
  207.97 seconds. Additional pipeline-finish and Snakemake flag checks pass.
- 100 seeded comparisons verify scalar/native-batch calls, ties, coordinates
  and excluded-span selection. Existing nucleotide/protein D tests pass.
- Actual CLI TSV comparison on shipped tandem examples preserves non-D fields
  across default, single-D and disabled-D modes. Clonotype finish verifies that
  correction does not re-enable tandem calls.
- Source ruff, whitespace checks and zero-warning ARDA Sphinx build pass.
- Optional vdjtools 5.9.0 integration against this ARDA source: three loci pass.
  For 1000 native draws per locus (seed 157), generated tandem histories /
  aligned tandem calls are TRB 46/0, TRD 292/50, IGH 124/4. Single-D mode has
  no D2 calls. Called spans partition the junction interior exactly.
- Companion vdjtools `model.data.annotate_reads` and `prepare` forward both
  controls. Default commands stay unchanged, retaining older ARDA compatibility.
  All four flag combinations and preparation forwarding pass (five tests).
  Companion zero-warning Sphinx build passes. Source, docs and skill are separate
  from ARDA; no dependency floor advanced before ARDA publication.

## Performance

Five repeats of 5000 IGH query interiors: median scalar-two-search time
1.067525 seconds; native-batch-once time 0.543230 seconds (1.97x D-selection
speedup). This is not a whole-annotation speedup.

Repeatable end-to-end benchmark: three fresh annotation calls per mode, median
wall seconds, 5000 raw sequencing observations per locus. Every mode preserves
all non-D/NP output fields.

| Locus | Tandem seconds | Single-D seconds | No-D seconds |
|---|---:|---:|---:|
| TRB | 1.854145 | 1.847468 | **1.833962** |
| TRD | 1.280892 | 1.332732 | **1.214757** |
| IGH | 4.662031 | 4.677154 | **4.638666** |

Run `RUN_BENCHMARK=1 python -m pytest tests/benchmark/test_tandem_controls.py -s`.
Set `ARDA_BENCH_FASTQ` to a compressed raw input for its first 5000 observations;
otherwise it uses 1000 queries from shipped tandem examples. Private raw input
provenance stays alongside temporary inputs outside Git.

## Audit and limitations

Correctness audit traced every D caller and CLI boundary: defaults, positional
compatibility, no-D precedence, protein coordinates, ambiguity, ordering and
E-value gates remain intact. Native batches retain allele order and convert
1-based spans back to the existing 0-based selection contract. No new parallel
layer or thread pool is introduced. Ponytail review found no new dependency or
speculative abstraction; existing native and scheduler helpers own the work.
An independent read-only reviewer inspected both commits, native coordinates,
defaults, downstream argv and main propagation paths, finding no actionable
correctness or overengineering issues. It did not independently rerun tests.
The final focused D/control/Snakemake run passes 67 tests in 0.51 seconds.

D alignment remains preliminary, greedy evidence, not an exhaustive history
posterior. Trimmed generated tandems often lack detectable second-D evidence;
missing D2 must not forbid tandem inference. Whole-annotation D-disable timing
gains are small on these inputs. Reported process lifetime peak RSS ranges
217–258 MB across raw-input benchmark modes; high-water marks are cumulative,
not independent mode-specific peaks. Optional OLGA tests need its extra.
No reference vocabulary, model fit, published asset, release, or remote state
was changed by this task.
