# Changelog

Notable changes per release. Earlier releases are described by their git tags
(`git tag --sort=-v:refname`); this file starts at 2.5.0.

## 2.5.0

### Behaviour changes for existing users

**`d_call` output moves.** Three independent fixes changed which records get a D and which D
they get. If you depend on 2.4.0's D calls, re-run rather than diff.

* The V..J interior is now bounded by the **per-allele junction anchors** shipped in
  `database/vdj/<org>/cdr3_anchors.tsv`, not by the mmseqs scaffold projection. A scaffold has
  a 9 nt N-pad where a read has 20–40 nt of N-D-N, so the projection collapsed the very window
  the D lives in — on real human IGH, an 11 nt interior where the truth is 37 nt. This, not
  "paralogous germlines + SHM", was the cause of the old 46–69 % IGH D concordance. IGH D-gene
  agreement with IgBLAST is now **86–98 %** across the five organisms (human 97 %, rat 98 %,
  rhesus 94 %, mouse 94 %, rabbit 86 %).
* The call is gated on a **Karlin–Altschul E-value** (`_D_MAX_EVALUE = 0.2`), replacing four
  hand-tuned per-locus score floors. The E-value ships as the AIRR `d_support` column so a
  consumer can re-threshold. Some records that got a D at 2.4.0 correctly get none now: the
  best hit simply does not clear the gate.
* **TRBD2 is never called with a TRBJ1.** TRBD2 lies 3′ of the entire TRBJ1 cluster and V(D)J
  joining deletes the intervening DNA, so the pair cannot be produced. Unenforced, TRBD2
  (16 nt) outscored TRBD1 (12 nt) on noise and took 17 % of real human TRB J1-cluster D calls.

**Output schemas gain columns** (additive; consumers reading by name are unaffected):

* `rnaseq correct` clonotype table: `d_call`, `d2_call`, `d_support`, `d2_support`.
* `rnaseq assemble` AIRR: the four above plus `np1`, `np2`, `np3`.
* `--seqtype aa` annotation now populates the D columns, which were previously always empty.

**`rnaseq correct` row order is now deterministic.** Clonotypes tied on
`(duplicate_count, consensus_count)` used to fall back to read order, which comes from a
threaded mmseqs search — the same FASTQ produced the same rows in a different sequence each
run. Ties now break on `(junction, v_call, j_call)`.

### Added

* **`arda.cdr3fix` + `arda markup`** — mark up and repair a bare `(CDR3 aa, V, J, species)`
  record with no read behind it: which residues each germline templates, where the submission
  disagrees, how far the disagreement extends, and a conservative repair. Everything is
  *junction space* (Cys104 … Phe/Trp118, both anchors included), which is what VDJdb's `cdr3`
  column holds and is **not** arda's `cdr3` field. Emits a VDJdb-compatible `cdr3fix` JSON
  blob and a human-readable `--report`.
* **`arda.annotate.dmap.map_d_junction`** — D and tandem D-D on a bare nucleotide junction,
  with no mmseqs pass.
* **`arda.dpost.posterior_d`** and `arda markup --d-posterior` — the D gene and its position
  inferred from the junction *length*, which pins `insVD + |D surviving| + insDJ`. Places the
  D to a median 1–3 nt even when the protein shows none of it. Shipped for human IGH/TRB/TRD
  and mouse TRB; every other pair returns nothing rather than guessing.
* **D mapping on amino-acid input**, against each D germline's three translated frames. A call
  on ~36 % of real IGH records, agreeing with the nucleotide call on 98 % of them; mostly
  silent for the TR loci, whose D is too short to survive trimming into protein.
* **New reference artifacts**, per organism: `cdr3_anchors.tsv` (per-allele conserved Cys104 /
  Phe-Trp118 anchors, with a `status` that is flagged, never guessed) and `d_prior.tsv`
  (generative-model summaries used by `dpost`; derived, not measured — see `SOURCES.md`).
* **`examples/`** — a runnable tour, every artifact derived from data already in the repo and
  rebuilt by `python examples/regenerate.py`: one mRNA per locus; the two human reads (of
  7,341 across five organisms) that carry a tandem D-D; seven VDJdb records covering every
  junction-repair outcome; and a 1,035-read FASTQ that runs the whole bulk RNA-seq pipeline in
  about six seconds.
* `--map-d/--no-map-d` on `rnaseq run`, `rnaseq assemble` and `rnaseq correct`.

### Fixed

* **`cdr3fix` repairs always land on a canonical junction.** `cdr3_repaired` is accepted only
  when it opens with Cys104 and closes with Phe/Trp118; otherwise the submission comes back
  untouched. Trimming flanking framework (`_MAX_TRIM = 3`) is now budgeted separately from
  inventing germline (`_MAX_FIX = 2`), and a trimmed residue costs something — a free trim
  could tie the untrimmed alignment, win the tie-break, and eat the conserved Phe of clean
  short IGK junctions. On the committed 250-row VDJdb fixture: repair reproduced **100/100**,
  zero novel rewrites, and `good`/`vCanonical`/`jCanonical` agree with VDJdb on every row.
* `Cdr3Error.applied` claimed an edit had been written to the output when a `Failed*` side had
  in fact discarded the whole repair.
* `Cdr3Markup.v_end` counted residues of the *submitted* junction, so a V-side trim or insert
  shifted it off `cdr3_repaired` — and `dpost` slices the non-templated middle with exactly
  those coordinates.
* `Cdr3Markup.v_canonical` / `j_canonical` described the submission rather than the repair,
  disagreeing with VDJdb's own booleans on 76 of 250 fixture rows.
* `dpost` reported `entropy = -0.0` for a degenerate posterior.
* `examples/example.airr.tsv` had been stale for four release rounds — written with 49 AIRR
  columns while the schema grew to 83. It is now regenerated and tested.

### Removed

* The private per-locus D score floors (`_D_MIN_SCORE`, `_D2_MIN_SCORE`), superseded by the
  E-value gate.
* The `family*01` rung of the allele-resolution ladder: reachable for **0** genes across all
  five shipped organisms.

### Docs and tests

* `arda.cdr3fix`, `arda.dpost` and `arda.annotate.dmap` are in the Sphinx API reference; the
  zero-warning gate still passes.
* Corrected two false claims: `rnaseq assemble` is implemented and on by default, and
  `rnaseq run` is three stages / four files.
* Test suite: 132 → 228 test functions, **268 passing** (1 skipped: optional `airr` dep).
