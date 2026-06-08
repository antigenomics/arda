# Discordance analysis & scaling

## Why arda ever disagrees with IgBLAST (resolved to ~0)
On 100 real IGH mRNA, region concordance is 99.7%. Categorized
(`scripts/analyze_discordance.py`):
- **582 exact**, **18 one-residue boundary** (FR1 start / CDR3 termini — both tools
  valid, just trim the partial terminal codon differently), **0 frameshift, 0 other**.

The earlier ~1-2% "garbage" cases had ONE root cause: when the mmseqs alignment
starts **mid-codon** (e.g. q1↔t3), the projected FR1 start is 1 nt off the true
codon boundary. The first fix (stop-free `detect_coding_frame`) failed on a record
(PX613029.1) where BOTH frames were coincidentally stop-free over the short V.
**Real fix** (transfer.py): derive the query reading frame from the **alignment
phase** — project the first in-phase target codon boundary (`p = tstart +
((t0-tstart)%3)`) to the query; that query position IS `coding_start`. Exact, no
heuristic. detect_coding_frame remains only as a gap fallback.

## V/J output
Yes — `v_call`/`j_call` come from the best-hit scaffold's allele set (comma-joined
ambiguous calls). V-gene agreement with IgBLAST was 100/100.

## Prefilter / bulk RNA-seq speed (measured, 16 threads, 150nt)
mmseqs k-mer prefilter rejects non-receptor reads before alignment, so mostly-junk
input is far faster: 100% receptor ~5.7k/s, 10% ~19k/s, 1% ~25k/s
(`scripts/bench_prefilter.py`). This matches blood bulk RNA-seq (~1-5% receptor).

## CDR3 / junction correctness (AIRR-critical)
CDR3 length is query-specific, so its end is **J-anchored** (= FR4 start − 1), NOT
taken from the fixed-length scaffold — otherwise long somatic CDR3s get truncated.
`cdr3_aa == junction_aa[1:-1]` holds **by construction** for every emitted junction.

**Out-of-frame junctions are reported** (not dropped). `transfer.py _junction_nt`:
the nt junction is the real query slice (Cys104 codon .. [FW]118 codon). For
translation, when V and J are in different frames (phase = (fwr4_start −
coding_start) % 3 ≠ 0), insert k=(3−phase)%3 N **after the V germline end**
(`v_sequence_end`, clamped to stay inside CDR3 so Cys/[FW] flanks survive); the
codon(s) containing inserted N render as `_`. FR4 is translated in its own J frame
(`translate(Q[fwr4_start-1:])`) so it reads `WGQG…` even for non-productive reads.
`productive` = in-frame AND no stop in V-side/junction.

**Extended scaffold markup**: `markup.tsv` now stores `v_sequence_end` /
`j_sequence_start` (scaffold nt); these are transferred to queries (point
projection) → AIRR `v_sequence_end`, `j_sequence_start`, locating the V/J split in
the junction. V and J can rarely overlap/cross — handled by clamping.

**Concordance scoring**: non-canonical junctions (not C…[FW], or containing `_`)
are reported but **excluded** from junction/cdr3 concordance metrics. On the
committed fixtures, productive-canonical junction & cdr3 match IgBLAST ~100% per
species; productive region concordance 98–99.7%.

## Trimmed inputs & region deletion
- V-only (FR1-FR3) fragments: annotate FR1-3 + CDR1-2 correctly, no CDR3/FR4.
- J-side (CDR3-FR4) fragments: coding frame is derived from the **alignment phase**
  (not FR1), so CDR3/FR4 aa are produced even without V.
- Deleting an internal region (e.g. CDR2) collapses it to ~0 residues; flanks and
  distal regions still round-trip exactly (nt and aa). Tests in tests/synthetic.

## Committed test fixtures (offline)
`tests/assets/realworld/<organism>.fasta.gz` + `<organism>.igblast.airr.tsv.gz`:
balanced ~7.3k GenBank mRNA across all 5 organisms × loci (IG all; TR human/mouse;
NCBI lacks 500 for rare groups like TRG so totals < 10k), gzipped, with IgBLAST AIRR
reference. `tests/realworld` runs offline (mmseqs + DB only), parametrized per
organism. Rebuild via `scripts/build_test_fixtures.py`.

## Multi-species concordance (productive records)
Region concordance vs IgBLAST: human 98.6%, mouse 99.6%, rat 98.0%, rabbit 99.7%,
rhesus 99.4%. junction/cdr3 vs IgBLAST ~99% all species. KEY: compare only on
IgBLAST-productive records and skip IgBLAST regions with stops — GenBank junk
(genomic/partial/pseudogene, e.g. mouse IGK AM0863xx with stops in FR3) otherwise
drags mouse to ~91%. On productive rearrangements arda ≈ IgBLAST everywhere.

## Junction emitted only with both flanks
A junction is emitted ONLY when both conserved residues are present (Cys before
cdr3 AND [FW] opening FR4). Truncated queries past the Cys get empty junction
rather than a partial one — so `cdr3_aa == junction_aa[1:-1]` holds for every
emitted junction. (transfer.py `_set_junction`.)

## 30M-read / 32-core estimate
~25k/s at 16 threads, 1% content → ~35-40k/s at 32 cores (mmseqs scales sublinearly)
→ **~10-20 min for 30M reads**. Same order of magnitude as a STAR genome-mapping
pass (STAR ~10^5 reads/s on many-core nodes, Dobin et al. Bioinformatics 2013 —
approximate, hardware/index dependent; STAR is faster per read but indexes the whole
genome). **Caveat**: current `annotate_file` loads all reads into memory (~10-15 GB
for 30M); needs streaming/sharding (Roadmap TODO) before running at that scale.
