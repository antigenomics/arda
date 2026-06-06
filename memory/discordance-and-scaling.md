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

## 30M-read / 32-core estimate
~25k/s at 16 threads, 1% content → ~35-40k/s at 32 cores (mmseqs scales sublinearly)
→ **~10-20 min for 30M reads**. Same order of magnitude as a STAR genome-mapping
pass (STAR ~10^5 reads/s on many-core nodes, Dobin et al. Bioinformatics 2013 —
approximate, hardware/index dependent; STAR is faster per read but indexes the whole
genome). **Caveat**: current `annotate_file` loads all reads into memory (~10-15 GB
for 30M); needs streaming/sharding (Roadmap TODO) before running at that scale.
