# arda design memory

Running notes on non-obvious design decisions and gotchas, so future work
(and future sessions) don't re-derive them. One topic per file.

- [scaffold-enumeration.md](scaffold-enumeration.md) — why V×J only (not V×D×J), dedup, N-spacer.
- [reading-frames.md](reading-frames.md) — V frame detection, J frame from aux, FR4 sanity check.
- [igblast-gotchas.md](igblast-gotchas.md) — dummy D db, IGDATA, aux files, AIRR coords.
- [markup-transfer.md](markup-transfer.md) — runtime projection, indel/strand semantics, D-segment TODO.
- [mmseqs-params.md](mmseqs-params.md) — tuned nt/aa params, all-loci single DB, caching, speedup.
- [discordance-and-scaling.md](discordance-and-scaling.md) — why arda≈IgBLAST, alignment-phase frame fix, bulk RNA-seq prefilter speed, 30M-read estimate.
- [cluster-slurm.md](cluster-slurm.md) — split/merge/slurm sharding design + how it's tested without a cluster.
