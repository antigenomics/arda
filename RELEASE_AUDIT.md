# Release audit — 17 October 2026

The structural annotation requirements are implemented on master `2c99626` (version 2.36.2): batched query annotation with transferred query-space FR/CDR markup, actual junction anchors, shared-locus allele membership, coherent content-validated reference indexes and exact searched-index provenance. Source discovery no longer mistakes an unrelated project's database for an Arda checkout.

Verification: **159 passed** across source-root, committed-index, concurrent-index, reference-invariant and committed real-world example checks, plus shared-anchor, locus-manifest and junction-boundary/repair regressions. Sphinx builds with `-W`. `arda info` resolves the intended database, MMseqs and native markup extension. The examples reproduce their biological fields exactly: calls, coordinates, regions and D support. Their five-locus and tandem-D output checks pass. These tests use the existing batched annotation APIs; there is no outer per-query MMseqs pool.

[Issue 149](https://github.com/antigenomics/arda/issues/149) records the two pre-existing documented roadmap gaps: TRDV × TRAJ junction recall and the full productivity contract beyond current frame/stop checks. They are distinct from the repaired 2.36.2 index/root capabilities, not newly discovered structural-annotation blockers. Live multi-node sharding remains the separately documented deployment gate.

This is a focused dependency/consumer audit, not a claim that every planned feature or a complete release matrix passed. The maintainer owns tag creation and publication. No tag or release was created.
