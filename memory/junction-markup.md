# Junction markup + repair (`arda.cdr3fix`)

Marks up a record that has **no read behind it**: a CDR3 amino acid, a V call, a J call, a
species. The VDJdb case.

## 1. The single biggest correctness trap: junction ≠ CDR3

VDJdb's `cdr3` column is the AIRR **junction** — Cys104 → Phe/Trp118, **both anchors
included**. arda's own `cdr3` field *excludes* both, and its `junction` includes them.
`junction_aa` is therefore two residues longer than `cdr3_aa`.

Everything in `cdr3fix`, `dmap` and `dpost` works in **junction space**. Conflating the two
silently corrupts every coordinate, and downstream corrupts Pgen, clustering and matching.
Assert it in tests; the invariant is `cdr3_aa == junction_aa[1:-1]`.

## 2. Anchors are derived per allele, and `markup.tsv` is not a valid source

`markup.tsv` is deduped, comma-joined and scaffold-gated, and already contains a P allele.
VDJdb cites pseudogene and truncated alleles. Derive per allele instead
(`database/vdj/<org>/cdr3_anchors.tsv`), keep every functionality (F / ORF / P), and **flag,
never guess** (`status` = `ok` / `truncated` / `no_anchor`).

Three rules that look right and are wrong:

* **"IMGT position 104 == gapped nt 310..312"** is false for mouse and rhesus — V-QUEST
  carries insertion positions. It silently produced 671 wrong rhesus anchors (first residue
  `Y`, not `C`). Use IgBLAST's `internal_data/<org>/<org>.ndm.imgt` col 11 (1-based FWR3
  stop) minus 3, with a conserved-FR3-aromatic motif fallback. 1183 agree, 0 disagree.
* **"take the last Cys"** has a 3.8 % error rate: the V CDR3 tail can hold a second Cys
  (`YYC-AC-DT` in TRDV2*01, `YYCC…` in IGLV2-11*01). Take the **5′-most** Cys with an
  aromatic at i−2 and a tail of 1..14 residues → 0 disagreements.
* **The aux `frame` column is unreliable.** For `TRAJ31*01` it contradicts its own
  `cdr3_stop` and yields `QCQTHV`. Frame is `anchor % 3` by construction. Now `NNNARLMF`.

V-anchor-not-Cys after all this: **0 / 3448**.

## 3. Detection and repair are different decisions

A mismatch inside the templated window is either a curation error or simply the V/N boundary
(exonuclease trims the germline, so `templated_aa` is an upper bound, not a promise).

Repairing everything corrupted 84/3000 VDJdb records: a single mismatch needs only two
flanking matches to outscore stopping (~1/400 by chance), so `CASSPRRY-N-L-QFF` was rewritten
to `…NEQFF` against TRBJ2-1's `SYNEQFF`. Fix: `_MAX_REPLACE = 1` — only anchor-adjacent edits
are **applied**; everything deeper is **reported** with `Cdr3Error.applied = False`.

Measured on 102,990 VDJdb records: `v_end` 98.66 %, `j_start` 96.52 %, idempotent 99.99 %,
reproduces VDJdb's own repair on 96.4 % of the 5,029 it flags. ~33k rec/s, pure Python.

## 4. The aligner

Semi-global Needleman–Wunsch anchored at the conserved residue, free end gaps toward the N
region. `match=1, mismatch=-1, gap=-2`.

* **Leading germline gaps are free** (`s[i][0] = 0`, `s[0][j] = j*gap`). Without this a
  truncated submission is unrepairable — `CASSRGSVRLGTTDPQ` missing its `YF` scored 0.
  94.1 % → 96.4 %.
* **The tie-break must be `max((s[i][j], j, -i))`** — prefer consuming query, then less
  germline. The obvious `max((s, i+j, i, j))` made the aligner *prepend* `CA` to
  `CYVPGDRGGYTDKLIF`, which already began with the conserved Cys. False repairs 0.63 % →
  0.35 %, under-reports 1 → 0.

**This bug was invisible to every VDJdb concordance check and only surfaced against OLGA
generative ground truth**, which knows the true `delV`/`delJ`. That is what the synthetic
tier is for.

## 5. Ground truth harness

`choose_random_recomb_events()` (not `gen_rnd_prod_CDR3()`, which discards the D and the
trims) returns `{V,D,J,delV,delDl,delDr,delJ,insVD,insDJ}`. Note `delV` can be **negative** —
palindromic P-nt *extend* the templated region. OLGA's germlines are not IMGT's (its
`TRBV3-1*01` anchor is 267 vs IMGT's 270), so filter on `germline_matches_imgt`.

`$ARDA_VDJREARM` supplies human TRD (with D) and TRG, which OLGA lacks — together, all 7
human chains plus mouse TRA/TRB.

Invariants that hold on 3,180 generated junctions: **never under-reports**; boundaries within
1 residue 98.3 %; false repairs 0.35 %; an injected anchor-adjacent typo is detected *and*
repaired 93.7 % of the time; a deep typo is reported 99.7 % of the time and **rewritten 0
times**.
