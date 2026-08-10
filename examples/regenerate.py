#!/usr/bin/env python
"""Rebuild every committed artifact under ``examples/``.

# Run from the repo root, in the `arda` conda env:
#     python examples/regenerate.py
#
# Everything here is DERIVED from data already committed to this repo -- no network, no
# hand-typed sequences. The inputs are the five real GenBank records in `example.fasta`
# and the realworld GenBank fixtures under `tests/assets/realworld/`. Provenance for both
# is recorded in SOURCES.md.
#
# `tests/realworld/test_examples.py` re-runs this logic and asserts the committed outputs
# still reproduce. The old `example.airr.tsv` went stale for four release rounds because
# nothing checked it: it was written with 49 AIRR columns when the schema had grown to 83.
#
# 2026-07-10

Artifacts, and what each one is for:

* ``example.airr.tsv``       -- one real mRNA per locus, annotated. The 30-second demo.
* ``dd.fasta`` / ``dd.airr.tsv`` -- the only two human reads in the realworld fixtures that
  carry a *tandem D-D* rearrangement. Rare and worth seeing.
* ``junctions.tsv`` / ``junctions.markup.tsv`` / ``junctions.report.txt`` -- real VDJdb
  records covering every repair outcome, including one arda reports and refuses to rewrite.
* ``rnaseq/reads.fq.gz`` / ``rnaseq/clones.tsv`` -- a 1712-read single-end FASTQ tiled from
  real mRNA, and the clonotype table ``arda rnaseq`` produces from it.
"""

from __future__ import annotations

import gzip
import io
import subprocess
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / "examples"
FIXTURES = ROOT / "tests" / "assets" / "realworld"

# The two human reads with a tandem D-D. Found by scanning all five realworld fixtures for a
# non-empty d2_call: human 3 (one of which is a duplicate of another), rabbit 2, rhesus 1,
# and -- a real finding, not a tooling gap -- mouse 0, rat 0.
DD_IDS = ["AM408133.1", "PX612894.1"]

# Real VDJdb records, one per repair outcome. `cdr3` here is the SUBMITTED junction
# (VDJdb's `cdr3_old`); VDJdb's own repaired value is in the comment.
JUNCTIONS = [
    # id                     cdr3                 v              j              species
    ("clean-trb",            "CASSPLGQAYEQYF",    "TRBV5-1*01",  "TRBJ2-7*01",  "HomoSapiens"),
    ("v-anchor-sub",         "FLVGPQGSSASKIIF",   "TRAV4*01",    "TRAJ3*01",    "HomoSapiens"),
    ("j-anchor-missing-F",   "CAIRDDKII",         "TRAV12-3*01", "TRAJ30*01",   "HomoSapiens"),
    ("j-anchor-extra-G",     "CATSSPGLASDEQFFG",  "TRBV15*01",   "TRBJ2-1*01",  "HomoSapiens"),
    ("reported-not-repaired", "CASSSPLLSSDTQYFG", "TRBV7-2*01",  "TRBJ2-3*01",  "HomoSapiens"),
    # Framework context on both flanks: `YF` before Cys104, `GAG` past Phe118. Both trimmed.
    ("flanking-fr3-trimmed", "YFCASPGGIQYFGAG",   "TRBV14*01",   "TRBJ2-4*01",  "HomoSapiens"),
    # A real IMGT allele that arda ships with no usable anchor (`status = no_anchor`). The
    # gene's *01 marks this junction up cleanly; the *02 is flagged, never guessed.
    ("bad-segment-refused",  "CAVRSMDSNYQLIW",    "TRAV1-2*02",  "TRAJ33*01",   "HomoSapiens"),
]

# rnaseq: tile 150 bp windows across real mRNA. Deterministic -- no RNG, fixed stride, and
# records chosen by sorted sequence_id. The loci come from the committed IgBLAST annotation of
# the same fixture, so the sample spans a D locus: tiling "the first N long records" instead
# picks TRA/IGK/IGL by file order and the clonotype table comes back with no `d_call` at all.
READ_LEN, STRIDE = 150, 4
RNASEQ_LOCI = {"IGH": 3, "TRB": 3, "TRD": 2, "TRA": 2, "IGK": 1, "IGL": 1}


def _sh(*args: str) -> None:
    print("  $", " ".join(args))
    subprocess.run(args, check=True, cwd=ROOT)


def annotate_examples() -> None:
    print("example.airr.tsv")
    _sh(sys.executable, "-m", "arda.cli", "annotate",
        "-i", "examples/example.fasta", "-o", "examples/example.airr.tsv",
        "--organism", "human", "--threads", "4")


def build_dd() -> None:
    print("dd.fasta + dd.airr.tsv")
    from arda.annotate.io import read_sequences

    want = set(DD_IDS)
    found = {sid: seq for sid, seq in read_sequences(FIXTURES / "human.fasta.gz") if sid in want}
    missing = want - set(found)
    if missing:
        raise SystemExit(f"D-D records missing from the fixture: {sorted(missing)}")
    with open(EX / "dd.fasta", "w") as fh:
        for sid in DD_IDS:                       # stable order
            fh.write(f">{sid}\n{found[sid]}\n")
    _sh(sys.executable, "-m", "arda.cli", "annotate",
        "-i", "examples/dd.fasta", "-o", "examples/dd.airr.tsv",
        "--organism", "human", "--threads", "4")


def build_junctions() -> None:
    print("junctions.tsv + junctions.markup.tsv + junctions.report.txt")
    pl.DataFrame(
        {"id": [r[0] for r in JUNCTIONS], "cdr3": [r[1] for r in JUNCTIONS],
         "v": [r[2] for r in JUNCTIONS], "j": [r[3] for r in JUNCTIONS],
         "species": [r[4] for r in JUNCTIONS]}
    ).write_csv(EX / "junctions.tsv", separator="\t")
    _sh(sys.executable, "-m", "arda.cli", "markup",
        "-i", "examples/junctions.tsv", "-o", "examples/junctions.markup.tsv",
        "--id-col", "id", "--d-posterior",
        "--report", "examples/junctions.report.txt", "--show-ok")


def _rnaseq_source_ids() -> list[str]:
    """Pick source mRNA per locus: deterministic, and known to carry a complete junction.

    GenBank holds genomic, partial and non-productive entries. Taking "the first N per locus"
    by accession lands on 65-84 nt fragments, RefSeq records with no V(D)J, and ESTs -- 94 %
    of the tiled reads then fail to map. So annotate first and select on the result:
    productive, canonical junction, long enough to tile.

    Within a locus, records that themselves carry a ``d_call`` sort first, so the example
    actually demonstrates D mapping. It is a deliberately *illustrative* sample, not a random
    one -- at the shipped E-value gate an ordinary TRB junction is a no-call more often than
    not (its interior is 11-21 nt and scores 5-6 against a min_score of 7).
    Ties break on sequence_id, so the choice is stable across runs.
    """
    from arda.annotate.io import read_sequences
    from arda.annotate.mapper import annotate_records

    recs = [r for r in read_sequences(FIXTURES / "human.fasta.gz") if len(r[1]) >= 400]
    by_locus: dict[str, list[tuple[int, str]]] = {}
    for r in annotate_records(recs, "human", "nt", threads=4):
        ja = r.get("junction_aa") or ""
        if (r.get("locus") in RNASEQ_LOCI and r.get("productive") == "T"
                and ja.startswith("C") and ja.endswith(("F", "W")) and len(ja) >= 10):
            rank = 0 if r.get("d_call") else 1        # D-bearing records first
            by_locus.setdefault(r["locus"], []).append((rank, r["sequence_id"]))
    chosen: list[str] = []
    for locus, n in RNASEQ_LOCI.items():
        picked = [sid for _, sid in sorted(by_locus.get(locus, []))][:n]
        if len(picked) < n:
            print(f"  warning: only {len(picked)}/{n} usable {locus} records")
        chosen.extend(picked)
    return chosen


def build_rnaseq() -> None:
    print("rnaseq/reads.fq.gz + rnaseq/clones.tsv")
    from arda.annotate.io import read_sequences

    (EX / "rnaseq").mkdir(exist_ok=True)
    want = set(_rnaseq_source_ids())
    seqs = {sid: s for sid, s in read_sequences(FIXTURES / "human.fasta.gz") if sid in want}
    reads: list[tuple[str, str]] = []
    for sid in sorted(seqs):                       # stable order
        seq = seqs[sid]
        if len(seq) < READ_LEN:
            continue
        for k, p in enumerate(range(0, len(seq) - READ_LEN + 1, STRIDE)):
            reads.append((f"{sid}_{k}", seq[p : p + READ_LEN]))

    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:   # mtime=0 => byte-stable
        for sid, seq in reads:
            gz.write(f"@{sid}\n{seq}\n+\n{'I' * len(seq)}\n".encode())
    (EX / "rnaseq" / "reads.fq.gz").write_bytes(buf.getvalue())
    print(f"  {len(reads)} reads, {len(buf.getvalue())} bytes gzipped")

    # `--exact` so the committed artifact is the shipped one-pass output, not the mode preset's:
    # the example is a REGRESSION fixture for the annotation, and --prefilter costs ~0.15 % of
    # mapped reads, which would make the committed table depend on a speed knob.
    _sh(sys.executable, "-m", "arda.cli", "rnaseq", "--exact",
        "--r1", "examples/rnaseq/reads.fq.gz", "-p", "ex", "-d", "examples/rnaseq",
        "--organism", "human", "--threads", "4")
    # The run report carries wall time and peak RSS -- not reproducible, so don't commit it.
    (EX / "rnaseq" / "ex.arda.json").unlink(missing_ok=True)
    for stale in ("ex.airr.tsv", "ex.assembled.airr.tsv"):       # large, and derivable
        (EX / "rnaseq" / stale).unlink(missing_ok=True)
    (EX / "rnaseq" / "ex.clones.tsv").rename(EX / "rnaseq" / "clones.tsv")


if __name__ == "__main__":
    annotate_examples()
    build_dd()
    build_junctions()
    build_rnaseq()
    print("\ndone. `git diff --stat examples/` shows what moved.")
