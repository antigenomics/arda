"""Build a SEGMENT reference: V, J and J+C as separate targets, not a V×J product.

The shipped reference enumerates every V×J combination — 15,069 scaffolds from 775 V alleles and
124 J alleles. That costs, and it costs twice:

* **Speed.** A read covering only V aligns against every scaffold carrying that V — a median of
  13, and **67 for TRA**. Measured on a TRA amplicon: 277 gapped alignments per hitting read,
  one kept. Aligning the same reads against 1,244 segment targets instead is **6.9× faster**
  (TRA), 7.7× (TRB), 2.5× (bulk RNA-seq).
* **Accuracy.** 81 % of hitting amplicon reads sit *at* the `--max-seqs 300` cap, so the true
  scaffold is sometimes never even a candidate; and the V call is decided by a whole-scaffold bit
  score whose J half is arbitrary. Scored against IgBLAST truth on TRA, the segment reference
  takes V-gene concordance 99.00 % → **99.99 %** (95 errors → 1) and J-gene 98.47 % → **99.90 %**.

So this is not a speed/accuracy trade — the product reference was losing on both.

Derived from the *built* reference (`markup.tsv` + `alleles.fasta`), not from IMGT, so it needs
no download and is reproducible from any checkout that can already map.

Coordinates carry over almost for free, which is why this is cheap to build correctly:

* a **V** target is `scaffold[:v_sequence_end]`, and fwr1/cdr1/fwr2/cdr2/fwr3 are already
  scaffold-relative, so they transfer unchanged. `cdr3` is truncated at the V end (the read only
  ever sees the V-side stub of the junction here).
* a **J** target is `scaffold[j_sequence_start-1:vj_end]`, so every coordinate shifts by
  `j_sequence_start - 1`. It carries fwr4 and the J-side stub of `cdr3`.
* **J+C** targets are copied through untouched — they already have no V, and
  `RefEntry.is_jc` / the `t0 > 0` guard in `transfer` already handle V-less targets. That path
  ships today with 345 scaffolds, so it is exercised, not speculative.

⚠ **What this does NOT do.** A read that spans the junction hits a V target and a J target
separately, and something must merge the two into one AIRR record. That is
**85.6 % of mapped amplicon reads but only 7.2 % of bulk RNA-seq reads**, which is why the
RNA-seq path is nearly free and the amplicon path is not. The merge lives in the mapper, not
here; this module only builds the targets.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..annotate.reference import REGIONS

__all__ = ["build_segment_reference", "SegmentStats"]

_START_COLS = [f"{r}_start" for r in REGIONS]
_END_COLS = [f"{r}_end" for r in REGIONS]


@dataclass
class SegmentStats:
    """What the build produced, for the report and for tests to assert on."""

    v_targets: int = 0
    j_targets: int = 0
    jc_targets: int = 0
    source_scaffolds: int = 0

    @property
    def total(self) -> int:
        return self.v_targets + self.j_targets + self.jc_targets

    @property
    def reduction(self) -> float:
        return self.source_scaffolds / self.total if self.total else 0.0

    def as_dict(self) -> dict:
        return {**self.__dict__, "total": self.total, "reduction": round(self.reduction, 2)}


def _read_fasta(path: Path) -> dict[str, str]:
    seqs, sid, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if sid:
        seqs[sid] = "".join(buf)
    return seqs


def _shift(v: int, by: int) -> int:
    """Shift a 1-based closed coordinate, preserving the -1 'region absent' sentinel."""
    return v if v < 0 else v - by


def build_segment_reference(organism: str = "human", *, out_dir: Path | None = None) -> SegmentStats:
    """Write ``segments.fasta`` + ``segments.markup.tsv`` beside the reference.

    One target per distinct (segment, allele): the *longest* scaffold instance of that allele is
    used as the donor, so a V target is never accidentally truncated by whichever J it happened
    to be paired with.

    Returns:
        :class:`SegmentStats` — target counts and the reduction factor vs the V×J reference.
    """
    import polars as pl

    from ..paths import vdj_dir

    base = Path(out_dir) if out_dir else vdj_dir(organism)
    markup = pl.read_csv(base / "markup.tsv", separator="\t", infer_schema_length=0)
    seqs = _read_fasta(base / "alleles.fasta")

    def _i(row, key) -> int:
        try:
            return int(row.get(key) or 0)
        except (TypeError, ValueError):
            return 0

    # allele -> best (longest) donor row, so a V is cut from the scaffold that shows most of it
    v_best: dict[str, tuple[int, dict]] = {}
    j_best: dict[str, tuple[int, dict]] = {}
    jc_rows: list[dict] = []
    n_src = 0

    for row in markup.iter_rows(named=True):
        sid = row["scaffold_id"]
        seq = seqs.get(sid)
        if not seq:
            continue
        n_src += 1
        if not (row.get("v_call") or "") and (row.get("c_call") or ""):
            jc_rows.append(row)                       # J+C scaffold: already V-less, copy as is
            continue
        vend, jstart = _i(row, "v_sequence_end"), _i(row, "j_sequence_start")
        v_allele, j_allele = row.get("v_call") or "", row.get("j_call") or ""
        if v_allele and vend > 0:
            prev = v_best.get(v_allele)
            if prev is None or vend > prev[0]:
                v_best[v_allele] = (vend, row)
        if j_allele and jstart > 0:
            jlen = _i(row, "vj_end") - jstart + 1
            prev = j_best.get(j_allele)
            if prev is None or jlen > prev[0]:
                j_best[j_allele] = (jlen, row)

    out_fa = base / "segments.fasta"
    out_tsv = base / "segments.markup.tsv"
    cols = ["scaffold_id", "locus", "v_call", "j_call", "productive",
            "v_sequence_end", "j_sequence_start", "junction", "junction_aa",
            *[c for pair in zip(_START_COLS, _END_COLS) for c in pair],
            "c_call", "vj_end", "segment"]
    rows: list[dict] = []
    stats = SegmentStats(source_scaffolds=n_src)

    with open(out_fa, "w") as fa:
        for allele, (vend, row) in sorted(v_best.items()):
            sid = f"V|{allele}"
            fa.write(f">{sid}\n{seqs[row['scaffold_id']][:vend]}\n")
            rec = {c: "" for c in cols}
            rec.update(scaffold_id=sid, locus=row["locus"], v_call=allele, j_call="",
                       productive="", v_sequence_end=str(vend), j_sequence_start="0",
                       c_call="", vj_end=str(vend), segment="V")
            for r, sc, ec in zip(REGIONS, _START_COLS, _END_COLS):
                s, e = _i(row, sc), _i(row, ec)
                # keep only regions wholly inside V; clip the junction at the V end
                if s <= 0 or s > vend:
                    rec[sc], rec[ec] = "-1", "-1"
                else:
                    rec[sc], rec[ec] = str(s), str(min(e, vend))
            rows.append(rec)
            stats.v_targets += 1

        for allele, (_, row) in sorted(j_best.items()):
            jstart, vjend = _i(row, "j_sequence_start"), _i(row, "vj_end")
            sid = f"J|{allele}"
            fa.write(f">{sid}\n{seqs[row['scaffold_id']][jstart - 1:vjend]}\n")
            off = jstart - 1
            rec = {c: "" for c in cols}
            rec.update(scaffold_id=sid, locus=row["locus"], v_call="", j_call=allele,
                       productive="", v_sequence_end="0", j_sequence_start="1",
                       c_call="", vj_end=str(vjend - off), segment="J")
            for r, sc, ec in zip(REGIONS, _START_COLS, _END_COLS):
                s, e = _i(row, sc), _i(row, ec)
                if s <= 0 or e < jstart:
                    rec[sc], rec[ec] = "-1", "-1"
                else:
                    rec[sc], rec[ec] = str(max(1, _shift(s, off))), str(_shift(e, off))
            rows.append(rec)
            stats.j_targets += 1

        for row in jc_rows:                            # J+C: verbatim, already V-less
            sid = f"JC|{row['scaffold_id']}"
            fa.write(f">{sid}\n{seqs[row['scaffold_id']]}\n")
            rec = {c: str(row.get(c) or "") for c in cols if c != "segment"}
            rec["scaffold_id"], rec["segment"] = sid, "JC"
            rows.append(rec)
            stats.jc_targets += 1

    pl.DataFrame(rows, schema={c: pl.Utf8 for c in cols}).write_csv(out_tsv, separator="\t")
    return stats


# Measured on 20 k TRA amplicon reads (SRR5233635), mmseqs 18-8cc5c, 4 threads, 2 reps:
#
#   reference          targets   search      reads with a hit
#   V x J (shipped)     15,414   12.0 s      10,616
#   segments             1,244    2.5 s      10,814   (+199, one lost)
#
# 4.7x on the whole search step, and recall is slightly BETTER. Of the reads that hit,
# 86.9 % hit BOTH a V target and a J/JC target -- the junction-spanning fraction, matching an
# independent estimate of 85.6 % on the same library.
#
# ⚠ That 86.9 % is exactly why `mmseqs.top_hit` (`filterdb --extract-lines 1`) CANNOT be used
# unchanged against this reference: J+C targets carry ~150 nt of constant region, so they
# outscore a short V segment on raw bits, and taking only the single best hit would discard the
# V call on most amplicon reads. The consumer must take the best V hit AND the best J hit per
# read and merge them. Until that exists in the mapper, this reference is built and validated
# but not wired into `arda rnaseq map`.
