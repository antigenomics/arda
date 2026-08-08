"""Export reference sequences with their FR/CDR markup.

The reference is arda's most valuable offline artifact: every in-frame V·J germline scaffold with
IgBLAST-quality FR1-4 / CDR1-3 coordinates, plus the per-segment markup and the per-allele CDR3
anchors. Until now it was only reachable by reading the build's TSVs by hand and re-joining them
against the FASTAs, which is exactly the kind of thing that goes wrong quietly — the coordinates
are 1-based closed, the aa reference has three frames per D allele, and a ``J + C`` scaffold has no
V at all, so a naive join produces plausible nonsense.

This module does the join once, correctly, and emits it in the formats people actually want:

* ``tsv``   — one row per record, sequence plus every region as its own column (nt and aa)
* ``fasta`` — the sequences alone, headed by their calls
* ``gff3``  — the regions as features on the sequence, for a genome browser or a plotting tool
* ``airr``  — the same rows shaped as an AIRR Rearrangement TSV, so a reference scaffold can be fed
  straight into a tool that expects arda's own output

⛔ Coordinates are **1-based closed** everywhere in arda, which is the AIRR convention. GFF3 is also
1-based closed, so those pass through unchanged; anything consuming 0-based half-open (BED, most
Python slicing) must convert, and the ``tsv`` output states the convention in a comment line.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO

from .annotate.reference import REGIONS, Reference, load_reference

__all__ = ["export_reference", "FORMATS", "KINDS"]

#: Output formats. `airr` is deliberately a *Rearrangement* shape, not a germline-set shape: the
#: point is to be able to round-trip a scaffold through anything that reads arda's own output.
FORMATS = ("tsv", "fasta", "gff3", "airr")

#: What to export. `scaffolds` is the V×J (and J+C) reference the mapper aligns against;
#: `segments` is the collapsed per-allele V / J / C reference the two-pass segment search uses;
#: `anchors` is the per-allele CDR3 anchor table (`germline_nt`, `templated_aa`, status).
KINDS = ("scaffolds", "segments", "anchors")


@dataclass(slots=True)
class _Record:
    seq_id: str
    locus: str
    v_call: str
    j_call: str
    c_call: str
    sequence: str
    #: region -> (start, end), 1-based closed; absent when the record has no such region.
    regions: dict[str, tuple[int, int]]
    junction: str
    junction_aa: str
    productive: str


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(parts)
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts)
    return out


def _records(ref: Reference, kind: str) -> Iterator[_Record]:
    base = ref.target_fasta.parent
    if kind == "segments":
        fasta, markup = base / "segments.fasta", base / "segments.markup.tsv"
        if not fasta.exists():
            raise FileNotFoundError(
                f"No segment reference at {fasta}. It is GENERATED, not shipped — run "
                f"`arda build-index --organism {ref.organism} --force`."
            )
    else:
        fasta = ref.target_fasta
        markup = base / ("markup.aa.tsv" if ref.seqtype == "aa" else "markup.tsv")
    seqs = _read_fasta(fasta)

    with open(markup) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["scaffold_id"]
            seq = seqs.get(sid)
            if seq is None:
                # A markup row with no sequence is a broken reference, not a record to emit.
                continue
            regions: dict[str, tuple[int, int]] = {}
            for r in REGIONS:
                try:
                    s, e = int(row[f"{r}_start"]), int(row[f"{r}_end"])
                except (KeyError, TypeError, ValueError):
                    continue
                if s > 0 and e >= s:
                    regions[r] = (s, e)
            yield _Record(
                seq_id=sid, locus=row.get("locus", ""), v_call=row.get("v_call", ""),
                j_call=row.get("j_call", ""), c_call=row.get("c_call", ""), sequence=seq,
                regions=regions, junction=row.get("junction", "") or "",
                junction_aa=row.get("junction_aa", "") or "",
                productive=row.get("productive", "") or "",
            )


def _write_tsv(recs: Iterator[_Record], fh: TextIO, seqtype: str) -> int:
    cols = (["sequence_id", "locus", "v_call", "j_call", "c_call", "productive",
             "junction", "junction_aa", "sequence_length", "sequence"]
            + [f"{r}_{k}" for r in REGIONS for k in ("start", "end", "seq")])
    fh.write(f"# arda reference export ({seqtype}); coordinates are 1-based CLOSED (AIRR/GFF3 "
             f"convention), an empty region means the record does not carry it\n")
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    n = 0
    for rec in recs:
        row = {"sequence_id": rec.seq_id, "locus": rec.locus, "v_call": rec.v_call,
               "j_call": rec.j_call, "c_call": rec.c_call, "productive": rec.productive,
               "junction": rec.junction, "junction_aa": rec.junction_aa,
               "sequence_length": len(rec.sequence), "sequence": rec.sequence}
        for r in REGIONS:
            se = rec.regions.get(r)
            row[f"{r}_start"] = se[0] if se else ""
            row[f"{r}_end"] = se[1] if se else ""
            row[f"{r}_seq"] = rec.sequence[se[0] - 1:se[1]] if se else ""
        w.writerow(row)
        n += 1
    return n


def _write_fasta(recs: Iterator[_Record], fh: TextIO) -> int:
    n = 0
    for rec in recs:
        calls = "|".join(x for x in (rec.v_call, rec.j_call, rec.c_call) if x)
        fh.write(f">{rec.seq_id} locus={rec.locus} {calls}\n")
        for i in range(0, len(rec.sequence), 60):
            fh.write(rec.sequence[i:i + 60] + "\n")
        n += 1
    return n


def _write_gff3(recs: Iterator[_Record], fh: TextIO) -> int:
    # GFF3 is 1-based closed, the same convention arda uses, so the coordinates pass through
    # unchanged. That is worth stating because it is the one export where an off-by-one would be
    # invisible: a browser would happily draw the shifted feature.
    fh.write("##gff-version 3\n")
    n = 0
    for rec in recs:
        fh.write(f"##sequence-region {rec.seq_id} 1 {len(rec.sequence)}\n")
        for r in REGIONS:
            se = rec.regions.get(r)
            if not se:
                continue
            attrs = f"ID={rec.seq_id}:{r};Name={r.upper()}"
            if rec.locus:
                attrs += f";locus={rec.locus}"
            fh.write(f"{rec.seq_id}\tarda\tregion\t{se[0]}\t{se[1]}\t.\t+\t.\t{attrs}\n")
        n += 1
    return n


def _write_airr(recs: Iterator[_Record], fh: TextIO) -> int:
    # The scaffold IS its own perfect alignment, so `sequence_alignment` == `germline_alignment` ==
    # the sequence and every *_start/_end is the region coordinate. That makes an exported scaffold
    # a valid input to anything that consumes arda's Rearrangement output.
    from .annotate.transfer import AIRR_COLUMNS
    w = csv.DictWriter(fh, fieldnames=AIRR_COLUMNS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    n = 0
    for rec in recs:
        row = dict.fromkeys(AIRR_COLUMNS, "")
        row.update(sequence_id=rec.seq_id, sequence=rec.sequence, locus=rec.locus,
                   v_call=rec.v_call, j_call=rec.j_call, c_call=rec.c_call,
                   junction=rec.junction, junction_aa=rec.junction_aa,
                   productive=rec.productive, rev_comp="F",
                   sequence_alignment=rec.sequence, germline_alignment=rec.sequence)
        for r in REGIONS:
            se = rec.regions.get(r)
            if se:
                row[f"{r}_start"], row[f"{r}_end"] = se[0], se[1]
                row[r] = rec.sequence[se[0] - 1:se[1]]
        w.writerow(row)
        n += 1
    return n


def _write_anchors(base: Path, fh: TextIO, loci: set[str] | None) -> int:
    src = base / "cdr3_anchors.tsv"
    if not src.exists():
        raise FileNotFoundError(f"No anchor table at {src}. Run `arda build-db`.")
    with open(src) as src_fh:
        rdr = csv.DictReader(src_fh, delimiter="\t")
        w = csv.DictWriter(fh, fieldnames=rdr.fieldnames, delimiter="\t")
        w.writeheader()
        n = 0
        for row in rdr:
            if loci and row.get("locus") not in loci:
                continue
            w.writerow(row)
            n += 1
    return n


def export_reference(organism: str = "human", *, kind: str = "scaffolds",
                     fmt: str = "tsv", seqtype: str = "nt",
                     loci: set[str] | None = None, out: Path | None = None) -> int:
    """Write the reference (or a locus subset) with markup. Returns the record count."""
    if fmt not in FORMATS:
        raise ValueError(f"format must be one of {FORMATS}, got {fmt!r}")
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
    ref = load_reference(organism, seqtype)
    fh = open(out, "w", newline="") if out else sys.stdout
    try:
        if kind == "anchors":
            # The anchor table is per-allele and has no sequence, so only its own shape makes
            # sense; refusing is better than emitting an empty FASTA.
            if fmt != "tsv":
                raise ValueError("kind='anchors' is a per-allele table and supports fmt='tsv' only")
            return _write_anchors(ref.target_fasta.parent, fh, loci)
        recs = _records(ref, kind)
        if loci:
            recs = (r for r in recs if r.locus in loci)
        if fmt == "tsv":
            return _write_tsv(recs, fh, seqtype)
        if fmt == "fasta":
            return _write_fasta(recs, fh)
        if fmt == "gff3":
            return _write_gff3(recs, fh)
        return _write_airr(recs, fh)
    finally:
        if out:
            fh.close()
