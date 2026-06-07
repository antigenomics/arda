"""Load the curated reference markup for runtime projection.

For nucleotide annotation we use ``markup.tsv`` + ``alleles.fasta``; for amino
acid annotation ``markup.aa.tsv`` + ``alleles.aa.fasta``. Both expose region
``*_start``/``*_end`` columns in the same coordinate space (nt or aa), so the
projection code is identical.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ..paths import vdj_dir

__all__ = ["REGIONS", "RefEntry", "Reference", "load_reference"]

# Canonical region order (matches build output and AIRR field grouping).
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")


@dataclass(slots=True)
class RefEntry:
    """Per-scaffold reference markup: region coords (in target space) + calls."""

    locus: str
    v_call: str
    j_call: str
    starts: list[int]   # one per REGIONS, 1-based closed (target coords)
    ends: list[int]
    v_sequence_end: int = 0    # scaffold nt position of V germline end (0 = unknown)
    j_sequence_start: int = 0  # scaffold nt position of J germline start


@dataclass
class Reference:
    """In-memory reference for one (organism, seqtype)."""

    organism: str
    seqtype: str
    target_fasta: Path
    entries: dict[str, RefEntry]
    d_germlines: dict[str, list[tuple[str, str]]]  # locus -> [(allele, nt_seq)]

    def get(self, scaffold_id: str) -> RefEntry | None:
        return self.entries.get(scaffold_id)


def _load_d_germlines(base: Path) -> dict[str, list[tuple[str, str]]]:
    """Load ``d_germlines.fasta`` (``>locus|allele``) grouped by locus.

    Used for runtime D-segment mapping in nucleotide space. Returns an empty
    mapping if the file is absent (older reference builds, or VJ-only species).
    """
    path = base / "d_germlines.fasta"
    out: dict[str, list[tuple[str, str]]] = {}
    if not path.exists():
        return out
    from ..refbuild.imgt import read_fasta

    for header, seq in read_fasta(path):
        locus, _, allele = header.partition("|")
        if allele and seq:
            out.setdefault(locus, []).append((allele, seq.upper()))
    return out


def load_reference(organism: str, seqtype: str = "nt") -> Reference:
    """Load reference markup + target FASTA path for an organism."""
    base = vdj_dir(organism)
    if not base.is_dir():
        raise FileNotFoundError(
            f"No reference DB for organism {organism!r} at {base}. Run `arda build-db`."
        )
    if seqtype == "aa":
        markup_path = base / "markup.aa.tsv"
        target_fasta = base / "alleles.aa.fasta"
    else:
        markup_path = base / "markup.tsv"
        target_fasta = base / "alleles.fasta"

    df = pl.read_csv(markup_path, separator="\t", infer_schema_length=0)
    start_cols = [f"{r}_start" for r in REGIONS]
    end_cols = [f"{r}_end" for r in REGIONS]
    has_vj = "v_sequence_end" in df.columns and "j_sequence_start" in df.columns

    def _int(v) -> int:
        try:
            return int(v)
        except (TypeError, ValueError):
            return 0

    entries: dict[str, RefEntry] = {}
    for row in df.iter_rows(named=True):
        entries[row["scaffold_id"]] = RefEntry(
            locus=row["locus"],
            v_call=row["v_call"],
            j_call=row["j_call"],
            starts=[int(row[c]) for c in start_cols],
            ends=[int(row[c]) for c in end_cols],
            v_sequence_end=_int(row["v_sequence_end"]) if has_vj else 0,
            j_sequence_start=_int(row["j_sequence_start"]) if has_vj else 0,
        )
    d_germlines = _load_d_germlines(base) if seqtype == "nt" else {}
    return Reference(organism, seqtype, target_fasta, entries, d_germlines)
