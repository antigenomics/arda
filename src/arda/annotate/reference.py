"""Load the curated reference markup for runtime projection.

For nucleotide annotation we use ``markup.tsv`` + ``alleles.fasta``; for amino
acid annotation ``markup.aa.tsv`` + ``alleles.aa.fasta``. Both expose region
``*_start``/``*_end`` columns in the same coordinate space (nt or aa), so the
projection code is identical.
"""

from __future__ import annotations

from dataclasses import dataclass, field
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
    starts: list[int]   # one per REGIONS, 1-based closed (target coords); -1 = region not present
    ends: list[int]
    v_sequence_end: int = 0    # scaffold nt position of V germline end (0 = unknown)
    j_sequence_start: int = 0  # scaffold nt position of J germline start
    c_call: str = ""           # constant genes; set on `J + C` scaffolds only
    # nt length of the V-J part: the scaffold length for a V-J scaffold, the J length for a `J + C`
    # scaffold. A hit with `tstart >= vj_end` lies wholly inside the constant region -- real receptor
    # mRNA carrying no V(D)J, hence no clonotype. 0 = unknown (reference built before this existed).
    vj_end: int = 0

    @property
    def is_jc(self) -> bool:
        """A constant-region scaffold: a J followed by the CH1 exon, with no V."""
        return not self.v_call and bool(self.c_call)


@dataclass
class Reference:
    """In-memory reference for one (organism, seqtype)."""

    organism: str
    seqtype: str
    target_fasta: Path
    entries: dict[str, RefEntry]
    # locus -> [(allele, seq)] in THIS reference's alphabet: one nt entry per allele, or three
    # translated-frame entries per allele when seqtype == "aa" (a trimmed D has no known frame).
    d_germlines: dict[str, list[tuple[str, str]]]
    anchors: dict = field(default_factory=dict)    # (segment, allele) -> cdr3fix.Anchor
    _jc_combos: dict[tuple[str, str], str] | None = field(default=None, repr=False)

    def jc_combinations(self) -> dict[tuple[str, str], str]:
        """``(j_call, c_call) -> J+C scaffold id``, the C-side twin of ``combinations.tsv``.

        The segment reference carries the J and the constant region as SEPARATE targets, because a
        constant sequence shared across a locus' J+C scaffolds is a cross-product and copying it
        through cost **76.4 %** of the segment search's alignments. So a J→C read now names its
        home the same way a V→J read does — by the pair it hit, resolved through this table.

        Derived from the loaded markup rather than from a file: the J+C scaffolds are already
        there, and a second generated artifact is one more thing that can go stale against it.
        """
        if self._jc_combos is None:
            self._jc_combos = {
                (e.j_call, e.c_call): sid
                for sid, e in self.entries.items() if e.is_jc and e.j_call and e.c_call}
        return self._jc_combos

    def segment_j_call(self, name: str) -> str:
        """J allele for a ``JC|`` segment target, which is named by SCAFFOLD id, not by allele.

        Feeding the raw target name into a (V, J) combination lookup silently fails for every
        J->C read; measured, that collapsed the two-pass fast path from 85.3 % to 0.1 %.

        Retained for references built before the constant region became its own ``C|`` target —
        a build and a mapper of different vintages must not silently mis-resolve a J call.
        """
        entry = self.entries.get(name)
        return entry.j_call if entry else name

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


def _load_d_germlines_aa(base: Path) -> dict[str, list[tuple[str, str]]]:
    """The same D set translated in all three reading frames, for aa annotation.

    A D segment is trimmed at both ends before joining, so its reading frame in the junction
    is not knowable from the germline: all three must be searched. Each allele therefore
    contributes three entries under one name, and ``transfer._best_d`` de-duplicates the
    allele list when two frames tie.
    """
    from ..refbuild.translate import translate

    out: dict[str, list[tuple[str, str]]] = {}
    for locus, alleles in _load_d_germlines(base).items():
        for allele, seq in alleles:
            for frame in (0, 1, 2):
                aa = translate(seq[frame:], 0)
                if aa:
                    out.setdefault(locus, []).append((allele, aa))
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

    start_cols = [f"{r}_start" for r in REGIONS]
    end_cols = [f"{r}_end" for r in REGIONS]

    def _int(v) -> int:
        try:
            return int(v)
        except (TypeError, ValueError):
            return 0

    def _load(path: Path, entries: dict[str, RefEntry]) -> None:
        df = pl.read_csv(path, separator="\t", infer_schema_length=0)
        has_vj = "v_sequence_end" in df.columns and "j_sequence_start" in df.columns
        for row in df.iter_rows(named=True):
            entries[row["scaffold_id"]] = RefEntry(
                locus=row["locus"],
                v_call=row["v_call"],
                j_call=row["j_call"],
                starts=[int(row[c]) for c in start_cols],
                ends=[int(row[c]) for c in end_cols],
                v_sequence_end=_int(row["v_sequence_end"]) if has_vj else 0,
                j_sequence_start=_int(row["j_sequence_start"]) if has_vj else 0,
                # absent from reference builds that predate the constant-region scaffolds
                c_call=row.get("c_call") or "",
                vj_end=_int(row.get("vj_end")),
            )

    entries: dict[str, RefEntry] = {}
    _load(markup_path, entries)
    # Segment targets (`V|allele`, `J|allele`, `JC|scaffold`) share the schema and live in the
    # same key space, so they load through the identical path. Without them every segment hit
    # resolves to None and `_annotate_chunk` drops the read as unmapped -- which is why searching
    # the segment reference produced 0 annotated reads against the scaffold reference's 278 on the
    # same input, despite the search finding them.
    #
    # Their coordinates are usable, and that is not an assumption: all 775 V segments carry
    # `fwr1..fwr3` coords IDENTICAL to their scaffolds' (a scaffold is `V + pad + J` with the V at
    # position 1), and FR4's offset *inside* the J is invariant on 15,390 of 15,414 scaffolds with
    # 0 differing. Generated by `build-index`, so absent on a reference built before 2.6.0 --
    # optional by design, never an error.
    seg_markup = markup_path.parent / "segments.markup.tsv"
    if seqtype == "nt" and seg_markup.exists():
        _load(seg_markup, entries)
    d_germlines = _load_d_germlines_aa(base) if seqtype == "aa" else _load_d_germlines(base)
    # Per-allele junction germlines: they pin `v_sequence_end` / `j_sequence_start` far
    # better than projecting the scaffold's N-pad boundaries (see transfer._anchored_vj_bounds).
    # Both alphabets need the anchors: nt reads them as `germline_nt`, aa as `templated_aa`.
    from ..cdr3fix import load_anchors
    anchors = load_anchors(organism)
    return Reference(organism, seqtype, target_fasta, entries, d_germlines, anchors)
