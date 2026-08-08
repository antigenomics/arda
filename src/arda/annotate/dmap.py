"""D-segment mapping on a bare nucleotide junction — no read, no mmseqs search.

``transfer._map_d`` already maps D into the V..J interior of a query, but it is fed
``v_sequence_end`` / ``j_sequence_start`` projected from an mmseqs scaffold hit. A
VDJdb-style record has no read to align: it has a junction and a V/J call. The
per-allele germlines shipped in ``database/vdj/<org>/cdr3_anchors.tsv`` close that
gap, so the interior can be located directly and the existing mapper reused.

**Junction space**, as everywhere in :mod:`arda.cdr3fix`: the input runs Cys104 ->
Phe/Trp118 inclusive.

Finding the interior. The V germline is exact at the junction's 5' end (V/D/J are
not somatically mutated in TCR, and IGH mutation is rare this close to the anchor),
so the V contribution is the longest common prefix of the junction and the V's
CDR3-region germline; the J contribution is the longest common suffix. Validated
against OLGA ground truth on 1300 junctions across human IGH/TRB/TRD and mouse TRB:
the prefix length is exact for 80-85 % of records and never underestimates (it can
overshoot by 1-2 nt when the first N-region base happens to match), and the derived
interior contains the whole true D segment in 94-99 % of records.

Only IGH, TRB and TRD have D germlines; VJ loci return an empty call.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..cdr3fix import load_anchors, resolve_allele, resolve_locus, resolve_species
from .reference import _load_d_germlines
from ..paths import vdj_dir
from .transfer import _map_d

__all__ = ["DCall", "map_d_junction"]


@dataclass
class DCall:
    """D mapping of one junction. Coordinates are 1-based closed, junction space."""

    locus: str = ""
    d_call: str = ""
    d_sequence_start: int = -1
    d_sequence_end: int = -1
    d_support: str = ""
    d2_call: str = ""
    d2_sequence_start: int = -1
    d2_sequence_end: int = -1
    d2_support: str = ""
    np1: str = ""
    np2: str = ""
    np3: str = ""
    v_sequence_end: int = -1     # last junction nt templated by the V germline
    j_sequence_start: int = -1   # first junction nt templated by the J germline
    extra: dict = field(default_factory=dict)   # d_germline_*/d_cigar when unambiguous

    @property
    def called(self) -> bool:
        return bool(self.d_call)

    @property
    def is_dd(self) -> bool:
        return bool(self.d2_call)

    def markup(self, junction_nt: str) -> list[tuple[str, str]]:
        """The junction cut into labelled parts, 5'->3'.

        ``[("V", ...), ("np1", ...), (d_call, ...), ("np2", ...), (d2_call, ...),
        ("np3", ...), ("J", ...)]`` for a tandem D-D, without the last two entries for a
        single D, and ``[("V", ...), ("N", ...), ("J", ...)]`` when no D was called.

        The parts concatenate back to ``junction_nt`` EXACTLY -- that is the contract this
        method exists to make checkable, and it is what a D-D markup consumer needs: the
        AIRR columns alone give it as a set of coordinates and three np strings that it
        must re-derive the D-observed sequence from.

        ⛔ The V-end / np / D-start boundaries INSIDE the junction are not identifiable
        from sequence -- exonuclease chew-back and N/P addition make the partition
        probabilistic. This is one consistent reading of the junction, not ground truth.
        Empty when the V/J split could not be located at all.
        """
        if self.v_sequence_end < 0 or self.j_sequence_start < 1:
            return []
        js = self.j_sequence_start - 1                      # 1-based closed -> slice index
        v, j = junction_nt[: self.v_sequence_end], junction_nt[js:]
        if not self.called:
            return [("V", v), ("N", junction_nt[self.v_sequence_end : js]), ("J", j)]
        d1 = junction_nt[self.d_sequence_start - 1 : self.d_sequence_end]
        parts = [("V", v), ("np1", self.np1), (self.d_call, d1)]
        if self.is_dd:
            d2 = junction_nt[self.d2_sequence_start - 1 : self.d2_sequence_end]
            parts += [("np2", self.np2), (self.d2_call, d2), ("np3", self.np3)]
        else:
            parts.append(("np2", self.np2))
        parts.append(("J", j))
        return parts


def _common_prefix(a: str, b: str) -> int:
    n = 0
    while n < len(a) and n < len(b) and a[n] == b[n]:
        n += 1
    return n


def _common_suffix(a: str, b: str) -> int:
    n = 0
    while n < len(a) and n < len(b) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


def _d_germlines(organism: str) -> dict[str, list[tuple[str, str]]]:
    return _load_d_germlines(vdj_dir(organism))


def map_d_junction(junction_nt: str, v_call: str, j_call: str,
                   species: str = "human", d_max_evalue: float | None = None) -> DCall:
    """Map D (and a tandem second D) into a bare nucleotide junction.

    ``d_max_evalue`` overrides the shipped E-value gate on the D call(s); see
    :func:`arda.annotate.transfer._map_d`. ``None`` keeps the shipped 0.2.
    """
    organism = resolve_species(species)
    junction_nt = (junction_nt or "").strip().upper()
    locus = resolve_locus(v_call, j_call)
    out = DCall(locus=locus)
    if not junction_nt:
        return out

    germlines = _d_germlines(organism).get(locus)
    if not germlines:
        return out                       # VJ locus, or no D germlines for this organism

    anchors = load_anchors(organism)
    v_id = resolve_allele(v_call.split(",")[0], "V", anchors)
    j_id = resolve_allele(j_call.split(",")[0], "J", anchors)
    va, ja = anchors.get(("V", v_id)), anchors.get(("J", j_id))
    if not va or not ja or va.status != "ok" or ja.status != "ok":
        return out

    v_end = _common_prefix(junction_nt, va.germline_nt.upper())
    j_len = _common_suffix(junction_nt, ja.germline_nt.upper())
    j_start = len(junction_nt) - j_len
    if j_start - v_end < 1:
        return out                       # V and J meet or overlap: no interior to search
    out.v_sequence_end, out.j_sequence_start = v_end, j_start + 1

    rec: dict = {}
    _map_d(rec, junction_nt, v_end, j_start + 1, germlines, j_call,
           d_max_evalue=d_max_evalue)
    if not rec.get("d_call"):
        return out

    out.d_call = rec["d_call"]
    out.d_sequence_start = rec["d_sequence_start"]
    out.d_sequence_end = rec["d_sequence_end"]
    out.d_support = rec.get("d_support", "")
    out.d2_call = rec.get("d2_call", "")
    out.d2_sequence_start = rec.get("d2_sequence_start", -1)
    out.d2_sequence_end = rec.get("d2_sequence_end", -1)
    out.d2_support = rec.get("d2_support", "")
    out.np1, out.np2, out.np3 = rec.get("np1", ""), rec.get("np2", ""), rec.get("np3", "")
    out.extra = {k: v for k, v in rec.items()
                 if k.endswith(("_germline_start", "_germline_end", "_cigar"))}
    return out
