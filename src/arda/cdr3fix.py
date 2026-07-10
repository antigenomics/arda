"""Markup and repair of bare ``(junction_aa, V, J)`` records — the VDJdb case.

**Coordinate convention.** Everything here is *junction* space: Cys104 through the
Phe/Trp118 that opens FR4, **both anchors included**. That is what VDJdb's ``cdr3``
column actually holds (``CASSARSGELFF`` with ``vEnd=4``, ``jStart=7``), and it is
NOT arda's ``cdr3`` (which excludes both anchors). Conflating the two silently
corrupts every coordinate emitted here.

The V and J germlines each template a known run of residues into the junction, and
``database/vdj/<organism>/cdr3_anchors.tsv`` ships them per allele. So marking up a
record needs no germline search: align the junction's 5' end against the V's
templated residues (anchored at Cys104) and its 3' end against the J's (anchored at
[FW]118), and read off the edit operations.

Both alignments are one semi-global Needleman-Wunsch anchored at the conserved
residue with free end gaps on the junction-interior side. The free end gap is what
makes the result honest: the germline templated run is an **upper bound** (V and J
are exonuclease-trimmed), so the alignment stops wherever germline agreement stops
paying for itself, and the untemplated N/D region is never scored. Concretely::

    germline CASS    vs CCSS...   -> sub at index 1, d=1 -> repaired to CASS...
    germline CASS    vs CGGS...   -> v_end = 1, no error (that is the V/N boundary)
    germline TNEKLFF vs ...NEKLF  -> deletion, d=0 -> repaired to ...NEKLFF
    germline TNEKLFF vs ...NNKLFF -> sub at index 8, d=4 -> REPORTED, not repaired

Detection and repair are deliberately separate: every germline disagreement is
reported (with its position, extent and distance from the anchor), but only edits
adjacent to the conserved anchor are applied. See ``_MAX_REPLACE`` for why.

Fix-type names mirror VDJdb's ``Cdr3Fixer`` so its ``cdr3fix`` JSON is directly
comparable; the per-position ``errors`` list is arda's addition.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import polars as pl

from .paths import vdj_dir
from .refbuild.loci import VDJDB_SPECIES

__all__ = [
    "Cdr3Error",
    "Cdr3Markup",
    "Anchor",
    "load_anchors",
    "markup_cdr3",
    "markup_records",
    "markup_batch",
    "to_frame",
    "format_report",
    "MARKUP_COLUMNS",
]

# Alignment scoring. Indels in a germline-templated run are rarer than
# substitutions, so a gap costs more than a mismatch; this is what makes the
# `CGGS` case stop at the Cys instead of forcing a 2-substitution alignment.
_MATCH, _MISMATCH, _GAP = 1, -1, -2

# Beyond this many edits in one germline run we refuse to repair: the record is
# more likely mis-assigned to the wrong allele than to hold that many typos.
_MAX_FIX = 2

# How far from the conserved anchor (Cys104 for V, [FW]118 for J) a mismatch may
# sit and still be *repaired*. This is the crux of the whole module.
#
# The germline templated run is an upper bound -- V and J are exonuclease-trimmed
# -- so a mismatch inside it is ambiguous: a curation typo, or simply the N/D
# region starting earlier than the germline could reach. The alignment cannot tell
# them apart, because a single mismatch only needs two flanking matches to score
# better than stopping, and two chance matches happen ~1/400 per opportunity.
# Repairing on that evidence rewrites real N-region residues: on 3000 VDJdb rows it
# silently "fixed" 84 records, e.g. CASSPRRY-N-L-QFF -> ...NEQFF against TRBJ2-1
# (`SYNEQFF`), where the L is N-region, not a typo.
#
# Adjacent to the conserved anchor the ambiguity collapses: the anchor is fixed, so
# a mismatch beside it cannot be explained away by trimming. VDJdb encodes the same
# prior as `max_replace_size = 1`. Mismatches further in are still REPORTED (with
# `applied=False`) -- the caller asked where the V/J mismatch is -- but never applied.
_MAX_REPLACE = 1

# VDJdb fix types, with its rank order (worst wins when several apply).
_RANK = {
    "NoFixNeeded": 0,
    "FixTrim": 1,
    "FixAdd": 2,
    "FixReplace": 3,
    "FailedBadSegment": 4,
    "FailedReplace": 5,
    "FailedNoAlignment": 6,
}
_GOOD = {"NoFixNeeded", "FixTrim", "FixAdd", "FixReplace"}

MARKUP_COLUMNS = [
    "cdr3", "cdr3_repaired", "v_call", "j_call", "locus", "species",
    "v_end", "j_start", "v_fix", "j_fix", "v_canonical", "j_canonical",
    "good", "fix_needed", "n_errors", "errors", "cdr3fix",
]


@dataclass(frozen=True)
class Anchor:
    """A germline segment's contribution to the junction."""

    locus: str
    segment: str          # "V" | "J"
    templated_aa: str     # V: starts at Cys104. J: ends at [FW]118.
    functionality: str
    status: str           # "ok" | "no_anchor"
    source: str           # "ndm" | "aux" | "motif" | "no_anchor"
    anchor_nt: int = -1   # 0-based offset of the anchor codon in the germline
    partial_nt: int = 0   # V: dangling 3' nt; J: dangling 5' nt (mid-codon)
    germline_nt: str = ""  # V: Cys104 -> 3' end. J: 5' end -> [FW]118 codon end.


@dataclass(frozen=True)
class Cdr3Error:
    """One edit between the observed junction and the germline-templated run.

    ``pos`` indexes the *observed* junction and ``length`` is how far the error
    extends. ``frm`` is what the record has, ``to`` what the germline says.
    ``dist`` is the distance from the conserved anchor, and ``applied`` says
    whether the repair was actually made (see ``_MAX_REPLACE``).
    """

    side: str      # "V" | "J"
    kind: str      # "sub" | "ins" | "del"
    pos: int
    length: int
    frm: str
    to: str
    dist: int = 0
    applied: bool = False

    def __str__(self) -> str:
        mark = "" if self.applied else " (reported, not repaired)"
        if self.kind == "sub":
            body = f"{self.side} sub@{self.pos} {self.frm}>{self.to}"
        elif self.kind == "del":
            body = f"{self.side} del@{self.pos} missing {self.to!r}"
        else:
            body = f"{self.side} ins@{self.pos} extra {self.frm!r}"
        return f"{body} d={self.dist}{mark}"


@dataclass
class Cdr3Markup:
    """Result of marking up one ``(junction_aa, V, J)`` record."""

    cdr3: str                      # as submitted (junction space)
    cdr3_repaired: str
    v_call: str = ""
    j_call: str = ""
    locus: str = ""
    species: str = ""
    v_end: int = -1                # count of V-templated residues
    j_start: int = -1              # index of the first J-templated residue
    v_fix: str = "FailedBadSegment"
    j_fix: str = "FailedBadSegment"
    errors: list[Cdr3Error] = field(default_factory=list)
    sequence_id: str = ""

    @property
    def v_canonical(self) -> bool:
        return self.cdr3.startswith("C")

    @property
    def j_canonical(self) -> bool:
        return self.cdr3.endswith(("F", "W"))

    @property
    def good(self) -> bool:
        return self.v_fix in _GOOD and self.j_fix in _GOOD

    @property
    def fix_needed(self) -> bool:
        return self.cdr3_repaired != self.cdr3

    def to_cdr3fix(self) -> dict:
        """The VDJdb ``cdr3fix`` JSON object, key-for-key."""
        return {
            "cdr3": self.cdr3_repaired,
            "cdr3_old": self.cdr3,
            "fixNeeded": self.fix_needed,
            "good": self.good,
            "jCanonical": self.j_canonical,
            "jFixType": self.j_fix,
            "jId": self.j_call,
            "jStart": self.j_start,
            "vCanonical": self.v_canonical,
            "vEnd": self.v_end,
            "vFixType": self.v_fix,
            "vId": self.v_call,
        }

    def explain(self) -> str:
        """One human-readable line: what happened, why, and where."""
        head = f"{self.sequence_id or self.cdr3}"
        state = "OK" if not self.fix_needed and self.good else (
            "FIXED" if self.fix_needed and self.good else "FAILED")
        parts = [
            f"[{state}] {head}",
            f"V={self.v_call or '?'} vEnd={self.v_end} {self.v_fix}",
            f"J={self.j_call or '?'} jStart={self.j_start} {self.j_fix}",
        ]
        if self.fix_needed:
            parts.append(f"{self.cdr3} -> {self.cdr3_repaired}")
        if self.errors:
            parts.append("; ".join(str(e) for e in self.errors))
        return " | ".join(parts)


# ---------------------------------------------------------------------------
# Anchor table
# ---------------------------------------------------------------------------

def _anchor_path(organism: str) -> Path:
    return vdj_dir(organism) / "cdr3_anchors.tsv"


@lru_cache(maxsize=8)
def load_anchors(organism: str) -> dict[tuple[str, str], Anchor]:
    """``{(segment, allele): Anchor}`` for one organism; ``{}`` if not built."""
    path = _anchor_path(organism)
    if not path.exists():
        return {}
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    out: dict[tuple[str, str], Anchor] = {}
    for r in df.iter_rows(named=True):
        out[(r["segment"], r["allele"])] = Anchor(
            locus=r["locus"], segment=r["segment"], templated_aa=r["templated_aa"] or "",
            functionality=r["functionality"], status=r["status"], source=r["source"],
            anchor_nt=int(r["anchor_nt"]), partial_nt=int(r["partial_nt"]),
            germline_nt=r["germline_nt"] or "")
    return out


def resolve_species(species: str) -> str:
    """VDJdb species name or arda organism -> arda organism."""
    s = (species or "").strip()
    return VDJDB_SPECIES.get(s.lower().replace("_", "").replace(" ", ""), s.lower())


def resolve_locus(v_call: str, j_call: str = "") -> str:
    """``TRBV6-1`` -> ``TRB``; ``TRAV29/DV5`` -> ``TRA`` (leading token wins)."""
    for call in (v_call, j_call):
        token = (call or "").split(",")[0].split("/")[0].strip()
        if len(token) >= 3:
            return token[:3].upper()
    return ""


def resolve_allele(call: str, segment: str, anchors: dict) -> str:
    """VDJdb's ``get_closest_id`` ladder: exact -> gene*01 -> family*01 -> any allele.

    Returns ``""`` when nothing resolves.
    """
    call = (call or "").strip()
    if not call:
        return ""
    if (segment, call) in anchors:
        return call
    gene = call.split("*")[0]
    if (segment, f"{gene}*01") in anchors:
        return f"{gene}*01"
    family = gene.split("-")[0]
    if (segment, f"{family}*01") in anchors:
        return f"{family}*01"
    for (seg, allele) in anchors:
        if seg == segment and allele.split("*")[0] == gene:
            return allele
    return ""


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def _align(germline: str, query: str) -> tuple[int, list[tuple[str, int, int]]]:
    """Semi-global align, anchored at index 0 of both, free end gaps.

    Returns ``(score, ops)`` where each op is ``(kind, gi, qi)`` and kind is
    ``M`` (match), ``X`` (substitution), ``I`` (extra residue in the query) or
    ``D`` (germline residue absent from the query).
    """
    n, m = len(germline), len(query)
    if not n or not m:
        return 0, []
    s = [[0] * (m + 1) for _ in range(n + 1)]
    # Germline residues nearest the anchor may be genuinely ABSENT from a truncated
    # submission (`CASSRGSVRLGTTDPQ` has lost its `YF`; `ASEGGNTIY` has lost its
    # Cys104), so those leading gaps are free -- charging for them made the whole
    # alignment score 0 and the record was silently left unrepaired. Query residues
    # sitting before the germline starts are still charged: those are extra
    # residues to trim, not missing ones to restore.
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + _GAP
    for i in range(1, n + 1):
        gi = germline[i - 1]
        for j in range(1, m + 1):
            diag = s[i - 1][j - 1] + (_MATCH if gi == query[j - 1] else _MISMATCH)
            s[i][j] = max(diag, s[i - 1][j] + _GAP, s[i][j - 1] + _GAP)

    # Free end gaps: stop at the best-scoring cell. Ties prefer (a) consuming more
    # query, so germline agreement is credited as far as it genuinely extends, then
    # (b) consuming LESS germline, so we never invent residues we did not have to.
    #
    # Tie-breaking on `i + j` instead is a trap: for `CYVPGDRGGYTDKLIF` against
    # TRDV2*03 (`CACDT`) it scores "skip the germline CA, then match the C" equal to
    # "match the C", picks the skip, and prepends `CA` to a junction that already
    # begins with the conserved Cys.
    score, bj, neg_i = max((s[i][j], j, -i) for i in range(n + 1) for j in range(m + 1))
    bi = -neg_i
    ops: list[tuple[str, int, int]] = []
    i, j = bi, bj
    while i > 0 and j > 0:
        gi, qj = germline[i - 1], query[j - 1]
        if s[i][j] == s[i - 1][j - 1] + (_MATCH if gi == qj else _MISMATCH):
            ops.append(("M" if gi == qj else "X", i - 1, j - 1))
            i, j = i - 1, j - 1
        elif s[i][j] == s[i - 1][j] + _GAP:
            ops.append(("D", i - 1, j))
            i -= 1
        else:
            ops.append(("I", i, j - 1))
            j -= 1
    while i > 0:
        ops.append(("D", i - 1, j))
        i -= 1
    while j > 0:
        ops.append(("I", i, j - 1))
        j -= 1
    ops.reverse()
    return score, ops


def _repair(germline: str, query: str, ops, max_replace: int) -> str:
    """Rebuild the aligned query run, applying only anchor-adjacent edits.

    ``qj`` counts residues from the conserved anchor (both alignments start there),
    so ``qj <= max_replace`` is exactly "adjacent to the anchor".
    """
    out: list[str] = []
    for kind, gi, qj in ops:
        near = qj <= max_replace
        if kind == "M":
            out.append(query[qj])
        elif kind == "X":
            out.append(germline[gi] if near else query[qj])
        elif kind == "D":
            if near:
                out.append(germline[gi])
        elif kind == "I":
            if not near:
                out.append(query[qj])
    return "".join(out)


def _errors(ops, germline: str, query: str, side: str, rev: bool,
            qlen: int, max_replace: int) -> list[Cdr3Error]:
    """Collapse runs of identical edit ops into ``Cdr3Error`` records."""
    def pos(qj: int) -> int:
        return (qlen - 1 - qj) if rev else qj

    out: list[Cdr3Error] = []
    run: list[tuple[str, int, int]] = []
    for op in list(ops) + [("M", -1, -1)]:
        if run and op[0] != run[0][0]:
            kind = run[0][0]
            gs = "".join(germline[g] for _, g, _ in run)
            qs = "".join(query[q] for _, _, q in run if 0 <= q < len(query))
            if rev:
                gs, qs = gs[::-1], qs[::-1]
            dist = min(q for _, _, q in run)          # residues from the anchor
            p = min(pos(q) for _, _, q in run)
            kw = dict(side=side, pos=p, length=len(run), dist=dist,
                      applied=dist <= max_replace)
            if kind == "X":
                out.append(Cdr3Error(kind="sub", frm=qs, to=gs, **kw))
            elif kind == "D":
                out.append(Cdr3Error(kind="del", frm="", to=gs, **kw))
            elif kind == "I":
                out.append(Cdr3Error(kind="ins", frm=qs, to="", **kw))
            run = []
        if op[0] in ("X", "I", "D"):
            run.append(op)
    out.sort(key=lambda e: e.pos)
    return out


def _fix_type(errs: list[Cdr3Error], aligned: bool) -> str:
    """Fix type reflects what was *applied*; reported-only errors leave it clean."""
    if not aligned:
        return "FailedNoAlignment"
    applied = [e for e in errs if e.applied]
    if not applied:
        return "NoFixNeeded"
    if sum(e.length for e in applied) > _MAX_FIX:
        return "FailedReplace"
    kinds = {e.kind for e in applied}
    if "sub" in kinds:
        return "FixReplace"
    if "del" in kinds:
        return "FixAdd"
    return "FixTrim"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def markup_cdr3(cdr3: str, v_call: str, j_call: str, species: str = "human", *,
                anchors: dict | None = None, sequence_id: str = "",
                max_replace: int = _MAX_REPLACE) -> Cdr3Markup:
    """Mark up and repair one junction. ``cdr3`` is junction space (C..[FW]).

    ``max_replace`` is how far from the conserved anchor an edit may sit and still
    be repaired; edits beyond it are reported with ``applied=False``. Raising it
    repairs more, at the cost of rewriting N-region residues that merely look like
    germline (see ``_MAX_REPLACE``).
    """
    organism = resolve_species(species)
    if anchors is None:
        anchors = load_anchors(organism)
    cdr3 = (cdr3 or "").strip().upper()
    rec = Cdr3Markup(cdr3=cdr3, cdr3_repaired=cdr3, species=organism,
                     locus=resolve_locus(v_call, j_call), sequence_id=sequence_id)
    if not cdr3 or not anchors:
        return rec

    v_id = resolve_allele(v_call.split(",")[0], "V", anchors)
    j_id = resolve_allele(j_call.split(",")[0], "J", anchors)
    rec.v_call, rec.j_call = v_id, j_id
    v_anchor = anchors.get(("V", v_id))
    j_anchor = anchors.get(("J", j_id))

    repaired = cdr3
    # ---- V side: anchored at Cys104 (index 0), free gap toward the N region.
    if v_anchor is None or v_anchor.status != "ok" or not v_anchor.templated_aa:
        rec.v_fix, rec.v_end = "FailedBadSegment", -1
    else:
        g = v_anchor.templated_aa
        q = cdr3[: len(g) + _MAX_FIX]
        score, ops = _align(g, q)
        consumed = sum(1 for k, _, _ in ops if k in "MXI")
        errs = _errors(ops, g, q, "V", False, len(cdr3), max_replace)
        rec.v_end = consumed if score > 0 else -1
        rec.v_fix = _fix_type(errs, score > 0)
        rec.errors.extend(errs)
        if rec.v_fix in _GOOD and any(e.applied for e in errs):
            repaired = _repair(g, q, ops, max_replace) + cdr3[consumed:]

    # ---- J side: anchored at [FW]118 (last index), free gap toward the N region.
    if j_anchor is None or j_anchor.status != "ok" or not j_anchor.templated_aa:
        rec.j_fix, rec.j_start = "FailedBadSegment", -1
    else:
        g = j_anchor.templated_aa[::-1]
        base = repaired
        q = base[::-1][: len(g) + _MAX_FIX]
        score, ops = _align(g, q)
        consumed = sum(1 for k, _, _ in ops if k in "MXI")
        errs = _errors(ops, g, q, "J", True, len(base), max_replace)
        rec.j_start = (len(base) - consumed) if score > 0 else -1
        rec.j_fix = _fix_type(errs, score > 0)
        rec.errors.extend(errs)
        if rec.j_fix in _GOOD and any(e.applied for e in errs):
            repaired = base[: len(base) - consumed] + _repair(g, q, ops, max_replace)[::-1]

    rec.cdr3_repaired = repaired
    rec.errors.sort(key=lambda e: (e.side, e.pos))
    return rec


def markup_records(df: pl.DataFrame, *, cdr3: str = "cdr3", v: str = "v", j: str = "j",
                   species: str = "species", sequence_id: str | None = None,
                   organism: str | None = None,
                   max_replace: int = _MAX_REPLACE) -> list[Cdr3Markup]:
    """Mark up a whole table. Anchors are loaded (and cached) once per organism."""
    out: list[Cdr3Markup] = []
    for r in df.iter_rows(named=True):
        org = organism or resolve_species(str(r.get(species) or "human"))
        out.append(markup_cdr3(
            str(r.get(cdr3) or ""), str(r.get(v) or ""), str(r.get(j) or ""), org,
            anchors=load_anchors(org), max_replace=max_replace,
            sequence_id=str(r.get(sequence_id) or "") if sequence_id else ""))
    return out


def to_frame(records: Iterable[Cdr3Markup]) -> pl.DataFrame:
    """Records -> a TSV-ready frame with the vdjdb-compatible ``cdr3fix`` column."""
    rows = [{
        "cdr3": m.cdr3, "cdr3_repaired": m.cdr3_repaired,
        "v_call": m.v_call, "j_call": m.j_call, "locus": m.locus, "species": m.species,
        "v_end": m.v_end, "j_start": m.j_start, "v_fix": m.v_fix, "j_fix": m.j_fix,
        "v_canonical": m.v_canonical, "j_canonical": m.j_canonical,
        "good": m.good, "fix_needed": m.fix_needed, "n_errors": len(m.errors),
        "errors": "; ".join(str(e) for e in m.errors),
        "cdr3fix": json.dumps(m.to_cdr3fix(), sort_keys=True),
    } for m in records]
    if not rows:
        return pl.DataFrame({c: [] for c in MARKUP_COLUMNS})
    return pl.DataFrame(rows, schema={c: None for c in MARKUP_COLUMNS})


def markup_batch(df: pl.DataFrame, **kw) -> pl.DataFrame:
    """``markup_records`` + ``to_frame``."""
    return to_frame(markup_records(df, **kw))


def format_report(records: Iterable[Cdr3Markup], *, show_ok: bool = False) -> str:
    """Human-readable log: a summary table, then a line per fixed/failed record.

    ``show_ok=True`` also lists the records that needed no repair.
    """
    records = list(records)
    lines: list[str] = []
    n = len(records)
    ok = [r for r in records if r.good and not r.fix_needed]
    fixed = [r for r in records if r.good and r.fix_needed]
    failed = [r for r in records if not r.good]

    lines.append(f"cdr3fix report: {n} records")
    lines.append(f"  correct (no fix needed) : {len(ok)}")
    lines.append(f"  repaired                : {len(fixed)}")
    lines.append(f"  failed                  : {len(failed)}")

    by_fix: dict[tuple[str, str], int] = {}
    for r in records:
        by_fix[(r.v_fix, r.j_fix)] = by_fix.get((r.v_fix, r.j_fix), 0) + 1
    lines.append("")
    lines.append(f"  {'vFixType':<20} {'jFixType':<20} {'count':>7}")
    for (vf, jf), c in sorted(by_fix.items(), key=lambda kv: (-kv[1], kv[0])):
        lines.append(f"  {vf:<20} {jf:<20} {c:>7}")

    errs: dict[str, int] = {}
    for r in records:
        for e in r.errors:
            errs[f"{e.side} {e.kind}"] = errs.get(f"{e.side} {e.kind}", 0) + 1
    if errs:
        lines.append("")
        lines.append("  errors by side/kind:")
        for k, c in sorted(errs.items()):
            lines.append(f"    {k:<10} {c:>7}")

    detail = (ok if show_ok else []) + fixed + failed
    if detail:
        lines.append("")
        lines.append("  --- records ---")
        lines.extend("  " + r.explain() for r in detail)
    return "\n".join(lines) + "\n"
