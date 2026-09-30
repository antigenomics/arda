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

Placement is one **gapless local alignment** per side (``_markup.d_local_align``, arda's C++ D
caller), and the outcome follows from where the hit starts in each sequence -- VDJdb's
``Cdr3Fixer`` table, unchanged::

    start_in_segment  start_in_cdr3
    0                 0              NoFixNeeded
    0                 > 0            FixTrim    -- framework outside the junction
    > 0               0              FixAdd     -- germline the submission cut
    > 0               <= max_replace FixReplace -- substitute at the anchor
    > 0               > max_replace  FailedReplace

Searching every offset is the point: ``CAMYLCASSLFGSPLHF`` against TRBV9 (``CASSV``) has a spurious
``CA`` at offset 0 and the real ``CASS`` at offset 5, and only the longest hit anywhere finds the
second. An engine anchored at offset 0 scores the first and calls the record clean.

**When it is not certain, flag it -- do not change it.** That is the whole policy, and it admits
exactly four edits, in this order:

1. **Re-call** the V or J, when another allele of the locus explains **2 more residues contiguously
   from its own anchor** (:data:`_RECALL_GAIN`, :func:`anchor_depth`) and does not break an anchor
   the submission already holds. Confirming the submitted call is as much the job as replacing it.
2. **Trim** framework sitting past the conserved anchor of the allele that won. There is no N region
   beyond an anchor, so those residues are not junction.
3. **Add** back germline the submission was cut inside of, when the surviving run clears
   ``_MIN_HIT``.
4. **Substitute the first or last residue only**, to the anchor the winning allele encodes, on a
   solid match.

Nothing else is ever rewritten. A residue that disagrees with germline **inside** the templated run
is reported (``mismatch``) and left exactly as submitted, because a curation error and an allele IMGT
does not record are indistinguishable from one junction -- and rewriting them collapsed distinct
records onto one string: ``CAAAETSYDKV{M,R,T,V}F`` on TRAJ50 are four junctions that all became
``CAAAETSYDKVIF``. Four residues observed at one templated position is a suspect call or an
unrecorded allele, not four independent typos.

**A contradicted call changes the allele, not the sequence** (:func:`guess_allele`), which is step 1
above. Measured over VDJdb's corpus, 74.7 % of anchor-adjacent substitutions into an
already-canonical junction were that case, against 1.1 % of untouched records.

**Repair always targets a canonical junction.** ``cdr3_repaired`` is only accepted when it opens
with Cys104 and closes with Phe/Trp118 (``_canonicalise``); otherwise the submission is returned
untouched. Legacy had no such rule. So ``good`` implies canonical, by construction.

**The verdict is a set of flags per side** (:data:`FLAGS`), because one worst-wins label cannot say
both "I trimmed a flank" and "I found a substitution I will not touch". ``v_fix``/``j_fix`` keep
VDJdb's names so its ``cdr3fix`` JSON stays comparable key-for-key.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
import logging
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import polars as pl

from ._markup import d_local_align
from .paths import vdj_dir
from .refbuild.loci import VDJDB_SPECIES
from .refbuild.translate import CODON_TABLE

__all__ = [
    "Cdr3Error",
    "Cdr3Markup",
    "Anchor",
    "load_anchors",
    "markup_cdr3",
    "boundary_nt",
    "markup_records",
    "markup_batch",
    "to_frame",
    "format_report",
    "MARKUP_COLUMNS",
]

# Shortest hit that counts as an alignment at all, and the depth at which germline agreement is
# EVIDENCE rather than the conserved residue being conserved. VDJdb's `min_hit_size`.
#
# One residue of agreement -- on the V side the Cys, on the J side the Phe, i.e. exactly the residue
# that is there by definition -- does not corroborate the call. `CASSQQQQQQQQQF` against TRBJ1-1
# (`NTEAFF`) shares only its terminal Phe; issue #141 raised 2,165 such keys.
#
# ⚠ It is NOT grounds for refusing the record, and that was measured. `CGGSARSGELFF` against TRBV9
# (`CASSV`) agrees on the Cys and nothing else -- structurally the SAME record -- and there it is
# the RIGHT answer, because nothing disagrees: the V is exonuclease-trimmed back to Cys104 and `GGS`
# is N region. No property of a junction separates the two. So a hit below this depth is reported as
# `shallow`, keeps its boundary and keeps `good`; legacy refused it outright. 1,986 keys (1.047 %),
# 1,425 V and 579 J, spread across alleles -- the top one, TRAV12-1*01, is 95 of the 1,425.
# ⛔ It is also the licence to write the CONSERVED ANCHOR -- substitute the junction's first residue
# to Cys104 or its last to the J's [FW]118, or restore a missing one. That is the one edit which
# changes what the record IS, because every consumer trusts `cdr3_repaired` and a junction made
# canonical against the wrong germline is corrupt in a way that looks clean. So arda does not force a
# junction canonical on a thin match: the hit must clear this floor, and the disagreement must be at
# the anchor itself. `guess_allele` runs first, so the allele being scored is the one the junction
# supports rather than the one submitted.
_MIN_HIT = 2

#: Extra agreeing residues another allele must show before it replaces the submitted call, per
#: side. Re-calling is the FIRST step of the repair: confirm the submitted call or replace it, and
#: only then measure the trim, the addition and the anchor against the allele that won.
#:
#: **2 is measured, not chosen.** Judged against ``isalgo/airr_control``'s nucleotide-called V/J
#: (``human.{tra,trb}.ntvj``, junctions seen in >= 3 donors carrying one V and one J throughout),
#: over the 3,059 human TRA and 426 human TRB VDJdb curation keys that appear there. Each cell is
#: "share of keys whose gene agrees with the nucleotide call (calls this margin moved the wrong
#: way)":
#:
#: ======  ===================  ===================  ==================  ===================
#: gain    TRA V                TRA J                TRB V               TRB J
#: ======  ===================  ===================  ==================  ===================
#: none    88.17 %              99.18 %              94.13 %             98.12 %
#: 1       90.62 % (27 broken)  99.54 % (11 broken)  96.48 % (6 broken)  99.77 % (1 broken)
#: **2**   **89.34 % (0)**      **99.84 % (0)**      **96.48 % (0)**     **100.00 % (0)**
#: 3       88.66 % (0)          99.80 % (0)          94.13 % (0)         99.77 % (0)
#: 4       88.17 % (0)          99.80 % (0)          94.13 % (0)         98.83 % (0)
#: ======  ===================  ===================  ==================  ===================
#:
#: 2 is the smallest margin at which NO re-call moves a call away from the nucleotide answer -- at
#: 1 it breaks 20 calls across the four columns -- and every larger margin only leaves errors
#: standing. Separate V and J thresholds are therefore available and, measured, not needed.
#:
#: What the margin counts is :func:`anchor_depth` -- contiguous agreement from the anchor, a
#: mis-read anchor forgiven -- and not the hit's score. That function documents why each of those
#: three words is there; measured on the same four columns, depth corrects more calls than the
#: score does (TRA V 53 against 49, TRB V 11 against 10) and still breaks none.
_RECALL_GAIN = {"V": 2, "J": 2}

#: The per-side verdict vocabulary. A side carries a **set** of these, not one of them, because a
#: single worst-wins label cannot say both "I trimmed a flank" and "I found a substitution I will
#: not touch" -- and collapsing those onto one string is how a declined repair came to read as a
#: clean record. See :attr:`Cdr3Markup.v_flags`.
#:
#: ==============  ========================================================================
#: flag            meaning
#: ==============  ========================================================================
#: ``ok``          the junction agrees with this side's germline; nothing was done
#: ``allele``      markup used an allele the submission did not name
#: ``sub``         residue(s) substituted to germline
#: ``add``         residue(s) restored from germline
#: ``trim``        residue(s) removed -- framework sitting outside the junction
#: ``mismatch``    disagrees with germline INSIDE the templated run; reported, never rewritten
#: ``shallow``     agrees, but below ``_MIN_HIT``: the call is uncorroborated
#: ``impossible``  no usable segment, no alignment, or an ANCHOR the rules cannot restore
#: ==============  ========================================================================
FLAGS = ("ok", "allele", "sub", "add", "trim", "mismatch", "shallow", "impossible")

# Junction residues standing where the anchor should be that a repair may replace with germline.
# This is legacy's `max_replace_size`, and with the default of 1 it is exactly the anchor residue:
# `CASSTQENTEAFL` on TRBJ1-1 (`NTEAFF`) ends in L where the germline says F, with `FAETN` agreeing
# behind it, and one residue is replaced.
#
# It does NOT bound residues INSIDE the germline run -- those are never rewritten at any setting
# (see the module docstring). 2.32.0 read it as a distance from the anchor instead, which admitted a
# substitution one residue in and produced 2,744 rewrites of already-canonical junctions that the
# authoritative release does not ship.
#
# Raising it lets a repair drop more than one junction residue before the germline run: 224 keys of
# VDJdb's corpus change at 2, `CAAPGAGSYQTF` on TRAJ28 among them.
_MAX_REPLACE = 1

# VDJdb fix types, with its rank order (worst wins when several apply). `TruncatedGermline` is
# arda's own, and sorts below every ordinary success: the boundary is real but it is a LOWER
# BOUND, read off a germline record IMGT ships incomplete (see `_MIN_TRUNCATED_AA`).
_RANK = {
    "NoFixNeeded": 0,
    "FixTrim": 1,
    "FixAdd": 2,
    "FixReplace": 3,
    "TruncatedGermline": 4,
    "FailedBadSegment": 5,
    "FailedReplace": 6,
    "FailedNoAlignment": 7,
}
_GOOD = {"NoFixNeeded", "FixTrim", "FixAdd", "FixReplace", "TruncatedGermline"}

# A `status = truncated` anchor is a germline record that stops inside the anchor region, so its
# `templated_aa` is short but correct as far as it goes. Place a boundary from one only when it
# carries at least this many residues.
#
# 3 is where the prefix starts to mean something. Over the 1,101 human V anchors of 2 residues or
# more, `CA` alone is 657 of them (59.7 %), so a 2-residue match is a coincidence rather than
# evidence; at 3 residues there are 62 distinct prefixes and the most common, `CAR`, is 283 of 989
# (28.6 %). This admits 38 of the 63 truncated human V anchors and 43 of the 53 mouse ones, and
# leaves the rest declining as `FailedBadSegment` -- which is the right answer for `C` and `CA`.
_MIN_TRUNCATED_AA = 3

MARKUP_COLUMNS = [
    "cdr3", "cdr3_repaired", "v_call", "j_call", "locus", "species",
    "v_end", "j_start", "v_end_nt", "j_start_nt", "v_fix", "j_fix",
    "v_flags", "j_flags", "v_canonical", "j_canonical",
    "good", "fix_needed", "n_errors", "errors", "cdr3fix",
]


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Anchor:
    """A germline segment's contribution to the junction."""

    locus: str
    segment: str          # "V" | "J"
    templated_aa: str     # V: starts at Cys104. J: ends at [FW]118.
    functionality: str
    status: str           # "ok" | "truncated" | "no_anchor"
    source: str           # "ndm" | "aux" | "motif" | "no_anchor"
    anchor_nt: int = -1   # 0-based offset of the anchor codon in the germline
    partial_nt: int = 0   # V: dangling 3' nt; J: dangling 5' nt (mid-codon)
    germline_nt: str = ""  # V: Cys104 -> 3' end. J: 5' end -> [FW]118 codon end.


@dataclass(frozen=True)
class Cdr3Error:
    """One edit between the observed junction and the germline-templated run.

    ``pos`` indexes the *observed* junction and ``length`` is how far the error
    extends. ``frm`` is what the record has, ``to`` what the germline says.
    ``dist`` is the distance from the conserved anchor.

    ``applied`` is true only when this edit was actually written into
    ``Cdr3Markup.cdr3_repaired``. Being within ``_MAX_REPLACE`` of the anchor makes an
    edit *eligible*; the whole side's repair is still discarded if its fix type comes
    back ``Failed*`` — no alignment, an edit further from the anchor than ``max_replace``, or a
    result that would not be canonical. An error can
    therefore be reported with ``applied=False`` and the junction left alone — which is
    the point: detection and repair are separate decisions.
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
    v_end_nt: int = -1             # junction NUCLEOTIDES the V germline templates
    j_start_nt: int = -1           # index of the first J-templated junction nucleotide
    v_fix: str = "FailedBadSegment"
    j_fix: str = "FailedBadSegment"
    v_flags: tuple[str, ...] = ("impossible",)
    j_flags: tuple[str, ...] = ("impossible",)
    errors: list[Cdr3Error] = field(default_factory=list)
    sequence_id: str = ""

    @property
    def v_canonical(self) -> bool:
        """Does the junction *as repaired* open with the conserved Cys104?

        Read off ``cdr3_repaired``, not the submission -- restoring the anchor is the whole
        point of the repair, and VDJdb's ``vCanonical``/``jCanonical`` mean the same thing.
        Reading the submission instead disagreed with VDJdb on 76 of 250 fixture rows.
        """
        return self.cdr3_repaired.startswith("C")

    @property
    def j_canonical(self) -> bool:
        """Does the junction *as repaired* close with the conserved Phe/Trp118?"""
        return self.cdr3_repaired.endswith(("F", "W"))

    @property
    def good(self) -> bool:
        """Neither side is ``impossible``, and the result carries both conserved anchors.

        A repair exists to produce a canonical junction. One that ends up without its
        Cys104 or its Phe/Trp118 has not repaired the record, it has invented a junction
        nobody submitted -- so it can never be ``good`` (see ``_canonicalise``).

        The anchor is guaranteed by ``_fix_side``, which refuses any repair that would not land
        on the called allele's own anchor residue -- so this needs no separate canonicality
        conjunct, and must not use one: ``v_canonical``/``j_canonical`` are VDJdb's literal
        ``C`` / ``[FW]`` tests, and ``TRBJ2-7*02``'s anchor decodes **V**, which would make a
        perfectly canonical junction for that allele read as broken.

        Read off the FLAGS, not off ``v_fix``/``j_fix``. The fix-type names are VDJdb's
        worst-wins vocabulary, and a declined repair projects onto ``NoFixNeeded`` there because
        that is what VDJdb's scanner returned for it -- its scan stopped at the anchor and never
        looked deeper in. arda does look deeper, so it knows the junction disagrees with the
        germline, and ``good`` has to mean what it says.
        """
        return "impossible" not in self.v_flags and "impossible" not in self.j_flags

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
            f"V={self.v_call or '?'} vEnd={self.v_end} {','.join(self.v_flags)}",
            f"J={self.j_call or '?'} jStart={self.j_start} {','.join(self.j_flags)}",
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


def _decisive(a: "Anchor") -> tuple:
    """The fields that DECIDE the junction, for conflict detection.

    Not the whole record: a TRAV/DV allele legitimately appears twice -- once from the TRA pass and
    once from TRD's `v_shared` -- differing only in `locus`. Comparing everything would flag 15
    benign human duplicates and train the reader to ignore the 3 mouse ones that matter.
    """
    return (a.anchor_nt, a.germline_nt, a.templated_aa, a.status)


@lru_cache(maxsize=8)
def load_anchors(organism: str) -> dict[tuple[str, str], Anchor]:
    """``{(segment, allele): Anchor}`` for one organism; ``{}`` if not built."""
    path = _anchor_path(organism)
    if not path.exists():
        return {}
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    out: dict[tuple[str, str], Anchor] = {}
    for r in df.iter_rows(named=True):
        key = (r["segment"], r["allele"])
        anchor = Anchor(
            locus=r["locus"], segment=r["segment"], templated_aa=r["templated_aa"] or "",
            functionality=r["functionality"], status=r["status"], source=r["source"],
            anchor_nt=int(r["anchor_nt"]), partial_nt=int(r["partial_nt"]),
            germline_nt=r["germline_nt"] or "")
        prev = out.get(key)
        # Compare only the fields that DECIDE the junction. A TRAV/DV allele legitimately appears
        # twice -- once from the TRA pass and once from TRD's `v_shared` -- differing only in
        # `locus`, and warning about those would be 15 lines of noise per human load that trains
        # the reader to ignore the 3 mouse cases that matter.
        if prev is not None and _decisive(prev) != _decisive(anchor):
            # Never: IMGT ships two accessions under one allele name (mouse `IGKV10-96*01` is both
            # AF441451/287 nt and M15520/286 nt; `IGLV2*01` is J00599 and M17529), so the anchor
            # table can carry two rows for one key with DIFFERENT templated_aa and germline_nt.
            # This loader used to be last-wins, which silently decided which junction germline the
            # Cys104 gate would score against -- on 3 mouse alleles, invisibly, depending on file
            # order. Prefer the usable row (`status == "ok"`), then the one that templates MORE
            # junction, and say so, so the choice is deterministic and auditable rather than
            # whatever `iter_rows` yielded last.
            better = max((prev, anchor),
                         key=lambda a: (a.status == "ok", len(a.germline_nt), a.germline_nt))
            if better is not prev:
                out[key] = better
            logger.warning(
                "cdr3_anchors.tsv has conflicting rows for %s %s in %s; keeping status=%s "
                "germline_nt=%s (IMGT ships two accessions under one allele name)",
                key[0], key[1], organism, better.status, better.germline_nt or "-")
            continue
        out[key] = anchor
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


def _leading_call(call: str) -> str:
    """The first candidate of a multi-candidate call: ``A,B`` and ``A+B`` both give ``A``."""
    return (call or "").replace("+", ",").split(",")[0].strip()


def _allele_moved(submitted: str, resolved: str) -> bool:
    """Did markup use an allele the submission did not name? -- the ``allele`` flag.

    This is a proofreading signal, so it has to fire on exactly the cases a curator must look at
    and stay quiet on the ones they need not.

    * a **bare gene** resolving to its ``*01`` is not a finding. VDJdb records genes without
      alleles as a matter of course, and ``*01`` is arda's documented default for them.
    * a **different gene** is always a finding: ``TRAV14 -> TRAV14/DV4``, ``TRBV20 -> TRBV20-1``,
      ``TRBJ2.1 -> TRBJ2-1``.
    * an **allele the submission named explicitly** that markup did not use is always a finding:
      ``TRBV9*99 -> TRBV9*01``. An allele IMGT does not mint should never have reached this module
      -- it is a curation defect, fixed upstream -- so arda records the substitution rather than
      making it silently.
    """
    if not resolved:
        return False
    if _norm_gene(submitted) != _norm_gene(resolved):
        return True
    return "*" in submitted and submitted.strip().upper() != resolved.upper()


def _norm_gene(name: str) -> str:
    """One spelling for one gene name: upper case, no allele suffix, ``.`` written as ``-``.

    The dot is VDJdb-era nomenclature (``TRBJ2.1`` for ``TRBJ2-1``) and appears in curation chunks
    to this day. Everything else here is whitespace and case.
    """
    return name.strip().upper().split("*")[0].replace(".", "-").replace(" ", "")


def _aliases_of(gene: str) -> set[str]:
    """Every spelling of one REFERENCE gene name that a submission might use.

    IMGT files the dual-use TCR genes under a slash -- ``TRAV14/DV4``, ``TRAV15D-1/DV6D-1`` -- and
    a submission may name either half, or write the slash as a dash. So one reference gene answers
    to up to four spellings, and none of them is a prefix of the gene name the earlier rungs
    look for.
    """
    n = _norm_gene(gene)
    out = {n, n.replace("/", "-")}
    out.update(part for part in n.split("/") if part)
    return out


def _alias_index(anchors: dict, segment: str) -> dict[str, frozenset[str]]:
    """``normalised spelling -> {reference gene names carrying it}`` for one segment."""
    idx: dict[str, set[str]] = {}
    for (seg, allele) in anchors:
        if seg != segment:
            continue
        gene = allele.split("*")[0]
        for alias in _aliases_of(gene):
            idx.setdefault(alias, set()).add(gene)
    return {k: frozenset(v) for k, v in idx.items()}


@lru_cache(maxsize=16)
def _cached_alias_index(organism: str, segment: str) -> dict[str, frozenset[str]]:
    """``_alias_index`` over an organism's shipped anchors. Per-process, never on disk.

    Keyed on the organism rather than on the anchors dict, because a cache keyed on ``id(obj)``
    is a correctness bug: CPython reuses an address after collection, so the next table to land
    there would inherit this one's index.
    """
    return _alias_index(load_anchors(organism), segment)


def _query_spellings(call: str) -> list[str]:
    """The spellings to try for a SUBMITTED call, most faithful first.

    Only the trailing ``-1`` is dropped, and only as a last resort. ``TRBV19-1`` and ``TRAJ42-1``
    are submissions naming a suffix IMGT does not use for those genes (the reference carries
    ``TRBV19`` and ``TRAJ42``), but ``TRBV12-1`` is a real mouse gene -- so this spelling is only
    ever reached when the literal name missed, and it still has to hit exactly one gene.
    """
    n = _norm_gene(call)
    out = [n]
    if n.endswith("-1") and len(n) > 2:
        out.append(n[:-2])
    return out


def resolve_allele(call: str, segment: str, anchors: dict, organism: str = "") -> str:
    """Exact -> ``gene*01`` -> first allele of the gene -> the family's one functional gene.

    VDJdb has a ``family*01`` rung in the middle of its ``get_closest_id`` ladder. That literal
    rung is dead: for a dashless gene (``TRBV9``) the family *is* the gene, so the rung above
    already fired; for a dashed one (``IGHV3-23``) it would need an allele literally named
    ``IGHV3*01``, which IMGT does not mint. Measured over all five shipped organisms it is
    reachable for **0** genes.

    What VDJdb's ladder does reach, and this one did not, is a call naming a **family** whose genes
    all carry a suffix: ``TRBV20`` for ``TRBV20-1``, ``TRBV24`` for ``TRBV24-1``. Nothing named
    ``TRBV20`` exists, so every earlier rung misses and the segment fails outright -- 6,291 human
    beta chains of VDJdb, whose CDR3s the germline places without a single edit once the call
    resolves.

    The last rung fills that, and it refuses rather than guesses: the family must contain **exactly
    one functional gene**. Functionality is the whole of the difference between refusing and
    answering for ``TRBV3``, whose family holds ``TRBV3-1`` (F) and ``TRBV3-2`` (P) -- a
    pseudogene cannot be the V of an expressed receptor, so the call is not ambiguous. Where
    several functional genes share the family (``TRBV6`` has five) there is nothing to resolve and
    the caller gets ``""``; VDJdb's ladder would have taken the lowest-numbered one, which is how
    ``TRAV6-7-DV9`` ends up marked up as ``TRAV6-1*01``.

    The last rung checks the call against **every gene name in the vocabulary**, because each rung
    above searches for a name built out of ``gene`` and so cannot reach a reference gene whose name
    the submission is not a prefix of. Measured over VDJdb's corpus, that was the largest single
    failure class: ``TRAV14`` (999 keys) is filed by IMGT as ``TRAV14/DV4``, ``TRAV15D-1-DV6D-1``
    (37) writes that slash as a dash, ``TRBV19-1*01`` (246) and ``TRAJ42-1*01`` (34) name a ``-1``
    suffix IMGT does not use for those genes, and ``TRBJ2.1`` (5) is VDJdb-era dot nomenclature.
    A spelling must land on exactly one gene or it is refused.

    ``organism`` is optional and only selects a cached alias index; the answer does not depend
    on it.

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
    # The rule for a call arda's reference does not carry: strip the allele and take `*01`, and
    # failing that the lowest allele the reference does carry -- `*02`, then `*03`. IMGT zero-pads
    # the suffix to two digits, so `min` over the names IS the lowest-numbered allele, and it is
    # deterministic in a way the previous `for (seg, allele) in anchors: ... return allele` was
    # not: that returned whatever `cdr3_anchors.tsv` happened to be sorted like the day it was
    # built, which put a reference-build detail into the gene call.
    #
    # `TRBV9*99` is the case, and it should never have reached this module -- an allele IMGT does
    # not mint is a curation defect, and VDJdb's proofreading is where it gets fixed. arda's job
    # is to RECORD it: `markup_cdr3` raises the `allele` flag whenever this rung or the `*01` one
    # moves a call the submission named explicitly, so the substitution lands on the curator's
    # worklist instead of happening silently.
    siblings = [allele for (seg, allele) in anchors
                if seg == segment and allele.split("*")[0] == gene]
    if siblings:
        return min(siblings)

    prefix = f"{gene}-"
    family = {
        allele.split("*")[0]
        for (seg, allele), anchor in anchors.items()
        if seg == segment
        and anchor.functionality == "F"
        and allele.split("*")[0].startswith(prefix)
    }
    if len(family) == 1:
        return resolve_allele(family.pop(), segment, anchors, organism)

    # Last rung: check the call against EVERY gene name in the vocabulary, under one normalisation
    # (`_norm_gene`) and one safety rule -- the spelling must land on exactly one gene. Everything
    # above searches for a name built out of `gene`, so none of it can reach a reference gene whose
    # name is not a prefix of the submission: `TRAV14` is filed as `TRAV14/DV4`, and
    # `TRAV15D-1-DV6D-1` writes that slash as a dash.
    #
    # The uniqueness rule is the whole of the difference between resolving and guessing, and it is
    # the same stance the family rung takes: `TRBV6` names five functional genes and `TRBJ2`
    # thirteen, so those stay refused. Ambiguity is reported, never broken by picking the
    # lowest-numbered candidate -- which is how `TRAV6-7-DV9` used to become `TRAV6-1*01`.
    idx = _cached_alias_index(organism, segment) if organism else _alias_index(anchors, segment)
    for spelling in _query_spellings(call):
        hit = idx.get(spelling)
        if hit and len(hit) == 1:
            resolved = next(iter(hit))
            if resolved != gene:                      # never recurse on the name we came in with
                return resolve_allele(resolved, segment, anchors, organism)
    return ""


# ---------------------------------------------------------------------------
# The fixer: legacy's positional logic, on arda's aligner
# ---------------------------------------------------------------------------
#
# This is a port of VDJdb's `Cdr3Fixer`, which arda 2.16.0 replaced with a semi-global
# Needleman-Wunsch anchored at the conserved residue. The NW engine is gone; what follows is
# legacy's decision table, driven by arda's C++ local aligner and arda's own reference.
#
# **What legacy scanned, and why `templated_aa` is it.** `Cdr3Fixer._load_segments_data` slices
# the germline with flanks -- `sequence[reference_point - 3:]` for V, `sequence[:reference_point + 4]`
# for J -- and translates. VDJdb's `reference_point` is offset differently from arda's `anchor_nt`,
# and the sequence those slices actually yield is the templated run itself: Cys104 onward for V,
# up to and including [FW]118 for J. The authoritative 2026-06-03 release proves it, because it
# trims the framework off BOTH ends completely --
#
#     YFCASSQSPGGVAFFGQG -> CASSQSPGGVAFF   vFixType=FixTrim vEnd=5  jFixType=FixTrim jStart=10
#
# -- which a segment carrying one flanking residue cannot do: it would keep that residue and hand
# back `FCASSQSPGGVAFFG`. arda ships the templated run per allele in `cdr3_anchors.tsv`
# (`Anchor.templated_aa`), with the anchor READ from the table rather than found by motif, which
# matters because TRBJ2-7*02's anchor decodes V and TRAJ35*01's decodes Cys.
#
# **What legacy had that the NW engine did not: placement anywhere in the junction.** `KmerScanner`
# takes the LONGEST common substring at any offset, so `CAMYLCASSLFGSPLHF` against TRBV9 (`CASSV`)
# finds `CASS` at query offset 5 and trims five residues (release: vEnd=4). Anchored at index 0 the
# NW engine instead scores the spurious `CA` at offset 0 and reports the record clean -- and no trim
# budget fixes that, because the wrong placement wins on score.
#
# `_markup.d_local_align` is that placement, as a real alignment rather than an exact scan: gapless
# local alignment, match +1 / mismatch -1, over every diagonal, returning the best segment's 0-based
# inclusive offsets in both sequences. Legacy's exact k-mer hit is its zero-mismatch special case.
# It is character-generic, so it runs on residues, and it is the same C++ the D caller uses.

@dataclass(frozen=True)
class Hit:
    """Where a germline segment best sits in a junction end.

    ``start_in_segment`` and ``start_in_cdr3`` are legacy's two offsets, and between them they
    decide the whole outcome (:func:`_fix_side`). ``size`` is the aligned length and ``score`` the
    match-minus-mismatch total, which is what ``_MIN_HIT`` is read against.
    """

    start_in_segment: int
    start_in_cdr3: int
    size: int
    score: int
    lead: int = 0


def _extend(segment: str, cdr3: str, off_s: int, off_q: int) -> tuple[int, int, int]:
    """Walk outward from a placement. Returns ``(size, score, lead)`` of the run it explains.

    A mismatch is **internal to a germline run** -- a typo -- only when the agreement after it is
    at least as long as the run of mismatches itself. Otherwise the germline has stopped and what
    follows is N region, which is the whole reason a templated run is an upper bound.

    That one rule keeps both documented cases. ``TNEKLFF`` against ``...NNKLFF`` has a single
    mismatched residue followed by a match, so it is reported (and, four residues from the anchor,
    not repaired). ``CASSV`` against ``CGGS...`` has TWO mismatches followed by one match, so the
    V simply ends at Cys104 and ``GGS`` is N region -- no error, ``v_end = 1``.
    """
    size = score = 0
    i = 0
    run: list[int] = []
    while off_s + i < len(segment) and off_q + i < len(cdr3):
        if segment[off_s + i] == cdr3[off_q + i]:
            run.append(i)
        else:
            # how far does agreement go after this mismatch run?
            j, matched = i, 0
            while (off_s + j < len(segment) and off_q + j < len(cdr3)
                   and segment[off_s + j] != cdr3[off_q + j]):
                j += 1
            bad = j - i
            k = j
            while (off_s + k < len(segment) and off_q + k < len(cdr3)
                   and segment[off_s + k] == cdr3[off_q + k]):
                k += 1
                matched += 1
            if matched < bad:
                break                      # the germline has stopped; the rest is N region
            i = j
            continue
        i += 1
    # `size` spans everything the run covers, mismatches included; `score` counts agreement only;
    # `lead` is the LEADING exact run, which is what legacy's `match_size` -- and therefore `vEnd`
    # and `len - jStart` -- actually reported. Counting a mismatched residue as templated moved the
    # boundary off the release on 2-4 % of records.
    size = i
    score = sum(1 for k in range(size) if segment[off_s + k] == cdr3[off_q + k])
    lead = 0
    while lead < size and segment[off_s + lead] == cdr3[off_q + lead]:
        lead += 1
    return size, score, lead


def scan(segment: str, cdr3: str) -> Hit | None:
    """Best placement of ``segment`` in ``cdr3``, or ``None`` if nothing aligns.

    Both strings are given anchor-first: the V side as submitted, the J side reversed. So offset 0
    is the conserved residue in each, and "the hit starts at 0 in the segment" means the germline's
    anchor is where the junction says it is.

    ``_markup.d_local_align`` supplies the placement -- gapless local alignment over every diagonal,
    which is legacy's longest-common-substring scan with mismatches allowed. ``CAMYLCASSLFGSPLHF``
    against ``CASSV`` needs exactly that: the real ``CASS`` is at offset 5 and only a search over
    offsets finds it.

    ⚠ **The ANCHORED placement wins a tie**, and that is not cosmetic. ``CASSQQQQQQQQQF`` against
    TRBJ1-1 (``NTEAFF``, reversed ``FFAETN``) scores 1 either way -- the junction's terminal Phe can
    match the germline's Phe118 or its FR4-side second F -- and taking the shifted one made the
    record a ``FixAdd`` that appended an F to a junction already ending in one. The conserved anchor
    is a fixed point, so when matching it explains the junction as well as anything else, it is the
    answer.
    """
    score, i_start, i_end, d_start, d_end = d_local_align(cdr3, segment)
    if i_start < 0:
        return None
    best: Hit | None = None
    if segment and cdr3:
        size0, score0, lead0 = _extend(segment, cdr3, 0, 0)
        if score0 > 0:
            best = Hit(0, 0, size0, score0, lead0)
    if d_start >= 0:
        size, sc, lead = _extend(segment, cdr3, d_start, i_start)
        hit = Hit(d_start, i_start, size, sc, lead)
        if best is None or sc > best.score:      # ties keep the anchored placement
            best = hit
    return best if best is not None and best.score >= 1 else None


def _canonicalise(candidate: str, side: str, anchor_aa: str = "") -> bool:
    """May this repaired junction be accepted for ``side``?

    A repair exists to restore the conserved anchor. One that hands back a junction without it has
    not repaired anything -- it has produced a string nobody submitted, and downstream that string
    will be trusted. Legacy had no such rule; it reported `vCanonical`/`jCanonical` after the fact
    and let a non-canonical repair ship.

    ⚠ ``anchor_aa`` is the residue THIS ALLELE actually encodes, and passing it matters: a
    conserved-motif check is not an anchor. ``TRBJ2-7*02`` templates ``SYEQYV`` -- its anchor codon
    decodes **V**, not [FW] -- so ``CASSKRGGYEQYV`` is a canonical junction for it, and the [FW]
    test calls it broken. VDJdb rewrote exactly that record's terminal V to F to force the
    functional *01. ``TRAJ35*01`` is the same trap on the other residue: its anchor decodes Cys.
    The motif is only the fallback for when no anchor is known.
    """
    if not candidate:
        return False
    if anchor_aa:
        return candidate.startswith(anchor_aa) if side == "V" else candidate.endswith(anchor_aa)
    return candidate.startswith("C") if side == "V" else candidate.endswith(("F", "W"))


def _fix_side(segment: str, cdr3: str, side: str, max_replace: int,
              ) -> tuple[str, tuple[str, ...], str, int, list[Cdr3Error]]:
    """One side, anchor-first. Returns ``(fixed, flags, legacy_name, templated, errors)``.

    Legacy's table, unchanged -- the two offsets are the whole answer, and no scoring, per-end
    floor or tie rule enters after the scan:

    ==================  =================  ==============================================
    ``start_in_segment``  ``start_in_cdr3``  outcome
    ==================  =================  ==============================================
    0                   0                  ``NoFixNeeded``
    0                   > 0                ``FixTrim``    -- drop the flank before the anchor
    > 0                 0                  ``FixAdd``     -- restore the germline that was cut
    > 0                 <= ``max_replace``  ``FixReplace`` -- substitute at the anchor
    > 0                 > ``max_replace``   ``FailedReplace``
    ==================  =================  ==============================================

    A mismatch INSIDE the hit is handled after the table, per residue, and only repaired within
    ``max_replace`` of the anchor.

    ``templated`` is the hit's LEADING EXACT run, which is what legacy reported: ``fix_both`` passes
    ``hit.match_size`` straight through as ``vEnd`` and as ``len(cdr3) - jStart``. So it counts
    residues the germline is OBSERVED to explain, and never the ones a ``FixAdd`` restored --
    ``CAISGEFGSGA`` -> ``CAISGEFGSGANVLTF`` reports ``jStart = 13``, the surviving ``SGA``, not
    ``8``. Counting restored residues as templated is defensible and is not what VDJdb's column
    means, and this value feeds ``boundary_nt`` and `dpost`'s slice of the non-templated middle.
    """
    hit = scan(segment, cdr3)
    # A junction that lost ONLY its conserved anchor is the commonest curation loss, and neither
    # legacy nor 2.32.0 repaired it: the surviving residues agree with the germline from position 1
    # on, so the best placement is a one-residue coincidence below `_MIN_HIT` and both engines
    # returned `FailedNoAlignment`. `AQGLLTGGGNKLTF` on TRAV29/DV5 (`CAAS`) is the shape of it --
    # the 2026-06-03 release ships it unrepaired and not good.
    #
    # Prepending the anchor is admissible exactly when doing so makes the germline agree from the
    # anchor outward: `CAQGLLTGGGNKLTF` places at 0/0 with `CA` agreeing, which clears `_MIN_HIT`.
    # `CALRPA` on TRAJ17 does not -- prepending its Phe still leaves a one-residue hit -- so that
    # record stays refused, as the release has it, and no germline is invented off a coincidence.
    added = ""
    if segment and cdr3 and cdr3[0] != segment[0] and (hit is None or hit.score < _MIN_HIT):
        probe = scan(segment, segment[0] + cdr3)
        if (probe is not None and probe.start_in_segment == 0 and probe.start_in_cdr3 == 0
                and probe.score >= _MIN_HIT):
            added, cdr3, hit = segment[0], segment[0] + cdr3, probe
    if hit is None:
        return cdr3, ("impossible",), "FailedNoAlignment", -1, []

    # Errors are collected anchor-first and converted at return, because `pos` cannot be resolved
    # until the junction's final length is known: a trim shortens it, so a mismatch's distance from
    # the anchor and its index in the output are not the same offset.
    #
    # ⚠ Two coordinate systems, and each is the only meaningful one for its kind. A `trim` reports
    # where the removed residues sat in the SUBMISSION -- they are not in the output at all -- and
    # everything else reports where it landed in the REPAIRED junction, which is what the caller
    # holds. `CASSSPLLSSDTQYFG` on TRBJ2-3 shows both at once: `J ins@15 extra 'G'` (submitted) and
    # `J sub@9 S>T d=5` (repaired, 15 residues long).
    raw: list[tuple[str, int, int, str, str, bool, bool]] = []

    def err(kind: str, pos: int, length: int, frm: str, to: str, applied: bool) -> None:
        raw.append((kind, pos, length, frm, to, applied, kind == "ins"))

    def build(final: str) -> list[Cdr3Error]:
        out = []
        for kind, off, length, frm, to, applied, submitted in raw:
            n = len(cdr3) if submitted else len(final)
            out.append(Cdr3Error(
                side=side, kind=kind,
                pos=(n - 1 - off) if side == "J" else off,
                length=length,
                frm=frm[::-1] if side == "J" else frm,
                to=to[::-1] if side == "J" else to,
                dist=off, applied=applied))
        return out

    def refuse() -> tuple[str, tuple[str, ...], str, int, list[Cdr3Error]]:
        return cdr3, ("impossible",), "FailedReplace", -1, _reported_only(build(cdr3))

    # ---- offsets first: legacy's table decides trim / add / replace / refuse.
    body = cdr3           # the junction with the offset edit applied, anchor-first
    flags: tuple[str, ...] = ()
    name = "NoFixNeeded"
    if added:
        # `cdr3` above is already the prepended junction, so the table sees a hit at 0/0 and has
        # nothing left to do. The restored residue is germline, so -- as in every `FixAdd` -- it is
        # not counted as templated: the boundary is measured from the residues that survived.
        err("del", 0, 1, "", added, True)
        flags, name = ("add",), "FixAdd"
    # Where the aligned run sits AFTER the offset edit -- in `body` and in `segment`. Each branch
    # moves the junction differently, so this cannot be assumed: reading the run at body[0:] while
    # comparing it to segment[start_in_segment:] invented mismatches and refused 100 of legacy's
    # own fixes.
    b_off, s_off = hit.start_in_cdr3, hit.start_in_segment
    if added:
        b_off = s_off = len(added)
    if hit.start_in_segment == 0 and hit.start_in_cdr3 > 0:
        # Framework sitting outside the junction. There is no N region beyond a conserved anchor,
        # so residues past it are always framework -- which is why legacy needed no trim budget to
        # take five residues off `CAMYLCASSLFGSPLHF`.
        err("ins", 0, hit.start_in_cdr3, cdr3[:hit.start_in_cdr3], "", True)
        body, flags, name = cdr3[hit.start_in_cdr3:], ("trim",), "FixTrim"
        b_off, s_off = 0, 0
    elif hit.start_in_segment > 0 and hit.start_in_cdr3 == 0:
        if hit.score < _MIN_HIT:
            # Restoring germline off a thin match is a guess, and legacy refused it too:
            # `min_hit_size = 2` is why the release leaves `CALRPA` (TRAJ17) alone at
            # `FailedNoAlignment` instead of handing back `CALRPAAGNKLTF`. The floor gates repairs
            # that INVENT residues, never the record itself -- a shallow hit that needs no edit
            # keeps its boundary and is reported `shallow`.
            return cdr3, ("impossible",), "FailedNoAlignment", -1, []
        # The junction was cut inside the germline: the anchor and everything up to the hit are
        # gone. Restoring them is the only route to a canonical junction, so legacy put no bound
        # on it and neither does this.
        err("del", 0, hit.start_in_segment, "", segment[:hit.start_in_segment], True)
        body, flags, name = segment[:hit.start_in_segment] + cdr3, ("add",), "FixAdd"
        b_off, s_off = hit.start_in_segment, hit.start_in_segment
    elif hit.start_in_segment > 0 and hit.start_in_cdr3 > 0:
        lead_g, lead_q = segment[:hit.start_in_segment], cdr3[:hit.start_in_cdr3]
        # ⛔ Legacy's gate, and the author's rule of 2026-09-30 is the same one: what may be
        # substituted is the ANCHOR, and `max_replace` bounds the junction side alone. With
        # `_MAX_REPLACE = 1` that is exactly one junction residue -- the one standing where the
        # anchor should be -- so no edit here can reach a residue inside the junction.
        #
        # It does NOT bound how much germline is restored, which is why `CAADNNARLF` on TRAJ31
        # (`NNNARLMF`) works: the hit places at segment 2 / query 1, so one junction residue is
        # dropped and two germline residues take its place, and the release ships `CAADNNARLMF`.
        # Bounding the germline side instead refused 100 of legacy's own fixes.
        #
        # Counting DISAGREEMENTS instead of the shift (2.32.0) is what admitted 2,758 substitutions
        # the authoritative release never made: `CCSSARSGELFF` on TRBV9 (`CASSV`) places on `SS` at
        # offset 2/2, and rewriting its position-1 C to the germline's A is a residue arda invents.
        # A residue that disagrees one in from the anchor may be a curation error or an allele IMGT
        # does not record, and nothing in the junction separates the two.
        if hit.score < _MIN_HIT:
            return cdr3, ("impossible",), "FailedNoAlignment", -1, []
        if hit.start_in_cdr3 > max_replace:
            err("sub", 0, hit.start_in_cdr3, lead_q, lead_g, False)
            return cdr3, ("impossible",), "FailedReplace", -1, build(cdr3)
        subs = [k for k, (x, y) in enumerate(zip(lead_g, lead_q)) if x != y]
        for k in subs:
            err("sub", k, 1, lead_q[k], lead_g[k], True)
        body = lead_g + cdr3[hit.start_in_cdr3:]
        flags = tuple(f for f, on in (("sub", bool(subs)),
                                      ("add", len(lead_g) > len(lead_q)),
                                      ("trim", len(lead_g) < len(lead_q))) if on)
        name = "FixReplace"
        b_off, s_off = hit.start_in_segment, hit.start_in_segment

    # ---- then the residues the run covers, one at a time. `_extend` walks past the first
    # mismatch when the germline agreement after it says the mismatch is internal, so `hit.size`
    # can span typos: `TNEKLFF` against `...NNKLFF` covers 7 residues with one disagreeing. Only
    # those within `max_replace` of the anchor are applied (see `_MAX_REPLACE`); the rest are
    # reported, which is what makes the side `impossible` rather than clean.
    declined = False       # the ANCHOR disagrees and the rules cannot restore it
    mismatch = False       # a disagreement inside the run: reported, left exactly as submitted
    out = list(body)
    for k in range(hit.size):
        qi, si = k + b_off, k + s_off
        if qi >= len(out) or si >= len(segment) or out[qi] == segment[si]:
            continue
        if qi != 0:
            # ⛔ Author's rule, 2026-09-30: only the first/last residue is ever substituted, and
            # only to the anchor. A residue that disagrees with germline INSIDE the templated run
            # may be a curation error or an allele IMGT does not record -- we cannot tell, so it is
            # reported and the submission stands. Rewriting these collapsed distinct records onto
            # one string: `CAAAETSYDKV{M,R,T,V}F` on TRAJ50 are four junctions that all became
            # `CAAAETSYDKVIF`, which is four residues observed at one templated position and reads
            # as a suspect J call or an unrecorded allele, not as four independent typos.
            err("sub", qi, 1, out[qi], segment[si], False)
            mismatch = True
        elif hit.score >= _MIN_HIT:
            # The anchor itself, on a solid match: restore it. This is legacy's FixReplace and the
            # third step of the rule -- `CASSTQENTEAFL` on TRBJ1-1 ends in L where the germline
            # says F, with `FAETN` agreeing behind it.
            err("sub", qi, 1, out[qi], segment[si], True)
            out[qi] = segment[si]
            if "sub" not in flags:
                flags = flags + ("sub",)
            name = "FixReplace"
        else:
            # The anchor, on a thin match: no allele evidence, so no edit and no verdict.
            err("sub", qi, 1, out[qi], segment[si], False)
            declined = True
    body = "".join(out)

    if not flags:
        flags = ("ok",) if hit.score >= _MIN_HIT else ("shallow",)
    if mismatch:
        # NOT a refusal. The junction is well formed and every residue is the curator's; what this
        # says is that one of them is not the germline's. `v_fix`/`j_fix` therefore keep legacy's
        # name -- its substring scanner never saw inside the run either -- and the flag is the
        # extra information. A consumer that wants legacy's strictness reads the flag.
        flags = tuple(f for f in flags if f != "ok") + ("mismatch",)
    if declined:
        if name == "NoFixNeeded":
            name = "FailedReplace"
        # Part of this side may have been repaired and part not. The repaired part is real and
        # named; the unrepaired disagreement still means the side does not agree with germline,
        # and that is what `impossible` says. Reading only the APPLIED edits is how a declined
        # repair came to report itself `NoFixNeeded` and `good` on 19,436 keys.
        flags = tuple(f for f in flags if f != "ok") + ("impossible",)
    # `segment[0]` IS the anchor residue: both runs are given anchor-first.
    anchor_aa = segment[0] if segment else ""
    if body != cdr3 and not _canonicalise(body[::-1] if side == "J" else body, side, anchor_aa):
        return refuse()

    # The boundary is the leading exact run measured on the REPAIRED body, not on the submission.
    # A substitution at the anchor itself makes the pre-repair run 0 -- `FLVGPQGSSASKIIF` against
    # TRAV4 (`CLVGD`) mismatches at offset 0 -- while the residue it writes is germline by
    # construction, so the junction that comes out does start with a templated run.
    lead = 0
    while (b_off + lead < len(body) and s_off + lead < len(segment)
           and body[b_off + lead] == segment[s_off + lead]):
        lead += 1
    return body, flags, name, lead, build(body)


def _reported_only(errs: list[Cdr3Error]) -> list[Cdr3Error]:
    """Clear ``applied`` on a side whose repair was discarded."""
    return [replace(e, applied=False) if e.applied else e for e in errs]


def _with_allele(flags: tuple[str, ...], moved: bool) -> tuple[str, ...]:
    """Add ``allele`` to a side's flags, keeping it exclusive with ``ok``.

    ``ok`` means nothing needed doing; a substituted call is something that was done.
    """
    if not moved:
        return flags
    return ("allele",) + tuple(f for f in flags if f != "ok")


def anchor_depth(segment: str, cdr3: str, start_in_cdr3: int = 0) -> int:
    """Residues the germline explains CONTIGUOUSLY from its anchor, a mis-read anchor forgiven.

    This is what :data:`_RECALL_GAIN` compares, and each of its three properties is there because
    the alternative was measured to re-call a record wrongly:

    * **Contiguous**, not match-minus-mismatch. ``_extend`` walks past a mismatch the agreement
      behind it pays for, so a score can be collected from an alternating pattern --
      ``CASSQQQQQQQQQF`` scores 3 against ``TRBJ1-5`` (``NQPQHF``) by matching every other residue,
      which is Q-richness, not germline evidence.
    * **From the anchor**, so a run found deeper inside the junction does not count: ``CSAR`` is
      ``TRBV20-1``'s run and it does sit inside ``CGGSARSGELFF``, at offset 3, where it is the N
      region reproducing three residues by chance.
    * **A mis-read anchor is forgiven** -- but only when something behind it agrees. Counting from
      offset 0 alone returns 0 for ``CAAAGDNNRKLIF`` against ``TRAJ38`` (``NAGNNRKLIW``), whose
      terminal residue is F where the germline says W, and that zero re-called a junction with six
      agreeing residues to an allele with two. The anchor is the one residue a solid match is
      allowed to rewrite (step three of the rule), so it cannot also be what decides the match.
    """
    n = min(len(segment), len(cdr3) - start_in_cdr3)
    run = 0
    while run < n and segment[run] == cdr3[start_in_cdr3 + run]:
        run += 1
    if run:
        return run
    run = 1
    while run < n and segment[run] == cdr3[start_in_cdr3 + run]:
        run += 1
    return run if run > 1 else 0


def _candidates(anchors: dict, segment: str, locus: str) -> tuple[tuple[str, str], ...]:
    """``(templated run, allele)`` per DISTINCT run for one side of one locus, best allele first.

    Two reductions, and re-calling every record is only affordable with both. The anchor table is
    the whole organism -- 4,000-odd alleles across every locus -- so filtering it per record in
    Python cost 19x the markup itself. And alleles of one locus share their anchor region heavily:
    58 human TRAJ alleles carry 51 distinct runs, 1,000 human IGHV alleles carry 205, so scanning
    per distinct run and mapping back is most of what is left.

    Ordered ``(functional first, then by name)``, so which allele represents a run it shares is
    fixed by the reference and not by how ``cdr3_anchors.tsv`` happened to be sorted.

    ⚠ ``status == "ok"`` only. An allele IMGT ships as a partial record has a SHORT templated run
    for a reason that is about the reference, not about the rearrangement, so it can neither win a
    re-call nor lose one: ``TRBV5-1*02`` carries ``CAS`` where ``*01`` carries ``CASSL``, and
    comparing their lengths would re-call every ``*02`` record to ``*01``.

    ⚠ **Pseudogenes are not candidates**; ORFs are. A junction in hand came from a rearrangement
    that produced a chain, so a pseudogene did not template it -- and a pseudogene's anchor is
    whatever its sequence happens to decode, which `_canonicalise` then accepts: ``TRBJ2-2P``
    (``GRLGG``) won a re-call on ``CASRPGAAGGRPELYF`` and trimmed it to ``CASRPGAAG``, canonical
    for that pseudogene and a junction for nothing. ORFs stay because ``TRBJ2-7*02`` is one and is
    the right answer for ``CASSKRGGYEQYV``, whose terminal V its germline templates.
    """
    by_run: dict[str, list[tuple[bool, str]]] = {}
    for (seg, allele), anchor in anchors.items():
        if seg != segment or anchor.locus != locus or anchor.status != "ok":
            continue
        if anchor.functionality == "P":
            continue
        run = anchor.templated_aa[::-1] if segment == "J" else anchor.templated_aa
        by_run.setdefault(run, []).append((anchor.functionality != "F", allele))
    return tuple((run, min(owners)[1]) for run, owners in by_run.items())


@lru_cache(maxsize=64)
def _cached_candidates(organism: str, segment: str, locus: str) -> tuple[tuple[str, str], ...]:
    """``_candidates`` over an organism's shipped anchors; keyed on the organism, never ``id()``."""
    return _candidates(load_anchors(organism), segment, locus)


@lru_cache(maxsize=64)
def _cached_candidate_index(organism: str, segment: str,
                            locus: str) -> dict[str, tuple[tuple[str, str], ...]]:
    """Candidates grouped by their SECOND templated residue, which is a necessary condition.

    A candidate has to reach ``anchor_depth >= 2`` to clear any margin >= 2, and both ways of
    reaching 2 require the residue at offset 1 to match: contiguously from the anchor it is the
    second of the two, and with the anchor forgiven it is the first one counted. So a junction only
    ever has to be scanned against the alleles whose run carries its own offset-1 residue -- three
    to five of them, against 60 for a whole locus, and the re-call is a per-record cost on every
    record in a corpus.
    """
    idx: dict[str, list[tuple[str, str]]] = {}
    for run, allele in _candidates(load_anchors(organism), segment, locus):
        if len(run) >= 2:
            idx.setdefault(run[1], []).append((run, allele))
    return {k: tuple(v) for k, v in idx.items()}


def guess_allele(cdr3: str, segment: str, locus: str, anchors: dict,
                 called: str = "", min_gain: int = 1, organism: str = "") -> str:
    """The allele of ``locus`` whose templated run best explains this junction end.

    ``cdr3`` is anchor-first (J reversed), as everywhere in this section.

    This is legacy's ``Cdr3Fixer.guess_id``, which recovered a segment id from the junction alone,
    rebuilt on the same scan as everything else instead of on a lookup table of CDR3 substrings.
    ⚠ Legacy's own V branch never worked: it puts ``return ""`` INSIDE the five-prime loop, so it
    tried one prefix length and gave up -- 3 non-empty V guesses in 4,000 sequences against 3,797
    for J, whose branch has the same statement correctly in a ``for ... else``.

    **This exists so a contradicted call changes the ALLELE, not the sequence.** A junction whose
    own anchor-side residues match a different allele better than the called one is evidence about
    the CALL; substituting residues to satisfy the call instead rewrites correct data.

    ``min_gain`` is how many MORE residues the winner must explain than the submitted call does
    (:data:`_RECALL_GAIN`). Returns ``""`` when nothing clears it, so the submitted call stands
    unless the junction argues against it by a margin -- confirming a call is as much the job as
    replacing one.
    """
    # ⚠ The comparison is the hit's LEADING EXACT RUN from the anchor (`_RECALL_GAIN`), and only a
    # hit that starts at the germline's OWN ANCHOR is evidence about the call. `CSAR` is
    # TRBV20-1's run and it does sit inside `CGGSARSGELFF` -- at offset 3, where it is the N region
    # reproducing three residues by chance, and TRBV20-1's junction would read `CSAR...`. Scoring
    # unanchored hits re-called that record off a coincidence. The same requirement is why a
    # junction cut INSIDE its germline (`start_in_segment > 0`, a `FixAdd`) is left with the call it
    # was submitted with: there is no anchored hit to compare, so there is nothing to compare.
    if not called or (segment, called) not in anchors:
        return ""
    if anchors[(segment, called)].status != "ok":
        return ""                        # a partial reference record cannot lose a length contest
    t = anchors[(segment, called)].templated_aa
    h = scan(t[::-1] if segment == "J" else t, cdr3)
    if h is None or h.start_in_segment != 0:
        return ""
    seg_t = t[::-1] if segment == "J" else t
    base = anchor_depth(seg_t, cdr3, h.start_in_cdr3)
    floor = base + min_gain
    # ⛔ If the submitted call already explains the junction's anchor residue, no candidate that
    # breaks it may displace it, however much more it explains behind. Otherwise the deeper allele
    # wins and step three then rewrites a junction that was canonical when it arrived:
    # `CAARLGNNYKLIW` is a clean TRAJ33 junction (`DSNYQLIW`, terminal W and all) that TRAJ12
    # outscored, whereupon its W became an F. The reverse case is untouched -- `CASSKRGGYEQYV` on
    # TRBJ2-7*01 (`SYEQYF`) does NOT explain its anchor, so *02 (`SYEQYV`) is free to take it.
    keep_anchor = bool(seg_t) and seg_t[0] == cdr3[h.start_in_cdr3:h.start_in_cdr3 + 1]
    # The offset-1 index is only a valid prefilter when the margin needs depth >= 2, which every
    # shipped setting does; a caller asking for less gets the full scan.
    q1 = cdr3[h.start_in_cdr3 + 1:h.start_in_cdr3 + 2]
    if organism and floor >= 2 and q1:
        cands = _cached_candidate_index(organism, segment, locus).get(q1, ())
    elif organism:
        cands = _cached_candidates(organism, segment, locus)
    else:
        cands = _candidates(anchors, segment, locus)
    # ⛔ Every candidate is scored AT THE CALLED ALLELE'S OWN OFFSET, not at a placement of its
    # own. Where the junction's germline run begins is a property of the JUNCTION -- a framework
    # flank in front of Cys104 is there whatever allele is named -- so scoring each allele at its
    # own best offset let one win by placing its anchor deeper inside the junction, where the trim
    # that follows then deletes everything before it: `TRBJ2-2P`'s `LRGAAG` matches `GAAG` seven
    # residues into `CASRPGAAGGRPELYF`, and winning there trimmed the record to `CASRPGAAG`. (IMGT
    # files that pseudogene as an ORF, so functionality does not catch it, and its anchor decodes
    # G, so `_canonicalise` accepts the result.)
    #
    # It is also what makes re-calling every record affordable: one `anchor_depth` per candidate,
    # which is a walk of at most the run's length, instead of a `scan` -- two `_extend` sweeps and
    # a C++ alignment. That was 1.4 M `_extend` calls over 20,000 records, 36 % of the markup.
    off = h.start_in_cdr3
    q_at = cdr3[off:off + 1]
    best, best_score = "", -1
    for t, allele in cands:
        if len(t) < base:
            continue                     # a run cannot agree on more residues than it has
        holds = t[:1] == q_at
        if keep_anchor and not holds:
            continue
        depth = anchor_depth(t, cdr3, off)
        # The margin is waived for a candidate that explains an anchor the call does NOT, at equal
        # depth: it accounts for the same residues PLUS the conserved one, so it explains strictly
        # more and asks for no edit. `CASSKRGGYEQYV` on TRBJ2-7*01 (`SYEQYF`) is that case --
        # `YQEYS` agrees either way and only *02 (`SYEQYV`) also accounts for the terminal V, which
        # is why VDJdb rewrote that V to an F to keep *01.
        if depth < (base if (holds and not keep_anchor) else floor):
            continue
        # ⚠ Functionality breaks a TIE and nothing more. A junction is evidence about which allele
        # it came from, so an ORF that explains it where no functional allele does still wins --
        # `CASSKRGGYEQYV` is a clean TRBJ2-7*02 (`SYEQYV`, ORF) junction. But where a functional
        # allele explains the junction just as well, it is the likelier rearrangement. `_candidates`
        # is ordered functional-first then by name and the comparison is strict, so which allele
        # represents a shared run does not depend on dict order.
        rank = (holds, depth, anchors[(segment, allele)].functionality == "F")
        if best_score == -1 or rank > best_score:
            best, best_score = allele, rank
    return "" if best == called else best


def _placeable(anchor: Anchor | None) -> bool:
    """Can this anchor place a junction boundary at all?

    ``status == "ok"`` always can. ``truncated`` can when it still carries
    ``_MIN_TRUNCATED_AA`` residues: IMGT ships those allele records as partial sequences that
    stop inside the anchor region, so the templated run is SHORT but correct as far as it goes,
    and refusing it throws away a boundary that is present. ``no_anchor`` never can.
    """
    if anchor is None or not anchor.templated_aa:
        return False
    if anchor.status == "ok":
        return True
    return (anchor.status == "truncated"
            and len(anchor.templated_aa) >= _MIN_TRUNCATED_AA)


# Codons per residue. `CODON_TABLE` is the standard code, nt -> aa; the boundary question
# below asks it backwards.
_CODONS_BY_AA: dict[str, list[str]] = {}
for _codon, _residue in CODON_TABLE.items():
    _CODONS_BY_AA.setdefault(_residue, []).append(_codon)

# Relative mass of a codon that agrees with the germline over exactly `e` nucleotides,
# summed over where the germline actually stopped inside that codon: sum(4**r for r in
# 0..e). A templated nucleotide is free; an inserted one that reproduces the germline by
# chance costs 1/4. Index is `e`; a boundary codon cannot reach 3, because the residue
# there differs from the germline's by construction.
_BOUNDARY_W = (1.0, 5.0, 21.0, 85.0)


def _codon_extension(residue: str, germ: str, five_prime: bool) -> tuple[int, float]:
    """Nucleotides of the boundary codon the germline still explains, and P(that many)."""
    codons = _CODONS_BY_AA.get(residue)
    if not codons or not germ:
        return 0, 1.0
    n = len(germ)
    mass = [0.0] * (n + 1)
    for x in codons:
        e = 0
        while e < n and (x[e] == germ[e] if five_prime else x[2 - e] == germ[n - 1 - e]):
            e += 1
        mass[e] += 1.0
    score = [m * _BOUNDARY_W[e] for e, m in enumerate(mass)]
    best = max(range(n + 1), key=score.__getitem__)      # ties go to the shorter run
    return best, score[best] / sum(score)


def boundary_nt(junction: str, residues: int, germline_nt: str,
                side: str) -> tuple[int, float]:
    """Where the germline stops in NUCLEOTIDES, from an amino-acid boundary.

    ``residues`` is the amino-acid answer -- ``Cdr3Markup.v_end`` (count of V-templated
    residues) or ``.j_start`` (index of the first J-templated residue). Returns the
    junction-nucleotide boundary and the probability of that value: for ``"V"`` a count,
    i.e. a half-open end, and for ``"J"`` the index of the first templated nucleotide.
    Both match the ``VEnd + 1`` / ``JStart`` an aligner reports off the observed sequence.

    **The amino acid cannot say this and the nucleotide can.** A germline run ends
    wherever the exonuclease stopped, which is not a codon boundary, so the last residue
    it touches is usually part germline and part N region -- and an alignment on the
    protein can only round that to a whole residue. What the protein DOES fix is which
    nucleotides are admissible: 17 of the 20 residues have the same first base in every
    one of their codons, so on the V side the answer is frequently forced. ``GGA`` (Gly)
    against an observed Glu can only be ``GAA``/``GAG``, both of which open with the
    germline's ``G``, so the germline demonstrably reaches one nucleotide further --
    whether it was templated or an insertion reproduced it, an aligner reads it as
    germline either way. Where the codons disagree, ``_BOUNDARY_W`` weighs them.

    Measured against ``isalgo/airr_control``'s ``human.trb.ntvj``, on the 8,133 VDJdb
    human TRB junctions whose boundary every control observation agrees on. VDJdb's
    ``v.end`` counts residues the V templates at least 2 nucleotides of, ``(nt + 1) // 3``:

    ======================  ==============  ==============
    ``v.end``               exact           within 1
    ======================  ==============  ==============
    k-mer scanner / 2.30.1  5,839 (71.8 %)  8,089 (99.5 %)
    this                    7,556 (92.9 %)  8,040 (98.9 %)
    ======================  ==============  ==============

    In nucleotides the same comparison is 28.96 % -> **80.39 %** on V and 30.27 % ->
    **74.47 %** on J. The two sides are limited by different things, and the numbers say
    which: the V-side residue count is right on 83.99 % of records and the extension on
    **95.71 %** of those, so V is bounded by the protein alignment; the J-side count is
    right on 97.97 % and the extension on **76.02 %**, because the nucleotide a J
    boundary turns on is the codon's third, the one the genetic code leaves free.

    ``j.start`` does not move: VDJdb defines it as the first FULLY J-templated residue,
    ``ceil(nt / 3)``, which is the same residue whether the germline reaches 1 or 2
    nucleotides into the one before it. 97.97 % before and after, by construction.
    """
    if residues < 0 or not germline_nt or side not in ("V", "J"):
        return -1, 0.0
    if side == "V":
        base = 3 * residues
        if base >= len(germline_nt) or residues >= len(junction):
            return min(base, len(germline_nt)), 1.0
        e, p = _codon_extension(junction[residues], germline_nt[base:base + 3], True)
        return base + e, p
    # J: `avail` is how far the germline reaches back past the templated run, so the
    # boundary codon is the three germline nucleotides ending there.
    base = 3 * residues
    avail = len(germline_nt) - 3 * (len(junction) - residues)
    if avail <= 0 or residues <= 0:
        return base, 1.0
    e, p = _codon_extension(junction[residues - 1], germline_nt[max(0, avail - 3):avail], False)
    return base - e, p


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

    # A call may name several candidates the curator could not separate -- comma-joined, or
    # `+`-joined as mouse `TRBV12-2+TRBV13-2` (532 keys) is. The leading one wins, which is what
    # this module has always done with a comma.
    v_sub, j_sub = _leading_call(v_call), _leading_call(j_call)
    v_id = resolve_allele(v_sub, "V", anchors, organism)
    j_id = resolve_allele(j_sub, "J", anchors, organism)
    rec.v_call, rec.j_call = v_id, j_id
    v_moved, j_moved = _allele_moved(v_sub, v_id), _allele_moved(j_sub, j_id)
    v_anchor = anchors.get(("V", v_id))
    j_anchor = anchors.get(("J", j_id))

    # Each side runs anchor-first -- the V as submitted, the J reversed -- so offset 0 is the
    # conserved residue in both and legacy's table reads the same way at either end. The J side
    # runs on the V-repaired junction, exactly as `Cdr3Fixer.fix_both` does.
    repaired = cdr3
    for side in ("V", "J"):
        anchor = v_anchor if side == "V" else j_anchor
        moved = v_moved if side == "V" else j_moved
        if not _placeable(anchor):
            flags, name, templated = _with_allele(("impossible",), moved), "FailedBadSegment", -1
        else:
            fwd = repaired if side == "V" else repaired[::-1]

            # ⛔ STEP ONE, before any edit: does another allele of this locus explain the submitted
            # residues better? A call the junction CONTRADICTS is a call to RE-ASSIGN, never a
            # sequence to rewrite -- and the allele that wins is the one the trim, the addition and
            # the anchor are all then measured against, so it has to be chosen first.
            #
            # Measured over VDJdb's corpus, 74.7 % of anchor-adjacent substitutions into an
            # ALREADY-canonical junction were records where another allele fit the submitted 3' end
            # better than the called one -- against 1.1 % of untouched records, a 68x enrichment.
            # `CASSLRGAATDTQYF` is the shape of it: a clean TRBJ2-3 junction called TRBJ2-1, which
            # 2.16.0-2.31.0 "repaired" to `CASSLRGAATDTQFF`, a string no germline supports.
            better = guess_allele(fwd, side, anchor.locus, anchors,
                                  called=v_id if side == "V" else j_id,
                                  min_gain=_RECALL_GAIN[side], organism=organism)
            if better:
                anchor = anchors[(side, better)]
                moved = True
                if side == "V":
                    rec.v_call = v_id = better
                else:
                    rec.j_call = j_id = better
            seg = anchor.templated_aa if side == "V" else anchor.templated_aa[::-1]
            fixed, flags, name, templated, errs = _fix_side(seg, fwd, side, max_replace)
            flags = _with_allele(flags, moved)
            repaired = fixed if side == "V" else fixed[::-1]
            rec.errors.extend(errs)
            # A germline record IMGT ships incomplete still places a boundary, but a LOWER-bound
            # one: residues past the record's end are unattributed, not known non-templated.
            if anchor.status == "truncated" and name in _GOOD:
                name = "TruncatedGermline"
        if side == "V":
            v_anchor = anchor
            rec.v_flags, rec.v_fix = flags, name
            rec.v_end = templated
        else:
            j_anchor = anchor
            rec.j_flags, rec.j_fix = flags, name
            rec.j_start = (len(repaired) - templated) if templated >= 0 else -1

    rec.cdr3_repaired = repaired
    # Read off the REPAIRED junction, which is what `v_end` / `j_start` index.
    if v_anchor is not None and rec.v_end >= 0:
        rec.v_end_nt = boundary_nt(repaired, rec.v_end, v_anchor.germline_nt, "V")[0]
    if j_anchor is not None and rec.j_start >= 0:
        rec.j_start_nt = boundary_nt(repaired, rec.j_start, j_anchor.germline_nt, "J")[0]
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
        "v_end": m.v_end, "j_start": m.j_start,
        "v_end_nt": m.v_end_nt, "j_start_nt": m.j_start_nt,
        "v_fix": m.v_fix, "j_fix": m.j_fix,
        "v_flags": ",".join(m.v_flags), "j_flags": ",".join(m.j_flags),
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
