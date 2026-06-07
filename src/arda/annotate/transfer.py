"""Project reference region markup onto a query via the C++ hot path.

Takes a parsed mmseqs hit plus the reference entry for the matched scaffold and
returns an AIRR-style record dict for the query.

Junction handling follows AIRR strictly: ``junction`` spans the conserved Cys104
through the [FW]118 that opens FR4; ``junction_aa`` starts with C and ends with
F/W for a canonical rearrangement. A junction is reported **even when not
canonical** (out-of-frame, missing the conserved residues). For an out-of-frame
junction (V and J in different frames) the amino-acid translation inserts 1-2 N
bases after the V germline end to restore the J frame; the codon that then
contains an inserted N is rendered as ``_``. The V/J split inside the junction is
located from the transferred ``v_sequence_end`` / ``j_sequence_start``.
"""

from __future__ import annotations

from .. import _markup
from ..refbuild.translate import translate, aa_coords_from_nt, detect_coding_frame
from .reference import RefEntry, REGIONS

__all__ = ["transfer_hit", "AIRR_COLUMNS"]

# Output column order (AIRR-compatible subset + locus).
AIRR_COLUMNS = (
    ["sequence_id", "sequence", "locus", "v_call", "j_call", "rev_comp", "productive",
     "v_sequence_start", "v_sequence_end", "j_sequence_start", "junction", "junction_aa"]
    + [c for r in REGIONS for c in (f"{r}_start", f"{r}_end", r, f"{r}_aa")]
)

_VSIDE = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")


def _empty_record(query_id: str, query_seq: str) -> dict:
    rec = {c: "" for c in AIRR_COLUMNS}
    rec["sequence_id"] = query_id
    rec["sequence"] = query_seq
    return rec


def _project_point(hit: dict, ref_pos: int) -> int:
    """Project a single scaffold (target) nt position onto query coords (or 0)."""
    if ref_pos <= 0:
        return 0
    qs = _markup.transfer_regions(
        hit["qaln"], hit["taln"], int(hit["qstart"]), int(hit["tstart"]),
        [ref_pos], [ref_pos])[0][0]
    return qs if qs > 0 else 0


def _junction_nt(query_seq, cs, f4, coding_start, v_end_q):
    """Build (junction_nt, junction_aa, cdr3_aa, phase) for nucleotide input.

    ``cs`` = CDR3 start, ``f4`` = FR4 start (query, 1-based). The nucleotide
    junction is the real query slice (no synthetic bases). For translation, if the
    rearrangement is out of frame, 1-2 N are inserted after the V germline end so
    the J side reads correctly; the codon containing inserted N becomes ``_``.
    """
    js, je = cs - 3, f4 + 2            # Cys104 codon start .. [FW]118 codon end
    if js < 1 or je > len(query_seq):
        return "", "", "", None
    junction_nt = query_seq[js - 1 : je]
    phase = (f4 - coding_start) % 3    # 0 => J in V frame (in-frame rearrangement)
    k = (3 - phase) % 3
    if k:
        # Insert N after the V germline end, but keep the bridge strictly inside
        # CDR3 so the conserved Cys (first codon) and [FW] (last codon) are
        # preserved; fall back to just before the [FW] codon if V end is unknown.
        ins_at = (v_end_q - js + 1) if v_end_q else (len(junction_nt) - 3)
        lo, hi = 3, len(junction_nt) - 3            # after Cys codon .. before [FW] codon
        ins_at = max(0, hi) if hi < lo else min(max(ins_at, lo), hi)
        corrected = junction_nt[:ins_at] + "N" * k + junction_nt[ins_at:]
        junction_aa = list(translate(corrected, 0))
        for ci in range(ins_at // 3, (ins_at + k - 1) // 3 + 1):   # codons holding N
            if ci < len(junction_aa):
                junction_aa[ci] = "_"
        junction_aa = "".join(junction_aa)
    else:
        junction_aa = translate(junction_nt, 0)
    cdr3_aa = junction_aa[1:-1]
    return junction_nt, junction_aa, cdr3_aa, phase


def transfer_hit(
    query_id: str,
    query_seq: str,
    hit: dict,
    ref: RefEntry,
    seqtype: str = "nt",
    rev_comp: bool = False,
) -> dict:
    """Build an AIRR record by projecting ``ref`` region coords onto the query."""
    coords = _markup.transfer_regions(
        hit["qaln"], hit["taln"], int(hit["qstart"]), int(hit["tstart"]),
        ref.starts, ref.ends)

    rec = _empty_record(query_id, query_seq)
    rec.update(locus=ref.locus, v_call=ref.v_call, j_call=ref.j_call,
               rev_comp="T" if rev_comp else "F", productive="")

    region_q: dict[str, tuple[int, int]] = {}
    for name, (qs, qe) in zip(REGIONS, coords):
        if qs < 0:
            continue
        region_q[name] = (qs, qe)
        rec[f"{name}_start"], rec[f"{name}_end"], rec[name] = qs, qe, query_seq[qs - 1 : qe]

    # CDR3 end is J-anchored (somatic length is query-specific): from the V-anchored
    # start up to just before the [FW] that opens FR4.
    if "cdr3" in region_q and "fwr4" in region_q:
        cs, ce = region_q["cdr3"][0], region_q["fwr4"][0] - 1
        if ce >= cs:
            region_q["cdr3"] = (cs, ce)
            rec["cdr3_start"], rec["cdr3_end"], rec["cdr3"] = cs, ce, query_seq[cs - 1 : ce]

    # Transfer the V germline end and J germline start (extended scaffold markup).
    v_end_q = _project_point(hit, ref.v_sequence_end)
    j_start_q = _project_point(hit, ref.j_sequence_start)
    if "fwr1" in region_q:
        rec["v_sequence_start"] = region_q["fwr1"][0]
    if v_end_q:
        rec["v_sequence_end"] = v_end_q
    if j_start_q:
        rec["j_sequence_start"] = j_start_q

    if seqtype == "nt":
        # V coding frame from the alignment phase (works even without FR1).
        t0 = ref.starts[0]
        tstart = int(hit["tstart"])
        p = tstart + ((t0 - tstart) % 3)
        pj = _project_point(hit, p)
        coding_start = pj or region_q.get("fwr1", (None,))[0]
        if pj == 0 and coding_start is not None:
            v_end = region_q.get("fwr3", (0, 0))[1] or len(query_seq)
            coding_start += detect_coding_frame(query_seq[coding_start - 1 : v_end])
        if coding_start is not None:
            protein = translate(query_seq[coding_start - 1:], 0)
            for name, (qs, qe) in region_q.items():       # V-side aa from V frame
                a_s, a_e = aa_coords_from_nt(qs, qe, coding_start)
                rec[f"{name}_aa"] = protein[max(1, a_s) - 1 : a_e]
            # FR4 reads in its own (J) frame regardless of productivity.
            if "fwr4" in region_q:
                f4s, f4e = region_q["fwr4"]
                rec["fwr4_aa"] = translate(query_seq[f4s - 1 : f4e], 0)
            # Junction (+ CDR3 aa) with out-of-frame N-bridging.
            phase = None
            if "cdr3" in region_q and "fwr4" in region_q:
                jnt, jaa, c3aa, phase = _junction_nt(
                    query_seq, region_q["cdr3"][0], region_q["fwr4"][0],
                    coding_start, v_end_q)
                if jaa:
                    rec["junction"], rec["junction_aa"], rec["cdr3_aa"] = jnt, jaa, c3aa
            vclean = all("*" not in rec.get(f"{r}_aa", "") for r in _VSIDE)
            jclean = "*" not in rec.get("junction_aa", "") and "_" not in rec.get("junction_aa", "")
            rec["productive"] = "T" if (phase == 0 and vclean and jclean) else "F"
    else:  # aa input: regions are already amino acids; no frame bridging needed.
        for name, (qs, qe) in region_q.items():
            rec[f"{name}_aa"] = query_seq[qs - 1 : qe]
        if "cdr3" in region_q and "fwr4" in region_q:
            cs, f4 = region_q["cdr3"][0], region_q["fwr4"][0]
            if cs >= 2 and f4 <= len(query_seq):
                rec["junction"] = query_seq[cs - 2 : f4]
                rec["junction_aa"] = rec["junction"]
                rec["cdr3_aa"] = rec["junction"][1:-1]
        if "fwr1" in region_q and "fwr4" in region_q:
            span = query_seq[region_q["fwr1"][0] - 1 : region_q["fwr4"][1]]
            rec["productive"] = "T" if "*" not in span else "F"
    return rec
