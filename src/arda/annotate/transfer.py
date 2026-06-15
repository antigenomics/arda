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

# Output column order (AIRR-compatible subset + locus). ``mmseqs2_score``/``_evalue``/
# ``_identity`` carry the alignment quality of the chosen scaffold hit. The score is the
# mmseqs2 bit score over the *whole* V+J scaffold (not a per-segment AIRR ``v_score``), so
# it is named after its source; callers use it to rank references and filter weak hits.
AIRR_COLUMNS = (
    ["sequence_id", "sequence", "locus", "v_call", "d_call", "d2_call", "j_call",
     "mmseqs2_score", "mmseqs2_evalue", "mmseqs2_identity",
     "rev_comp", "productive",
     "v_sequence_start", "v_sequence_end",
     "d_sequence_start", "d_sequence_end", "d2_sequence_start", "d2_sequence_end",
     "j_sequence_start", "np1", "np2", "np3", "junction", "junction_aa"]
    + [c for r in REGIONS for c in (f"{r}_start", f"{r}_end", r, f"{r}_aa")]
)

_VSIDE = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

# D-segment mapping. D germlines are short (~8-31 nt) and trimmed on both ends, so
# they are mapped by gapless local alignment (``_markup.d_local_align``, net score
# match=+1/mismatch=-1) of every locus D allele against the V..J interior of the
# junction — not via the mmseqs scaffold DB. Thresholds below balance sensitivity
# against chance matches in a short interior.
_D_MIN_SCORE = 6        # minimum net score to call a D segment
_D2_MIN_SCORE = 7       # stricter threshold for a second (D-D fusion) segment
_DD_LOCI = {"IGH", "TRD"}   # loci where tandem D-D fusions occur


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


def _best_d(interior, d_germlines, min_score, exclude=None):
    """Best-aligning D against ``interior`` as ``(score, length, allele, s, e)``.

    ``s``/``e`` are 0-based inclusive offsets within ``interior``. ``exclude`` is an
    optional ``(s, e)`` span the match must not overlap (used to find a second,
    non-overlapping D). Returns ``None`` if nothing scores at least ``min_score``.
    """
    best = None
    for allele, dseq in d_germlines:
        score, s, e = _markup.d_local_align(interior, dseq)
        if score < min_score or s < 0:
            continue
        if exclude is not None:
            xs, xe = exclude
            if not (e < xs or s > xe):       # overlaps the excluded span
                continue
        cand = (score, e - s + 1, allele, s, e)
        if best is None or cand[:2] > best[:2] or (cand[:2] == best[:2] and allele < best[2]):
            best = cand
    return best


def _map_d(rec, query_seq, locus, v_end_q, j_start_q, d_germlines):
    """Map D segment(s) into the V..J interior and populate d_call/np regions.

    Coordinates emitted are AIRR (1-based closed, query space). For D-D loci a
    second non-overlapping D is sought; the two are then ordered 5'->3' as
    ``d_call`` / ``d2_call`` with ``np1``/``np2``/``np3`` between V, the D(s), and J.
    """
    if not d_germlines or not v_end_q or not j_start_q:
        return
    i_lo, i_hi = v_end_q + 1, j_start_q - 1       # 1-based interior bounds (query)
    if i_hi < i_lo:
        return
    interior = query_seq[i_lo - 1 : i_hi]
    d1 = _best_d(interior, d_germlines, _D_MIN_SCORE)
    if d1 is None:
        return

    segs = [d1]
    if locus in _DD_LOCI:
        d2 = _best_d(interior, d_germlines, _D2_MIN_SCORE, exclude=(d1[3], d1[4]))
        if d2 is not None:
            segs.append(d2)
    segs.sort(key=lambda c: c[3])                 # order 5'->3' by interior start

    def q(off):                                   # interior 0-based offset -> query 1-based
        return i_lo + off

    _, _, a1, s1, e1 = segs[0]
    rec["d_call"] = a1
    rec["d_sequence_start"], rec["d_sequence_end"] = q(s1), q(e1)
    if len(segs) == 2:
        _, _, a2, s2, e2 = segs[1]
        rec["d2_call"] = a2
        rec["d2_sequence_start"], rec["d2_sequence_end"] = q(s2), q(e2)
        rec["np1"] = query_seq[v_end_q : q(s1) - 1]
        rec["np2"] = query_seq[q(e1) : q(s2) - 1]
        rec["np3"] = query_seq[q(e2) : j_start_q - 1]
    else:
        rec["np1"] = query_seq[v_end_q : q(s1) - 1]
        rec["np2"] = query_seq[q(e1) : j_start_q - 1]


def transfer_hit(
    query_id: str,
    query_seq: str,
    hit: dict,
    ref: RefEntry,
    seqtype: str = "nt",
    rev_comp: bool = False,
    d_germlines: list[tuple[str, str]] | None = None,
) -> dict:
    """Build an AIRR record by projecting ``ref`` region coords onto the query."""
    coords = _markup.transfer_regions(
        hit["qaln"], hit["taln"], int(hit["qstart"]), int(hit["tstart"]),
        ref.starts, ref.ends)

    rec = _empty_record(query_id, query_seq)
    rec.update(locus=ref.locus, v_call=ref.v_call, j_call=ref.j_call,
               rev_comp="T" if rev_comp else "F", productive="")

    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return ""

    # Alignment quality of the scaffold hit (mmseqs2 bit score / E-value / % identity).
    rec["mmseqs2_score"] = _num(hit.get("bits"))
    rec["mmseqs2_evalue"] = _num(hit.get("evalue"))
    rec["mmseqs2_identity"] = _num(hit.get("pident"))

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
        # D-segment mapping (VDJ loci only; gated by presence of D germlines).
        _map_d(rec, query_seq, ref.locus, v_end_q, j_start_q, d_germlines)
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
