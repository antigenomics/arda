"""Project reference region markup onto a query via the C++ hot path.

Takes a parsed mmseqs hit plus the reference entry for the matched scaffold and
returns an AIRR-style record dict for the query.
"""

from __future__ import annotations

from .. import _markup
from ..refbuild.translate import translate, aa_coords_from_nt, detect_coding_frame
from .reference import RefEntry, REGIONS

__all__ = ["transfer_hit", "AIRR_COLUMNS"]

# Output column order (AIRR-compatible subset + locus).
AIRR_COLUMNS = (
    ["sequence_id", "sequence", "locus", "v_call", "j_call", "rev_comp", "productive",
     "v_sequence_start", "v_sequence_end", "junction", "junction_aa"]
    + [c for r in REGIONS for c in (f"{r}_start", f"{r}_end", r, f"{r}_aa")]
)


def _empty_record(query_id: str, query_seq: str) -> dict:
    rec = {c: "" for c in AIRR_COLUMNS}
    rec["sequence_id"] = query_id
    rec["sequence"] = query_seq
    return rec


def transfer_hit(
    query_id: str,
    query_seq: str,
    hit: dict,
    ref: RefEntry,
    seqtype: str = "nt",
) -> dict:
    """Build an AIRR record by projecting ``ref`` region coords onto the query.

    Args:
        query_id: Query sequence id.
        query_seq: Full query sequence (nt or aa per ``seqtype``).
        hit: Parsed mmseqs row with qstart/qend/tstart/tend/qaln/taln.
        ref: Reference entry for the matched scaffold.
        seqtype: ``"nt"`` or ``"aa"``.
    """
    coords = _markup.transfer_regions(
        hit["qaln"], hit["taln"],
        int(hit["qstart"]), int(hit["tstart"]),
        ref.starts, ref.ends,
    )

    rec = _empty_record(query_id, query_seq)
    rec.update(locus=ref.locus, v_call=ref.v_call, j_call=ref.j_call,
               rev_comp="F", productive="")

    region_q: dict[str, tuple[int, int]] = {}
    for name, (qs, qe) in zip(REGIONS, coords):
        if qs < 0:
            continue
        region_q[name] = (qs, qe)
        rec[f"{name}_start"] = qs
        rec[f"{name}_end"] = qe
        rec[name] = query_seq[qs - 1 : qe]

    # V span (FR1 start .. FR3 end) and coding frame origin.
    if "fwr1" in region_q:
        rec["v_sequence_start"] = region_q["fwr1"][0]
    if "fwr3" in region_q:
        rec["v_sequence_end"] = region_q["fwr3"][1]

    coding_start = region_q.get("fwr1", (None, None))[0]

    # Amino acid regions + productivity.
    if seqtype == "nt":
        if coding_start is not None:
            # The alignment may start mid-codon, so the projected FR1 start is not
            # necessarily a codon boundary. Recover the true reading frame from the
            # V span (exactly one frame is stop-free for a real V), then translate.
            v_end = region_q.get("fwr3", (0, 0))[1] or region_q.get("cdr3", (0, 0))[0]
            v_end = v_end or min(len(query_seq), coding_start + 200)
            coding_start += detect_coding_frame(query_seq[coding_start - 1 : v_end])
            protein = translate(query_seq[coding_start - 1:], 0)
            for name, (qs, qe) in region_q.items():
                a_s, a_e = aa_coords_from_nt(qs, qe, coding_start)
                a_s = max(1, a_s)
                rec[f"{name}_aa"] = protein[a_s - 1 : a_e]
            rec["productive"] = "T" if "*" not in protein else "F"
    else:  # aa input: region seqs already amino acids
        for name, (qs, qe) in region_q.items():
            rec[f"{name}_aa"] = query_seq[qs - 1 : qe]
        if "fwr1" in region_q and "fwr4" in region_q:
            span = query_seq[region_q["fwr1"][0] - 1 : region_q["fwr4"][1]]
            rec["productive"] = "T" if "*" not in span else "F"

    # Junction = CDR3 extended by one conserved codon/residue each side.
    if "cdr3" in region_q:
        cs, ce = region_q["cdr3"]
        flank = 3 if seqtype == "nt" else 1
        js, je = cs - flank, ce + flank
        if js >= 1 and je <= len(query_seq):
            rec["junction"] = query_seq[js - 1 : je]
            rec["junction_aa"] = (
                translate(rec["junction"], 0) if seqtype == "nt" else rec["junction"]
            )
    return rec
