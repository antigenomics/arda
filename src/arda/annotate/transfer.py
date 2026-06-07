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
    rev_comp: bool = False,
) -> dict:
    """Build an AIRR record by projecting ``ref`` region coords onto the query.

    Args:
        query_id: Query sequence id.
        query_seq: Full query sequence (nt or aa per ``seqtype``). For a
            reverse-strand nt hit this is already the reverse complement (coding
            strand), with ``hit`` coordinates remapped onto it by the caller.
        hit: Parsed mmseqs row with qstart/qend/tstart/tend/qaln/taln.
        ref: Reference entry for the matched scaffold.
        seqtype: ``"nt"`` or ``"aa"``.
        rev_comp: Whether the input was reverse-complemented to align.
    """
    coords = _markup.transfer_regions(
        hit["qaln"], hit["taln"],
        int(hit["qstart"]), int(hit["tstart"]),
        ref.starts, ref.ends,
    )

    rec = _empty_record(query_id, query_seq)
    rec.update(locus=ref.locus, v_call=ref.v_call, j_call=ref.j_call,
               rev_comp="T" if rev_comp else "F", productive="")

    region_q: dict[str, tuple[int, int]] = {}
    for name, (qs, qe) in zip(REGIONS, coords):
        if qs < 0:
            continue
        region_q[name] = (qs, qe)
        rec[f"{name}_start"] = qs
        rec[f"{name}_end"] = qe
        rec[name] = query_seq[qs - 1 : qe]

    # CDR3 length is query-specific (somatic D / N insertions), so its end CANNOT
    # come from the fixed-length scaffold (V-end + short N spacer + J-start). Anchor
    # the CDR3 end to FR4: CDR3 runs from the V-anchored start up to just before the
    # conserved [FW] that opens FR4. This makes long CDR3s come out at full length.
    if "cdr3" in region_q and "fwr4" in region_q:
        cs = region_q["cdr3"][0]
        ce = region_q["fwr4"][0] - 1
        if ce >= cs:
            region_q["cdr3"] = (cs, ce)
            rec["cdr3_start"], rec["cdr3_end"], rec["cdr3"] = cs, ce, query_seq[cs - 1 : ce]

    # V span (FR1 start .. FR3 end) and coding frame origin.
    if "fwr1" in region_q:
        rec["v_sequence_start"] = region_q["fwr1"][0]
    if "fwr3" in region_q:
        rec["v_sequence_end"] = region_q["fwr3"][1]

    coding_start = region_q.get("fwr1", (None, None))[0]

    # Amino acid regions + productivity.
    if seqtype == "nt":
        # Derive the query coding frame from the ALIGNMENT PHASE: the scaffold reads
        # in frame 0 from its FR1 target start (t0), so the query position aligning
        # to the first in-phase target codon boundary is the frame origin. This works
        # even when FR1 is absent (e.g. a CDR3-FR4 / J-side fragment).
        t0 = ref.starts[0]
        tstart = int(hit["tstart"])
        p = tstart + ((t0 - tstart) % 3)  # first pos >= tstart with (p-t0)%3==0
        pj = _markup.transfer_regions(
            hit["qaln"], hit["taln"], int(hit["qstart"]), tstart, [p], [p])[0][0]
        if pj > 0:
            coding_start = pj
        elif coding_start is not None:  # boundary in a gap; stop-free fallback
            v_end = region_q.get("fwr3", (0, 0))[1] or len(query_seq)
            coding_start += detect_coding_frame(query_seq[coding_start - 1 : v_end])
        if coding_start is not None:
            protein = translate(query_seq[coding_start - 1:], 0)
            cdr3_aa_start = None
            for name, (qs, qe) in region_q.items():
                a_s, a_e = aa_coords_from_nt(qs, qe, coding_start)
                a_s = max(1, a_s)
                rec[f"{name}_aa"] = protein[a_s - 1 : a_e]
                if name == "cdr3":
                    cdr3_aa_start = a_s
            rec["productive"] = "T" if "*" not in protein else "F"
            _set_junction(rec, query_seq, region_q, seqtype="nt",
                          protein=protein, cdr3_aa_start=cdr3_aa_start)
    else:  # aa input: region seqs already amino acids
        for name, (qs, qe) in region_q.items():
            rec[f"{name}_aa"] = query_seq[qs - 1 : qe]
        if "fwr1" in region_q and "fwr4" in region_q:
            span = query_seq[region_q["fwr1"][0] - 1 : region_q["fwr4"][1]]
            rec["productive"] = "T" if "*" not in span else "F"
        _set_junction(rec, query_seq, region_q, seqtype="aa")
    return rec


def _set_junction(rec, query_seq, region_q, *, seqtype, protein=None, cdr3_aa_start=None):
    """Set junction / junction_aa with the exact AIRR semantics.

    junction spans the conserved Cys104 through the [FW]118 that opens FR4. The
    nucleotide ``junction`` is sliced from the query; ``junction_aa`` is built as
    ``Cys + cdr3_aa + FW`` so that ``cdr3_aa == junction_aa[1:-1]`` holds *by
    construction* — never off by one — even for non-productive sequences.
    """
    if "cdr3" not in region_q or "fwr4" not in region_q:
        return
    cs = region_q["cdr3"][0]
    f4 = region_q["fwr4"][0]
    flank = 3 if seqtype == "nt" else 1
    js, je = cs - flank, f4 + (flank - 1)
    cdr3_aa = rec.get("cdr3_aa", "")
    fw = rec.get("fwr4_aa", "")[:1]
    if seqtype == "nt":
        cys = protein[cdr3_aa_start - 2] if (cdr3_aa_start and cdr3_aa_start >= 2) else ""
    else:
        cys = query_seq[cs - 2] if cs >= 2 else ""
    # A valid AIRR junction needs BOTH conserved flanks (Cys104 and [FW]118). If
    # the query is truncated past the Cys, emit no junction rather than a partial
    # one — this keeps the invariant cdr3_aa == junction_aa[1:-1] always true.
    if not cys or not fw or js < 1 or je > len(query_seq):
        return
    rec["junction"] = query_seq[js - 1 : je]
    rec["junction_aa"] = cys + cdr3_aa + fw
