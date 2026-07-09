"""Per-segment AIRR CIGAR strings from the mmseqs scaffold alignment.

arda aligns a query to a `V + N*pad + J [+ C]` scaffold, not to each germline segment separately.
AIRR wants a CIGAR per segment (``v_cigar``/``j_cigar``/``c_cigar``) whose *reference* is that
segment's germline; since the scaffold's V part IS the V germline (target position == germline
position), its J part IS the J allele, and its C part the CH1 exon, each segment's CIGAR is the
sub-walk of the one query->scaffold alignment whose target falls in that segment's range.

CIGAR operators follow the AIRR spec (SAM subset):
  * ``S`` -- query positions before the alignment starts (query 5' offset). Required, precedes N.
  * ``N`` -- reference positions before the alignment starts (germline 5' offset). Required.
  * ``M`` / ``I`` (gap in reference) / ``D`` (gap in query) -- the aligned body.
  * trailing ``S`` (query 3' remainder) is emitted; trailing ``N`` (germline 3' remainder) is
    optional per the spec and omitted (arda does not always know the full germline length, e.g. the
    C-region CH1 exon is longer than the shipped stub).

``segment_cigars`` builds all three in a SINGLE pass over the aligned strings.

Correcting cigars for CONTIGS (Stage 3). A contig is just a long query, so BOTH ways to get its
cigars end in ``segment_cigars`` and produce the same record (see :mod:`arda.annotate.contig`):

  * RE-ANNOTATE the assembled contig through ``mapper.annotate_records`` -- one mmseqs alignment,
    then ``segment_cigars``. No cigar arithmetic; ``check_cigar`` validates it.
  * MERGE the reads' existing alignments column-by-column into the contig's (C++
    ``_markup.merge_alignment``), skipping the alignment pass.

Both are built and proven byte-for-byte equal (``tests/unit/test_contig_merge.py`` on 29 real
GenBank contigs). MEASURED (arda-benchmark ``scripts/bench_contig_cigars.py``): at scRNA-seq scale
(~10^5 contigs/sample) merge is ~9x faster -- the whole gap is mmseqs; the C++ stitch is ~3 % of
merge's wall and barely grows with read depth. Prefer merge when the assembly layout is available
(the reads carry their scaffold + offset); re-annotate is the fallback when it is not.
"""

from __future__ import annotations

import re

_CIGAR_RE = re.compile(r"(\d+)([MIDNS=X])")


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """``"57S291M1054S"`` -> ``[(57,"S"), (291,"M"), (1054,"S")]``. Inverse of :func:`build_cigar`."""
    return [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]


def cigar_query_length(cigar: str) -> int:
    """Query (read/contig) bases the CIGAR spans -- M/I/S/=/X; D and N are reference-side."""
    return sum(n for n, op in parse_cigar(cigar) if op in "MIS=X")


def cigar_reference_length(cigar: str) -> int:
    """Reference (germline) bases the CIGAR spans -- M/D/N/=/X; I and S are query-side."""
    return sum(n for n, op in parse_cigar(cigar) if op in "MDN=X")


def check_cigar(cigar: str, query_len: int) -> bool:
    """A CIGAR is consistent with a query of ``query_len`` iff its query-side ops sum to it.

    This is the invariant a corrected/re-annotated sequence (a read OR an assembled contig -- a
    contig is just a long query) must satisfy: ``v_cigar``/``j_cigar``/``c_cigar`` each lay over the
    WHOLE sequence, soft-clipping the parts outside their own segment. Use it to validate a cigar
    after correcting or re-deriving it.
    """
    return cigar_query_length(cigar) == query_len


def _classify(tpos: int, t_vend: int, t_jstart: int, t_vjend: int) -> str | None:
    """Which germline segment target position ``tpos`` (1-based, on the scaffold) belongs to.

    V is ``[1, t_vend]``, J is ``[t_jstart, t_vjend]``, C is ``> t_vjend``. The N-pad between the V
    end and the J start belongs to no germline -- it is the np region. A zero boundary means that
    segment is absent (a `J + C` scaffold has ``t_vend == 0``).
    """
    if t_vend and tpos <= t_vend:
        return "v"
    if t_jstart and t_vjend and t_jstart <= tpos <= t_vjend:
        return "j"
    if t_vjend and tpos > t_vjend:
        return "c"
    return None


def build_cigar(q_lead: int, g_lead: int, ops: list[str], q_trail: int) -> str:
    """Assemble one AIRR CIGAR: ``{q_lead}S {g_lead}N <body> {q_trail}S`` (parts of length 0 are
    dropped). ``ops`` is the per-column M/I/D sequence of the aligned body; consecutive equal
    operators are run-length encoded. Trailing germline ``N`` is intentionally omitted (optional)."""
    runs: list[list] = []
    for op in ops:
        if runs and runs[-1][1] == op:
            runs[-1][0] += 1
        else:
            runs.append([1, op])
    body = "".join(f"{n}{op}" for n, op in runs)
    return ((f"{q_lead}S" if q_lead > 0 else "")
            + (f"{g_lead}N" if g_lead > 0 else "")
            + body
            + (f"{q_trail}S" if q_trail > 0 else ""))


def _germline_pos(key: str, t: int, t_jstart: int, t_vjend: int) -> int:
    """1-based position of scaffold target ``t`` within its own germline (V: ==t; J and C: offset)."""
    if key == "v":
        return t
    if key == "j":
        return t - t_jstart + 1
    return t - t_vjend                                    # C: germline starts one past the V-J end


def segment_cigars(qaln: str, taln: str, qstart: int, tstart: int, qlen: int,
                   t_vend: int, t_jstart: int, t_vjend: int) -> dict[str, str]:
    """Return ``{"v_cigar":…, "j_cigar":…, "c_cigar":…}`` (only the segments that have a body).

    ``qaln``/``taln`` are the mmseqs aligned strings (``-`` for gaps), ``qstart``/``tstart`` their
    1-based query/target start, ``qlen`` the full query length. Boundaries are 1-based scaffold
    positions; pass 0 for an absent segment.
    """
    body: dict[str, list[str]] = {"v": [], "j": [], "c": []}
    q_first: dict[str, int] = {}
    q_last: dict[str, int] = {}
    g_first: dict[str, int] = {}
    q, t = qstart, tstart
    for qa, ta in zip(qaln, taln):
        cq, ct = qa != "-", ta != "-"
        op = "M" if (cq and ct) else ("I" if cq else "D")
        key = _classify(t, t_vend, t_jstart, t_vjend)     # for an insertion, t is the next base
        if key:
            body[key].append(op)
            if cq:
                q_first.setdefault(key, q)
                q_last[key] = q
            if ct:
                g_first.setdefault(key, _germline_pos(key, t, t_jstart, t_vjend))
        if cq:
            q += 1
        if ct:
            t += 1

    out = {}
    for key, ops in body.items():
        if key not in q_first:                            # no query aligned to this segment
            continue
        cig = build_cigar(q_first[key] - 1, g_first.get(key, 1) - 1, ops, qlen - q_last[key])
        if cig:
            out[f"{key}_cigar"] = cig
    return out
