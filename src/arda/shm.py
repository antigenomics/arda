"""Keep the SHM fields OUT of the junction.

``v_mutations`` / ``j_mutations`` / ``v_identity`` are scoped BY SEGMENT — a position counts as V
if it falls in the scaffold's V germline range and as J if it falls in the J one. That is not the
same as being outside the junction. A rearranged junction is ``V 3' tail + N/P + J 5' head``, so
the templated tails of both germlines are *inside* it, and every chew-back or non-templated base
there reads as a substitution against a germline that does not template it.

Measured on SRR5233636 — a **TRA** amplicon, where TCRs do not hypermutate, so every entry is
spurious by construction: **1.046 V and 1.658 J entries per read**, with **86.2 %** of the J
entries at J germline position ≤ 10. The high-frequency ones are junction-internal too
(``TRAV8-6*01`` 281/282 at 0.88 against its anchor at 270; ``TRAJ8*01`` position 1 at 0.67 against
its anchor at 26) — an allele difference in the templated tail, not somatic mutation.

The fix needs exactly the two numbers arda already emits per read:

* ``v_anchor_nt`` — 0-based offset of the Cys104 codon in the called V allele's germline
* ``j_anchor_nt`` — 0-based offset of the [FW]118 codon in the called J allele's germline

so a ``v_mutations`` entry at 1-based position *p* is junction-internal iff ``p > v_anchor_nt``,
and a ``j_mutations`` entry iff ``p <= j_anchor_nt + 3``.

Everything here is a pure function of one AIRR record, which is why the same code serves both
``transfer_hit`` (at emission) and the standalone ``arda shm`` stage (recounting a TSV that was
written before this existed) — there is one rule, not two implementations of it.

⚠ IGH/IGK/IGL are where SHM is real. Scoping is applied to every locus anyway, because the defect
is not IG-specific; on the TR loci the entries that survive are allele mismatches in the templated
framework, not hypermutation.

⚠ A read whose called allele has no usable anchor (the reference flags it, or the call is a tie
list whose first entry has none) is left **unscoped** rather than silently emptied. That is the
honest failure: arda does not know where the junction starts in that germline.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["SHM_MODES", "FULL_COLUMNS", "scope_record", "recount_airr"]

#: ``framework`` = scope in place (the default). ``both`` = scope in place AND keep the old,
#: junction-inclusive values in :data:`FULL_COLUMNS`. ``off`` = emit no SHM fields at all.
SHM_MODES = ("framework", "both", "off")

#: The legacy, junction-INCLUSIVE values, emitted only under ``--shm both``. They are appended
#: after every shipped column (the ``junction_quality`` rule), so a consumer reading the shipped
#: set by position is unaffected.
#:
#: ⛔ ``v_identity`` / ``v_mutations`` / ``j_mutations`` mean the SAME thing in every mode —
#: framework-scoped. ``both`` ADDS the old numbers under new names rather than swapping the
#: meaning of a column based on a flag, which would be unreadable downstream.
FULL_COLUMNS = ("v_identity_full", "v_mutations_full", "j_mutations_full")

_SHM_FIELDS = ("v_identity", "v_mutations", "j_mutations")


def _position(entry: str) -> int | None:
    """1-based germline position of a ``G45A`` mutation entry, or None if it does not parse."""
    core = entry[1:-1]
    return int(core) if core.isdigit() else None


def _keep(muts: str, lo: int, hi: int) -> str:
    """Keep the entries of a ``G45A,C112T`` list whose position falls in ``[lo, hi]``."""
    return ",".join(e for e in muts.split(",")
                    if (p := _position(e)) is None or lo <= p <= hi)


def _anchor(rec: dict, col: str) -> int | None:
    v = rec.get(col)
    if v is None or v == "":
        return None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def scope_record(rec: dict, mode: str = "framework", identity_fn=None) -> dict:
    """Restrict ``rec``'s SHM fields to framework positions, in place. Returns ``rec``.

    ``identity_fn(qaln, taln, tstart, t_lo, t_hi)`` recomputes ``v_identity`` over a target range;
    pass :func:`arda.annotate.transfer._aln_identity`. Without it the identity is left alone (the
    mutation lists are still scoped) — a record with no alignment strings cannot be re-measured.
    """
    if mode not in SHM_MODES:
        raise ValueError(f"shm mode must be one of {SHM_MODES}, got {mode!r}")
    if mode == "off":
        for f in _SHM_FIELDS:
            if f in rec:
                rec[f] = ""
        return rec

    if mode == "both":
        for f in _SHM_FIELDS:
            rec[f"{f}_full"] = rec.get(f, "")

    v_anchor = _anchor(rec, "v_anchor_nt")
    j_anchor = _anchor(rec, "j_anchor_nt")

    if v_anchor is not None:
        if rec.get("v_mutations"):
            rec["v_mutations"] = _keep(rec["v_mutations"], 1, v_anchor)
        # The scaffold's V part IS the V germline verbatim (target pos == V germline pos), so the
        # framework cut is the anchor offset in target coordinates with no lookup.
        if identity_fn is not None and rec.get("v_identity") != "" and rec.get("sequence_alignment"):
            ts = _anchor(rec, "mmseqs2_tstart")
            t_vend = _anchor(rec, "mmseqs2_t_vend")
            if ts is not None and t_vend and ts <= min(t_vend, v_anchor):
                rec["v_identity"] = identity_fn(
                    rec["sequence_alignment"], rec["germline_alignment"],
                    ts, ts, min(t_vend, v_anchor))
    if j_anchor is not None and rec.get("j_mutations"):
        # A J germline is short and its 5' head is junction; anything past [FW]118 is FR4.
        rec["j_mutations"] = _keep(rec["j_mutations"], j_anchor + 4, 1 << 30)
    return rec


def _count_mutations(rows: list[dict]) -> int:
    """Total ``v_mutations`` + ``j_mutations`` entries across ``rows`` — the before/after number."""
    return sum(len(v.split(","))
               for r in rows
               for v in (r.get("v_mutations"), r.get("j_mutations"))
               if v)


def recount_airr(input: str | Path, output: str | Path, mode: str = "framework") -> dict:
    """Re-scope the SHM fields of an existing AIRR TSV. Returns a small report dict.

    Needs no reference and no re-map: ``v_anchor_nt`` / ``j_anchor_nt`` are already in the file,
    and so are the alignment strings ``v_identity`` is measured over. That is what makes this
    usable on tables arda wrote before 2.16.0 — the columns it needs have shipped since 2.14.0.

    ⚠ A file written by an arda older than that has no anchor columns, so nothing can be scoped.
    This RAISES rather than copying the input through unchanged: silently returning a
    byte-identical file with a success message is the exact failure `resolve_airr` shipped.
    """
    from .annotate.airr_out import read_airr
    from .annotate.transfer import _aln_identity

    df = read_airr(input)
    missing = [c for c in ("v_anchor_nt", "j_anchor_nt") if c not in df.columns]
    if missing:
        raise ValueError(
            f"{input} has no {'/'.join(missing)} column, so the junction boundary is unknown — "
            "it predates arda 2.14.0. Re-run `arda map` instead of recounting.")

    rows = df.to_dicts()
    before = _count_mutations(rows)
    for r in rows:
        scope_record(r, mode, identity_fn=_aln_identity)
    after = _count_mutations(rows)

    import polars as pl

    cols = list(df.columns) + [c for c in FULL_COLUMNS
                               if mode == "both" and c not in df.columns]
    pl.DataFrame(rows).select(cols).write_csv(output, separator="\t", quote_style="never")
    return {"rows": len(rows), "mutations_in": before, "mutations_out": after,
            "removed": before - after, "mode": mode}
