"""Tie lists — the germlines a read's alignment cannot distinguish from the one it was called on.

Why this exists
---------------
A read aligned to germline ``G`` over ``[gstart, gend)`` is explained *exactly as well* by any
other germline containing that same stretch: the read carries no base that separates them, so
naming only ``G`` is a claim the data does not support. IgBLAST says so — 11.68 % of its V calls on
a Ramos library are multi-gene tie lists — and arda, on the same library, emitted **0 in 504**.

That silence is not an accuracy win, it is a coin flip. Measured on the 104 reads where arda said
``IGLV2-23`` and IgBLAST said ``IGLV2-14``: aligning each read to **both** germlines, **59 of 60**
fit identically, typically at identity **1.0000 over 63–70 nt**. One read favoured IgBLAST's call,
none favoured arda's. Neither tool was right; both were overconfident, and they disagreed by
accident.

⛔ Why arda lost the ties in the first place, and why this does not undo that
----------------------------------------------------------------------------
``top_hit`` runs before ``convertalis`` deliberately: aligning every exactly-tied allele made the
alignment TSV **2.88× larger** — the same shape as the 194 MB → 877 MB regression ``top_hit``
exists to prevent — *and* collapsing to one hit fixed allele-level agreement, **.9735 → .9956**, by
making both code paths break ties by one rule instead of two. Restoring tie lists by aligning more
candidates would give both problems back.

This does not align anything. Given the span the read **already** aligned over, a tie is a
**string comparison against the reference**: any germline containing ``G[gstart:gend]`` verbatim.
So it costs neither the memory nor the alignment time ``top_hit`` was protecting.

Cost control
------------
The search is in C++ (:func:`arda._denoise.containing`) and the result is memoised on
``(allele, gstart, gend)``. On a primer-anchored amplicon that collapses hundreds of thousands of
reads onto a few thousand distinct spans, so the search runs once per span, not once per read.

⛔ Off by default. This changes ``v_call``/``j_call`` on every library, and a downstream consumer
that splits on ``,`` and takes ``[0]`` sees no change while one that treats the field as a single
gene sees a new shape. Turn it on with ``--tie-lists``.
"""

from __future__ import annotations

from functools import lru_cache

__all__ = ["TieResolver", "rank_ties", "resolve_airr"]

try:                                              # optional: built by scikit-build-core
    from .. import _denoise as _cpp
except ImportError:                               # pragma: no cover - source checkout without ext
    _cpp = None


def _containing_py(segment: str, candidates: list[str]) -> list[int]:
    """Reference implementation of :func:`arda._denoise.containing`."""
    if not segment:
        return []
    return [i for i, c in enumerate(candidates) if segment in c]


class TieResolver:
    """Expand a single call into every germline the aligned span cannot rule out.

    ``germlines`` maps allele name -> ungapped nucleotide sequence, i.e. exactly what
    ``refbuild.imgt.load_functional_alleles`` returns, restricted to one segment type.
    """

    #: A span shorter than this is not evidence of anything — nearly every allele of a family
    #: contains a 20 nt stretch of its neighbours, so the "tie list" would be the whole family and
    #: would say less than the single call it replaced. Reads this short do not get a tie list.
    MIN_SPAN = 30

    def __init__(self, germlines: dict[str, str], *, max_ties: int = 16):
        # Sorted once, so the emitted call string does not depend on dict order — `correct`'s
        # clonotype key is built from these strings and this project has already shipped a
        # nondeterministic Stage 2 for exactly this class of reason.
        self._names = sorted(germlines)
        self._seqs = [germlines[n] for n in self._names]
        self._index = {n: i for i, n in enumerate(self._names)}
        self._max_ties = max_ties
        self._expand = lru_cache(maxsize=200_000)(self._expand_uncached)

    def _expand_uncached(self, allele: str, gstart: int, gend: int) -> tuple[str, ...]:
        i = self._index.get(allele)
        if i is None or gend - gstart < self.MIN_SPAN:
            return ()
        segment = self._seqs[i][gstart:gend]
        if len(segment) < self.MIN_SPAN:
            return ()
        hits = (_cpp.containing(segment, self._seqs) if _cpp is not None
                else _containing_py(segment, self._seqs))
        names = tuple(self._names[j] for j in hits)
        # ⛔ A runaway tie list is worse than no tie list: it turns one wrong-but-usable call into
        # an unusable one, and it inflates every downstream string. Above the cap the call is left
        # exactly as it was.
        return () if len(names) > self._max_ties else names

    def expand(self, call: str, gstart, gend) -> str:
        """``call`` widened to every indistinguishable germline, or ``call`` unchanged.

        ``call`` may already be a comma-joined list; the first element is the one that aligned, and
        anything the caller already put there is preserved and comes first. Coordinates are the
        AIRR 1-based closed ``*_germline_start`` / ``*_germline_end``; empty or malformed
        coordinates return the call untouched, because a tie list computed over a span that is not
        known is a guess.
        """
        if not call:
            return call
        try:
            s, e = int(gstart), int(gend)
        except (TypeError, ValueError):
            return call
        if s <= 0 or e < s:
            return call
        first = call.split(",")[0].strip()
        names = self._expand(first, s - 1, e)          # AIRR is 1-based closed
        if len(names) <= 1:
            return call
        existing = [x.strip() for x in call.split(",") if x.strip()]
        merged = existing + [n for n in names if n not in existing]
        return ",".join(merged)

    def cache_info(self):
        return self._expand.cache_info()


def rank_ties(calls: list[str], scores: list[float] | None = None,
              evidence: list[str] | None = None) -> list[str]:
    """Second pass: reorder each tie list so the allele the WHOLE LIBRARY supports leads.

    A tie list says the read cannot choose. The library usually can: the same gene is seen on other
    reads that *do* reach a discriminating base, and those reads are evidence this read does not
    have. So rank the members by that evidence and put the winner first, leaving the rest as the
    honest statement of what this read alone could not rule out.

    ⛔ **Only UNAMBIGUOUS reads vote.** A read whose own call is ``A,B`` cannot be evidence for A
    over B — counting it would let a common tie bootstrap itself, and the more often two alleles
    are confusable the more confidently the pair would elect one of them. Ties are ranked, never
    ranking. If no member was ever seen unambiguously, the summed score over all reads naming it
    breaks the tie, and if that is level too the order is lexicographic — so the output never
    depends on row order.

    ``scores`` is an optional per-row weight (mmseqs2 bit score is the natural one); ``None``
    counts reads. Returns a new list; the input is not modified. The leading element is all that
    moves — membership is unchanged, so a consumer taking ``[0]`` gets the better answer and one
    reading the whole field still sees everything that was indistinguishable.
    """
    if scores is not None and len(scores) != len(calls):
        raise ValueError("scores must be the same length as calls")
    # ⛔ The votes come from `evidence`, which must be the calls BEFORE tie expansion. Ranking on
    # the expanded calls is self-defeating: expansion makes every read ambiguous, so the
    # unambiguous-reads-only rule has nothing left to count and the whole thing degenerates to
    # lexicographic order. Caught on the real IGLV2-14/IGLV2-23 pair, where the expanded ranking
    # elected `IGLV2-14*01` purely because it sorts first.
    src = calls if evidence is None else evidence
    if evidence is not None and len(evidence) != len(calls):
        raise ValueError("evidence must be the same length as calls")
    solo: dict[str, float] = {}
    total: dict[str, float] = {}
    for i, c in enumerate(src):
        if not c:
            continue
        parts = [x.strip() for x in c.split(",") if x.strip()]
        w = 1.0 if scores is None else float(scores[i] or 0.0)
        for p in parts:
            total[p] = total.get(p, 0.0) + w
        if len(parts) == 1:
            solo[parts[0]] = solo.get(parts[0], 0.0) + w

    def key(name: str):
        # Unambiguous evidence first, then total, then the name -- a TOTAL order, so two runs over
        # the same data cannot disagree.
        return (-solo.get(name, 0.0), -total.get(name, 0.0), name)

    out = []
    for c in calls:
        if not c or "," not in c:
            out.append(c)
            continue
        parts = [x.strip() for x in c.split(",") if x.strip()]
        best = min(parts, key=key)
        if best == parts[0]:
            out.append(c)
        else:
            out.append(",".join([best] + [p for p in parts if p != best]))
    return out


def resolve_airr(path, out, *, organism: str = "human", segments: tuple[str, ...] = ("v", "j"),
                 rank: bool = True, echo=None) -> dict:
    """Add tie lists to an AIRR TSV's ``v_call``/``j_call``, then rank them library-wide.

    Two passes over one file, which is why this is a separate step rather than something ``map``
    does inline: the ranking needs every read before it can order any of them.

    ⛔ Membership is decided per read, from the span that read aligned; only the ORDER is decided
    library-wide. Nothing is added or removed by the second pass, so a consumer taking the first
    element gets a better answer and one reading the whole field still sees every germline the read
    could not rule out.
    """
    import polars as pl

    from ..refbuild import imgt
    from ..refbuild.loci import IMGT_SPECIES_DIR, loci_for

    species_dir = IMGT_SPECIES_DIR[organism]
    df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
    report = {"rows": df.height, "expanded": {}, "reranked": {}}

    for seg in segments:
        call_col, gs, ge = f"{seg}_call", f"{seg}_germline_start", f"{seg}_germline_end"
        if call_col not in df.columns or gs not in df.columns:
            continue
        germ: dict[str, str] = {}
        for locus in loci_for():
            stem = getattr(locus, seg, None)
            if not stem:
                continue
            try:
                germ.update(imgt.load_functional_alleles(species_dir, locus.group, stem))
            except Exception:                      # a locus whose germline files are absent
                continue
        if not germ:
            continue
        res = TieResolver(germ)
        calls = df[call_col].to_list()
        starts, ends = df[gs].to_list(), df[ge].to_list()
        widened = [res.expand(c or "", a, b) for c, a, b in zip(calls, starts, ends)]
        report["expanded"][seg] = sum(1 for a, b in zip(calls, widened) if (a or "") != b)
        if rank:
            scores = (df["mmseqs2_score"].to_list() if "mmseqs2_score" in df.columns else None)
            sc = None
            if scores is not None:
                sc = []
                for x in scores:
                    try:
                        sc.append(float(x))
                    except (TypeError, ValueError):
                        sc.append(0.0)
            ranked = rank_ties(widened, sc, evidence=[c or "" for c in calls])
            report["reranked"][seg] = sum(1 for a, b in zip(widened, ranked) if a != b)
            widened = ranked
        df = df.with_columns(pl.Series(call_col, widened))
        if echo:
            echo(f"[arda] {seg}_call: +{report['expanded'].get(seg, 0)} tie lists, "
                 f"{report['reranked'].get(seg, 0)} reordered")
    df.write_csv(out, separator="\t", quote_style="never")
    return report
