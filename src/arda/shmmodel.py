"""A sequence-context model of somatic hypermutation, fitted from arda's own output.

``arda.shm`` reports WHERE a read differs from its germline. This module turns a pile of those
reports into a model of **where a difference is expected**: ``P(substitution | germline context,
region)``, which is what a probabilistic reading of a mutated IG sequence needs and what the
V -> N1 -> D -> N2 -> J model in :mod:`arda.scenarios` is missing before it can be pointed at IGH.

Why context and a region multiplier, and NOT a per-allele-per-position table
---------------------------------------------------------------------------

``ROADMAP.md`` item 10 asked for per-allele-per-position. Measured on four bulk IGH libraries
from two donors (benchmark ``results/round33``), that object does not survive the donor change:

======================================  ==================  ==================
quantity                                within one donor    between donors
======================================  ==================  ==================
per-allele per-position profile         r .8578 - .9897     r .2557 - .5601
pooled positional profile               r .9708 / .9599     r .5994 - .6224
**5-mer context**                       r .9718 / .9785     **r .7424 - .7854**
======================================  ==================  ==================

A per-position table repeats at r ~ .99 inside one donor because the same expanded clones carry
the same mutations at the same positions -- it is a portrait of that donor's clonal history, not
a property of the allele.

⚠ **And context is the only one of the two that reaches the positions the consumer needs.** The
junction model wants ``P(observed nt | V germline)`` for the V 3' tail *inside* the junction, and
those are exactly the positions :mod:`arda.shm` scopes out as unidentifiable (a difference there
is hypermutation, chew-back or N-addition, and sequence cannot say which). A position has no rate
there; a 5-mer does, estimated from every allele that carries it in clean framework sequence.

The region multiplier is earned separately: fitting contexts on one donor and predicting the
other's per-region counts leaves a residual that is itself reproducible -- **FWR1 0.630 / 0.626
and CDR2 1.389 / 1.291** across the two directions. Context carries about half of the CDR:FWR
contrast (4.65x raw, 2.2x residual) and the region term carries the rest.

Never: the fitted **scale** is a per-sample number and is written for provenance, never applied.
The two donors differ 1.6x in overall substitution rate (.031 against .050) with the same shape,
so a shipped rate would be wrong for everyone; a consumer rescales to the sample in front of it.

Never: nothing in the annotation path calls this, exactly as ``arda scenarios`` shipped before
``arda markup --d-prior`` existed to read it. Fitting a model and adopting one are separate
decisions, and this module only does the first.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from .annotate.reference import load_reference

__all__ = [
    "V_REGIONS", "DEFAULT_K", "ShmModel", "germline_map", "read_airr", "estimate",
    "write_table", "load_model",
]

#: The V regions a read can be scored in. ``cdr3`` and ``fwr4`` are junction and J, not V.
V_REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

#: 5-mer, i.e. two germline bases either side. Odd, so the scored position is the centre.
DEFAULT_K = 5

#: A context seen fewer times than this falls back to the sample's overall rate rather than to a
#: ratio over single digits. Contexts are 4**k cells and the tail of that distribution is long.
MIN_CONTEXT = 50

_COLUMNS = ("locus", "v_call", "v_germline_start", "v_germline_end", "v_mutations", "v_anchor_nt")


@dataclass(frozen=True, slots=True)
class ShmModel:
    """``P(substitution)`` as a function of germline context and region."""

    k: int
    #: Overall substitution rate of the sample this was fitted on. Provenance, not a parameter.
    scale: float
    context: dict[str, float]
    region: dict[str, float]
    observations: int = 0
    covered: int = 0
    substitutions: int = 0
    organism: str = "human"
    locus: str = ""
    min_context: int = MIN_CONTEXT

    def rate(self, germline: str, pos: int, region: str = "") -> float:
        """``P(substitution)`` at 1-based ``pos`` of ``germline``.

        ``region`` is one of :data:`V_REGIONS`; anything else (including the junction-internal
        tail, which belongs to no V region) takes the context term alone.
        """
        half = self.k // 2
        i = pos - 1
        if i < half or i + half >= len(germline):
            p = self.scale
        else:
            p = self.context.get(germline[i - half:i + half + 1], self.scale)
        return p * self.region.get(region, 1.0)


@dataclass
class _Counts:
    """Covered and substituted tallies, keyed however the caller wants."""

    covered: dict[str, float] = field(default_factory=dict)
    mutated: dict[str, float] = field(default_factory=dict)

    def add(self, key: str, hit: bool) -> None:
        self.covered[key] = self.covered.get(key, 0.0) + 1
        if hit:
            self.mutated[key] = self.mutated.get(key, 0.0) + 1

    def rate(self, key: str) -> float:
        c = self.covered.get(key, 0.0)
        return self.mutated.get(key, 0.0) / c if c else 0.0


def germline_map(organism: str = "human", locus: str = "") -> dict[str, tuple[str, list]]:
    """``{v_call: (V germline nt, [(region, start, end)])}`` in V-germline coordinates.

    A scaffold is ``V + N-pad + J``, so a V region's scaffold coordinates ARE its germline
    coordinates and ``v_sequence_end`` is where the V germline stops. One entry per ``v_call``:
    every scaffold of a V carries the same V markup, so the first one answers for all of them.
    """
    from .refexport import _read_fasta          # the same dict[id, seq] reader export-ref uses

    ref = load_reference(organism)
    targets = _read_fasta(ref.target_fasta)
    out: dict[str, tuple[str, list]] = {}
    for sid, entry in ref.entries.items():
        if "|" in sid or not entry.v_call or entry.v_call in out:
            continue
        if locus and entry.locus != locus:
            continue
        seq = targets.get(sid, "")
        if not seq or entry.v_sequence_end <= 0:
            continue
        spans = [(name, entry.starts[i], entry.ends[i])
                 for i, name in enumerate(V_REGIONS)
                 if entry.starts[i] > 0 and 0 < entry.ends[i] <= entry.v_sequence_end]
        out[entry.v_call] = (seq[:entry.v_sequence_end], spans)
    return out


def _positions(entry: str) -> int | None:
    """1-based germline position of a ``G45A`` mutation entry, or None."""
    core = entry[1:-1]
    return int(core) if core.isdigit() else None


def read_airr(path: str | Path, locus: str = "IGH", weight: str = "unique"
              ) -> list[tuple[str, int, int, str]]:
    """``(v_call, start, capped end, mutations)`` per observation.

    ⚠ **The denominator is capped at ``v_anchor_nt``.** ``v_mutations`` is framework-scoped by
    :mod:`arda.shm` -- junction-internal entries are dropped -- so counting every covered position
    would divide a scoped numerator by an unscoped denominator and report a rate too low by
    however far the read runs past Cys104.

    ⚠ **A tie list is dropped, not resolved.** ``IGHV3-23*01,IGHV3-23D*01`` does not say which
    germline the read came from, so its differences cannot be attributed to one of them.

    ``weight="unique"`` counts a distinct ``(call, span, mutations)`` tuple once. A bulk read
    carries no junction, so there is no clonotype key to collapse on, and one expanded clone
    would otherwise vote once per read. ``weight="reads"`` is the other reading.
    """
    df = pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)
    missing = [c for c in _COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: no {', '.join(missing)} column. Annotate with --shm framework (the "
            "default) so the mutation entries and the Cys104 anchor are both present.")
    df = df.filter(
        (pl.col("locus") == locus)
        & (pl.col("v_call") != "")
        & (~pl.col("v_call").str.contains(","))
        & (pl.col("v_germline_start") != "")
        & (pl.col("v_germline_end") != "")
        & (pl.col("v_anchor_nt") != "")
    )
    rows: list[tuple[str, int, int, str]] = []
    for call, start, end, muts, anchor in df.select(
            "v_call", "v_germline_start", "v_germline_end", "v_mutations", "v_anchor_nt"
    ).iter_rows():
        cap = min(int(end), int(anchor))
        if cap >= int(start):
            rows.append((call, int(start), cap, muts or ""))
    return sorted(set(rows)) if weight == "unique" else rows


def estimate(rows, *, organism: str = "human", locus: str = "IGH", k: int = DEFAULT_K,
             min_context: int = MIN_CONTEXT) -> ShmModel:
    """Fit context rates and region multipliers from :func:`read_airr` rows.

    The region term is ``observed / expected under context alone``, normalised so the
    coverage-weighted mean multiplier is 1 -- the overall level belongs to ``scale``, which is
    the sample's own and is not a parameter of the model.

    A region carrying no observed substitution, or less than ``min_context`` covered nt, is left
    OUT of the table rather than given a multiplier of 0.
    """
    if k % 2 == 0:
        raise ValueError(f"k must be odd so the scored position is the centre of the context: {k}")
    germlines = germline_map(organism, locus)
    half = k // 2
    ctx, reg = _Counts(), _Counts()
    cells: list[tuple[str, str, bool]] = []

    for call, start, cap, muts in rows:
        found = germlines.get(call)
        if found is None:
            continue
        seq, spans = found
        mutated = {p for p in map(_positions, filter(None, muts.split(","))) if p is not None}
        for pos in range(max(start, half + 1), min(cap, len(seq) - half) + 1):
            context = seq[pos - 1 - half:pos + half]
            hit = pos in mutated
            ctx.add(context, hit)
            for name, lo, hi in spans:
                if lo <= pos <= hi:
                    reg.add(name, hit)
                    cells.append((name, context, hit))
                    break

    covered = sum(ctx.covered.values())
    substitutions = sum(ctx.mutated.values())
    scale = substitutions / covered if covered else 0.0
    rates = {c: ctx.rate(c) for c, n in sorted(ctx.covered.items()) if n >= min_context}

    expected: dict[str, float] = {}
    for name, context, _hit in cells:
        expected[name] = expected.get(name, 0.0) + rates.get(context, scale)
    # Never: a region with NO observed substitution is UNESTIMATED, never a multiplier of 0.
    # Zero would not say "thin evidence here", it would say "a substitution in FWR1 is
    # impossible" -- and `rate()` would return 0 for every position in it, which is the silent,
    # plausible-looking failure this project keeps finding. An omitted region takes `rate()`'s
    # default of 1.0, i.e. the context term alone.
    multiplier = {
        name: reg.mutated[name] / exp
        for name, exp in expected.items()
        if exp > 0 and reg.mutated.get(name, 0.0) > 0 and reg.covered.get(name, 0.0) >= min_context
    }
    weighted = sum(multiplier[n] * reg.covered[n] for n in multiplier)
    total = sum(reg.covered[n] for n in multiplier)
    if weighted and total:
        norm = weighted / total
        multiplier = {n: v / norm for n, v in multiplier.items()}

    return ShmModel(
        k=k, scale=scale, context=rates, region=dict(sorted(multiplier.items())),
        observations=len(rows), covered=int(covered), substitutions=int(substitutions),
        organism=organism, locus=locus, min_context=min_context,
    )


def write_table(model: ShmModel, path: str | Path) -> None:
    """Write the model as a TSV, sorted, with one ``#`` provenance line above the header."""
    lines = [
        f"# arda shm-model: organism={model.organism} locus={model.locus} k={model.k} "
        f"observations={model.observations} covered={model.covered} "
        f"substitutions={model.substitutions} scale={model.scale:.6g} "
        f"min_context={model.min_context}",
        "kind\tkey\tvalue",
        f"scale\tall\t{model.scale:.6g}",
    ]
    lines += [f"context\t{c}\t{v:.6g}" for c, v in sorted(model.context.items())]
    lines += [f"region\t{r}\t{v:.6g}" for r, v in sorted(model.region.items())]
    Path(path).write_text("\n".join(lines) + "\n")


def load_model(path: str | Path, *, organism: str = "human") -> ShmModel:
    """Read a table :func:`write_table` wrote.

    Never: comments, blank lines and the header are skipped by **what they are**, never by where
    they sit. ``d_prior.tsv`` shipped a loader that dropped line 1 by position, and the day a
    generator put a provenance line above the header, ``float("value")`` raised on the one file
    the docs called a drop-in.
    """
    context: dict[str, float] = {}
    region: dict[str, float] = {}
    scale, k, min_context = 0.0, DEFAULT_K, MIN_CONTEXT
    locus = ""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                for field_ in line[1:].split():
                    name, _, value = field_.partition("=")
                    if name == "k" and value.isdigit():
                        k = int(value)
                    elif name == "min_context" and value.isdigit():
                        min_context = int(value)
                    elif name == "locus":
                        locus = value
                continue
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 3 or parts[0] == "kind":
                continue
            kind, key, raw = parts
            try:
                value = float(raw)
            except ValueError:
                continue
            if kind == "context":
                context[key] = value
            elif kind == "region":
                region[key] = value
            elif kind == "scale":
                scale = value
    return ShmModel(k=k, scale=scale, context=context, region=region, organism=organism,
                    locus=locus, min_context=min_context)
