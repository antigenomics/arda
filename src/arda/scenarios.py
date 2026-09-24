"""Recombination scenarios from a nucleotide junction, and the counts they imply.

`arda.dpost` places a D from an amino-acid junction by marginalising a generative model, and
every number in that model is borrowed from OLGA/vdjrearm and shipped as
``database/vdj/<org>/d_prior.tsv``. This module is how arda estimates its own: read
``(junction_nt, v_call, j_call)`` records, count recombination scenarios, and write a table of
the same ``locus / kind / key / value`` shape, so it is a drop-in for the shipped file.

**Junction space**, as in :mod:`arda.cdr3fix` and :mod:`arda.annotate.dmap`: Cys104 -> [FW]118
inclusive. ``Vg`` and ``Jg`` are the **untrimmed** junction-region germlines already shipped per
allele in ``cdr3_anchors.tsv``; nothing here reads a FASTA.

A scenario reproduces the junction exactly::

    junction == Vg[:len(Vg)-del_v] + N1 + Dg[del_dl:len(Dg)-del_dr] + N2 + Jg[del_j:]
                                     ins_vd                           ins_dj

Never: **a scenario is not identifiable from sequence.** Chew-back and N/P addition make several
tuples explain one junction exactly -- a first non-templated base that happens to match germline
is indistinguishable from one less nucleotide of trimming. Counting one reading (which is what
``dmap.map_d_junction`` gives, and it is the right thing for *annotation*) assigns weight 1 to one
member of that set and 0 to the rest, biasing every distribution toward less trimming and shorter
inserts. So the tallies here are **expected counts under the current model**, summed over
scenarios -- which, with renormalisation, is EM.

See ``project/design-scenarios.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .cdr3fix import load_anchors, resolve_allele, resolve_locus, resolve_species
from .annotate.reference import _load_d_germlines
from .paths import vdj_dir

__all__ = ["Scenario", "Germlines", "germlines_for", "enumerate_scenarios",
           "SufficientStats", "estimate", "accumulate", "lattice", "PRIOR_COLUMNS"]

#: The long table both this module and ``dpost`` speak. Same idiom as ``stats.py``.
PRIOR_COLUMNS = ("locus", "kind", "key", "value")

#: Kinds ``dpost.load_d_prior`` actually consumes, plus the two trimming distributions nothing
#: consumed before because nothing produced them. ``beta`` is a fitted temperature, not a count,
#: and is deliberately absent -- see ``project/design-scenarios.md``.
KINDS = ("insVD", "insDJ", "dlen", "d_marginal", "d_given_j", "delV", "delJ")


@dataclass(frozen=True, slots=True)
class Scenario:
    """One exact reading of a junction. ``d_call`` is ``""`` at a VJ locus."""

    del_v: int = 0
    ins_vd: int = 0
    d_call: str = ""
    del_dl: int = 0
    del_dr: int = 0
    ins_dj: int = 0
    del_j: int = 0


@dataclass(frozen=True, slots=True)
class Germlines:
    """The three germline strings one ``(V, J)`` pair needs, in junction space."""

    locus: str
    v_allele: str
    j_allele: str
    v_nt: str
    j_nt: str
    d_germlines: tuple[tuple[str, str], ...] = ()


@lru_cache(maxsize=8)
def _d_germlines(organism: str) -> dict[str, list[tuple[str, str]]]:
    # Cached for the reason `dmap._d_germlines` is: `_load_d_germlines` re-parses the FASTA on
    # every call, and this runs once per record.
    return _load_d_germlines(vdj_dir(organism))


@lru_cache(maxsize=65536)
def germlines_for(v_call: str, j_call: str, species: str = "human") -> Germlines | None:
    """Resolve one ``(V, J)`` call pair to its junction-space germlines, or ``None``.

    Cached on the call strings: a repertoire has tens of thousands of clonotypes over a few
    hundred distinct pairs, so this resolves each pair once.
    """
    organism = resolve_species(species)
    locus = resolve_locus(v_call, j_call)
    if not locus:
        return None
    anchors = load_anchors(organism)
    v_id = resolve_allele((v_call or "").split(",")[0], "V", anchors)
    j_id = resolve_allele((j_call or "").split(",")[0], "J", anchors)
    va, ja = anchors.get(("V", v_id)), anchors.get(("J", j_id))
    if not va or not ja or va.status != "ok" or ja.status != "ok":
        return None
    if not va.germline_nt or not ja.germline_nt:
        return None
    return Germlines(
        locus=locus, v_allele=v_id, j_allele=j_id,
        v_nt=va.germline_nt.upper(), j_nt=ja.germline_nt.upper(),
        d_germlines=tuple((n, s.upper()) for n, s in (_d_germlines(organism).get(locus) or [])),
    )


def _common_prefix(a: str, b: str) -> int:
    n = 0
    while n < len(a) and n < len(b) and a[n] == b[n]:
        n += 1
    return n


def _common_suffix(a: str, b: str) -> int:
    n = 0
    while n < len(a) and n < len(b) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


def _d_placements(junction: str, d_germlines) -> list[tuple[str, int, int, int, int]]:
    """Every ``(allele, start, length, del_dl, del_dr)`` a D germline can occupy exactly.

    Never: ``length == 0`` is emitted, once per allele, and it is not a failure -- it is the modal
    outcome. The shipped table's largest single entry is ``IGHD1-1*01:0 = 0.8759``. A VDJ locus
    always rearranged a D; dropping the records where none survived would truncate ``dlen`` at 1
    and inflate every insertion distribution by the length of the D that was really there.
    Its ``del_dl``/``del_dr`` split is unidentifiable, so a zero-length placement reports
    ``(0, len(germline))`` and contributes to ``dlen`` only.
    """
    seen: dict[tuple[str, int, int], tuple[int, int]] = {}
    n = len(junction)
    for allele, germ in d_germlines:
        g = len(germ)
        # Offsets indexed by first base: three quarters of the (position, offset) pairs cannot
        # match at all, and testing them was the single largest cost in the estimator once the
        # per-span fix below removed the accumulation overhead.
        by_base: dict[str, list[int]] = {}
        for off, ch in enumerate(germ):
            by_base.setdefault(ch, []).append(off)
        for p in range(n):
            for off in by_base.get(junction[p], ()):
                m = 1
                while p + m < n and off + m < g and junction[p + m] == germ[off + m]:
                    m += 1
                # Never: keyed on (allele, start, LENGTH) -- two germline offsets giving the same
                # surviving string at the same place are ONE state of this model, whose D term is
                # `dlen` (surviving length) and carries no `del_dl`/`del_dr` joint. Counting them
                # separately would weight a D with internal repeats by its repeat count.
                if m and (allele, p, m) not in seen:
                    seen[(allele, p, m)] = (off, g - off - m)
    return [(a, p, m, dl, dr) for (a, p, m), (dl, dr) in seen.items()]


def _dist(sparse: dict[int, float], n: int) -> list[float]:
    """A length-``n`` distribution from a sparse map, uniform when the map is empty."""
    if not sparse:
        return [1.0 / n] * n if n > 0 else []
    out = [sparse.get(i, 0.0) for i in range(n)]
    total = sum(out)
    # A feasible value the shipped model gives zero mass cannot be ruled out by a model that was
    # fitted on someone else's data; a flat floor keeps it reachable so EM can move it.
    if total <= 0:
        return [1.0 / n] * n if n > 0 else []
    floor = 1e-9
    out = [max(x, floor) for x in out]
    total = sum(out)
    return [x / total for x in out]


@dataclass
class SufficientStats:
    """Expected counts, keyed exactly as ``d_prior.tsv`` keys them."""

    counts: dict[tuple[str, str, str], float]     # (locus, kind, key) -> expected count
    records: int = 0
    skipped: int = 0
    log_likelihood: float = 0.0

    def add(self, locus: str, kind: str, key: str, w: float) -> None:
        if w:
            k = (locus, kind, key)
            self.counts[k] = self.counts.get(k, 0.0) + w

    def rows(self) -> list[tuple[str, str, str, float]]:
        """Normalised within each ``(locus, kind)`` group, or within its allele for keyed kinds."""
        groups: dict[tuple[str, str, str], float] = {}
        for (locus, kind, key), v in self.counts.items():
            # `dlen` normalises per D allele and `d_given_j` per J allele -- those are conditional
            # distributions, and normalising them globally would make every allele's curve depend
            # on how often that allele was used.
            sub = key.rpartition(":")[0] if kind == "dlen" else (
                key.partition("|")[2] if kind == "d_given_j" else
                key.rpartition(":")[0] if kind in ("delV", "delJ") else "")
            groups[(locus, kind, sub)] = groups.get((locus, kind, sub), 0.0) + v
        out = []
        for (locus, kind, key), v in self.counts.items():
            sub = key.rpartition(":")[0] if kind == "dlen" else (
                key.partition("|")[2] if kind == "d_given_j" else
                key.rpartition(":")[0] if kind in ("delV", "delJ") else "")
            total = groups[(locus, kind, sub)]
            if total > 0:
                out.append((locus, kind, key, v / total))
        return sorted(out)


class _Model:
    """The current parameters, as flat lookup tables. Initialised from the shipped prior."""

    def __init__(self, organism: str, max_del: int = 24, max_ins: int = 30):
        from .dpost import load_d_prior

        self.max_del, self.max_ins = max_del, max_ins
        self.prior = load_d_prior(organism)
        self.ins_vd: dict[str, list[float]] = {}
        self.ins_dj: dict[str, list[float]] = {}
        self.dlen: dict[str, dict[str, list[float]]] = {}
        self.d_given_j: dict[str, dict[str, dict[str, float]]] = {}
        self.del_v: dict[str, dict[str, list[float]]] = {}
        self.del_j: dict[str, dict[str, list[float]]] = {}

    def _ins(self, table: dict, locus: str, shipped: str) -> list[float]:
        if locus not in table:
            p = self.prior.get(locus)
            src = getattr(p, shipped, None) if p else None
            sparse = {i: v for i, v in enumerate(src or []) if v > 0}
            table[locus] = _dist(sparse, self.max_ins + 1)
        return table[locus]

    def ins_vd_of(self, locus: str) -> list[float]:
        return self._ins(self.ins_vd, locus, "ins_vd")

    def ins_dj_of(self, locus: str) -> list[float]:
        return self._ins(self.ins_dj, locus, "ins_dj")

    def dlen_of(self, locus: str, allele: str, g: int) -> list[float]:
        per = self.dlen.setdefault(locus, {})
        if allele not in per:
            p = self.prior.get(locus)
            src = (p.dlen.get(allele) if p else None) or []
            per[allele] = _dist({i: v for i, v in enumerate(src) if v > 0}, g + 1)
        return per[allele]

    def d_given_j_of(self, locus: str, j_allele: str, alleles) -> dict[str, float]:
        per = self.d_given_j.setdefault(locus, {})
        if j_allele not in per:
            p = self.prior.get(locus)
            table = (p.d_given_j.get(j_allele) or p.d_marginal) if p else {}
            vals = {a: max(float(table.get(a, 0.0)), 1e-9) for a in alleles}
            total = sum(vals.values()) or 1.0
            per[j_allele] = {a: v / total for a, v in vals.items()}
        return per[j_allele]

    def _del(self, table: dict, locus: str, allele: str, germline_len: int) -> list[float]:
        """P(deletion) for one allele, **indexed by DELETION COUNT**, length ``germline_len+1``.

        Never: indexed by deletion, not by surviving length. They are mirror images only when the
        germline is exactly as long as the cap, so indexing one table as the other silently reads
        `P(del=k)` as `P(del=len-k)` on every allele of a different length -- the distributions
        come out plausible, reversed, and per-allele wrong.

        Never: uniform, not geometric. Nothing shipped a `delV`/`delJ` distribution because
        nothing consumed one, so there is no prior to start from -- and inventing a parametric
        shape would put a guess into the estimate that EM only partly removes. Uniform over the
        feasible range is the one start that asserts nothing.
        """
        per = table.setdefault(locus, {})
        if allele not in per:
            per[allele] = [1.0 / (germline_len + 1)] * (germline_len + 1)
        return per[allele]

    def del_v_of(self, locus: str, allele: str, germline_len: int) -> list[float]:
        return self._del(self.del_v, locus, allele, germline_len)

    def del_j_of(self, locus: str, allele: str, germline_len: int) -> list[float]:
        return self._del(self.del_j, locus, allele, germline_len)


def _templated(p_del: list[float], germline_len: int, t_max: int, max_del: int) -> list[float]:
    """Re-index a deletion distribution by templated length, zeroing what the cap forbids."""
    out = []
    for t in range(t_max + 1):
        d = germline_len - t
        out.append(p_del[d] if 0 <= d < len(p_del) and d <= max_del else 0.0)
    return out


#: P(one non-templated base). Never: **an insertion must cost its own sequence, not just its
#: length.** A scenario explains the junction's 5' end either as templated V or as an insertion
#: that happens to match V germline, and with a length term alone the second is FREE -- so EM
#: walks straight into the degenerate corner where everything is insertion. Measured here before
#: this term existed: three iterations on 503 real TRB junctions moved `insVD` mass onto 10-11 nt
#: and the log-likelihood rose monotonically the whole way, which is EM doing exactly what it was
#: asked. P(sequence | length) = 0.25^length restores the 2-bits-per-base cost that makes an
#: accidental 10-nt germline match 10^-6 as likely as the templated reading.
#: Uniform rather than a fitted first-order Markov (IGoR's choice): composition is a refinement,
#: the cost per base is what breaks the degeneracy, and a fitted chain here would be one more
#: thing initialised from nothing.
_P_BASE = 0.25


def _side(p_len: list[float], p_ins: list[float], span: int, max_side: int
          ) -> tuple[float, list[tuple[int, int, float]]]:
    """``(total, [(germline_len, ins_len, weight)])`` for one flank covering ``span`` nt.

    The flank contributes ``germline_len`` templated nt and ``ins_len = span - germline_len``
    non-templated ones. This is the factor that makes the whole enumeration cheap: the V side and
    the J side are conditionally independent given the D placement, so each is summed once per
    span rather than once per (V, D, J) triple.
    """
    parts: list[tuple[int, int, float]] = []
    total = 0.0
    for glen in range(0, min(max_side, span) + 1):
        ins = span - glen
        if ins >= len(p_ins) or glen >= len(p_len):
            continue
        w = p_len[glen] * p_ins[ins] * (_P_BASE ** ins)
        if w > 0:
            parts.append((glen, ins, w))
            total += w
    return total, parts


def lattice(junction_nt: str, v_call: str, j_call: str, model: _Model,
            species: str = "human"):
    """``(g, L, terms, left, right)`` for one junction, or ``None``.

    **This is the forward-backward pass.** ``terms`` is every state of the semi-Markov chain
    V -> N1 -> D -> N2 -> J that emits this junction, each already carrying its unnormalised
    posterior weight; ``left``/``right`` are the flank factors, summed once per span because the
    two flanks are conditionally independent given the D placement. Their sum is the partition
    function, i.e. ``P(junction | V, J)`` under the current parameters.

    :mod:`arda.hmm` is this function read as inference and :func:`accumulate` is it read as an
    E-step. They are the same recursion, which is why there is one implementation.
    """
    g = germlines_for(v_call, j_call, species)
    junction = (junction_nt or "").strip().upper()
    if g is None or not junction:
        return None
    locus, L = g.locus, len(junction)
    v_max = _common_prefix(junction, g.v_nt)
    j_max = _common_suffix(junction, g.j_nt)
    # `del_*_of` is indexed by DELETION; `_side` wants an index of TEMPLATED LENGTH. Templated
    # `t` means `len(germline) - t` deleted, and `max_del` caps the deletion -- so short templated
    # lengths are what the cap removes, not long ones.
    pv = _templated(model.del_v_of(locus, g.v_allele, len(g.v_nt)), len(g.v_nt),
                    v_max, model.max_del)
    pj = _templated(model.del_j_of(locus, g.j_allele, len(g.j_nt)), len(g.j_nt),
                    j_max, model.max_del)
    v_cap, j_cap = len(pv) - 1, len(pj) - 1
    ins_vd, ins_dj = model.ins_vd_of(locus), model.ins_dj_of(locus)

    # (v_templated, v_ins, w) x (j_templated, j_ins, w) per span, computed once.
    left = [_side(pv, ins_vd, p, v_cap) for p in range(L + 1)]

    terms: list[tuple[float, int, int, str, int]] = []           # (w, p, q, allele, m)
    right: list = []
    if g.d_germlines:
        right = [_side(pj, ins_dj, q, j_cap) for q in range(L + 1)]
        alleles = [n for n, _ in g.d_germlines]
        pdj = model.d_given_j_of(locus, g.j_allele, alleles)
        lens = {n: len(s) for n, s in g.d_germlines}
        places = _d_placements(junction, g.d_germlines)
        # The surviving-nothing case, at every position it could have sat: `dlen == 0` is the
        # modal outcome, not a failure (see `_d_placements`).
        places += [(a, p, 0, 0, lens[a]) for a in alleles for p in range(L + 1)]
        for allele, p, m, _dl, _dr in places:
            q = L - p - m
            if q < 0:
                continue
            lt, rt = left[p][0], right[q][0]
            if lt <= 0 or rt <= 0:
                continue
            w = lt * rt * pdj.get(allele, 0.0) * model.dlen_of(locus, allele, lens[allele])[m]
            if w > 0:
                terms.append((w, p, q, allele, m))
    else:
        # VJ locus: one non-templated run, no D and no insDJ. Never skipped -- TRA/TRG/IGK/IGL
        # carry most of the trimming evidence, and `delV`/`delJ` are the point of collecting it.
        for v_t in range(0, min(v_cap, L) + 1):
            for j_t in range(0, min(j_cap, L - v_t) + 1):
                ins = L - v_t - j_t
                if ins >= len(ins_vd):
                    continue
                w = pv[v_t] * pj[j_t] * ins_vd[ins] * (_P_BASE ** ins)
                if w > 0:
                    terms.append((w, v_t, j_t, "", -1))

    return g, L, terms, left, (right if g.d_germlines else [])


def accumulate(junction_nt: str, v_call: str, j_call: str, model: _Model,
               stats: SufficientStats, weight: float = 1.0, species: str = "human") -> bool:
    """Add one record's expected counts. Returns ``False`` when it could not be read."""
    built = lattice(junction_nt, v_call, j_call, model, species)
    if built is None:
        return False
    g, L, terms, left, right = built
    locus = g.locus
    z = sum(t[0] for t in terms)
    if z <= 0:
        return False
    stats.records += 1
    stats.log_likelihood += weight * _log(z)

    if not g.d_germlines:
        for w, v_t, j_t, _a, _m in terms:
            f = weight * w / z
            stats.add(locus, "delV", f"{g.v_allele}:{len(g.v_nt) - v_t}", f)
            stats.add(locus, "delJ", f"{g.j_allele}:{len(g.j_nt) - j_t}", f)
            stats.add(locus, "insVD", str(L - v_t - j_t), f)
        return True

    # Never: the flank weight is accumulated PER SPAN, then distributed once -- never per term.
    # A flank's split distribution depends only on its span, and the contribution is linear in
    # the term's posterior, so summing the posterior into `wl[p]` first is exactly equivalent and
    # collapses `terms x splits` adds into `terms + spans x splits`. Measured on 2,496 real
    # records: 73.5 M `add` calls and 17.7 s per EM iteration before, and `add` plus its `dict.get`
    # were 45 % of the whole run -- not the germline matching anyone would profile first.
    wl = [0.0] * (L + 1)
    wr = [0.0] * (L + 1)
    for w, p, q, allele, m in terms:
        f = weight * w / z
        stats.add(locus, "dlen", f"{allele}:{m}", f)
        stats.add(locus, "d_marginal", allele, f)
        stats.add(locus, "d_given_j", f"{allele}|{g.j_allele}", f)
        wl[p] += f
        wr[q] += f
    for p, f in enumerate(wl):
        if not f:
            continue
        lt, lparts = left[p]
        for glen, ins, lw in lparts:
            share = f * lw / lt
            stats.add(locus, "delV", f"{g.v_allele}:{len(g.v_nt) - glen}", share)
            stats.add(locus, "insVD", str(ins), share)
    for q, f in enumerate(wr):
        if not f:
            continue
        rt, rparts = right[q]
        for glen, ins, rw in rparts:
            share = f * rw / rt
            stats.add(locus, "delJ", f"{g.j_allele}:{len(g.j_nt) - glen}", share)
            stats.add(locus, "insDJ", str(ins), share)
    return True


def _log(x: float) -> float:
    from math import log
    return log(x) if x > 0 else float("-inf")


def enumerate_scenarios(junction_nt: str, v_call: str, j_call: str,
                        species: str = "human", *, max_del: int = 24) -> list[Scenario]:
    """Every scenario that reproduces ``junction_nt`` exactly, unweighted.

    For inspection and for tests. The estimator does **not** call this -- it sums the same set in
    factorised form, because the V and J flanks are conditionally independent given the D
    placement and enumerating the product explicitly is what makes a naive version slow.
    """
    g = germlines_for(v_call, j_call, species)
    junction = (junction_nt or "").strip().upper()
    if g is None or not junction:
        return []
    L = len(junction)
    v_max = min(_common_prefix(junction, g.v_nt), max_del)
    j_max = min(_common_suffix(junction, g.j_nt), max_del)
    out: list[Scenario] = []
    if not g.d_germlines:
        for v_t in range(v_max + 1):
            for j_t in range(min(j_max, L - v_t) + 1):
                out.append(Scenario(del_v=len(g.v_nt) - v_t, ins_vd=L - v_t - j_t,
                                    del_j=len(g.j_nt) - j_t))
        return out
    lens = {n: len(s) for n, s in g.d_germlines}
    places = list(_d_placements(junction, g.d_germlines))
    places += [(a, p, 0, 0, lens[a]) for a in lens for p in range(L + 1)]
    for allele, p, m, dl, dr in places:
        q = L - p - m
        if q < 0:
            continue
        for v_t in range(min(v_max, p) + 1):
            for j_t in range(min(j_max, q) + 1):
                out.append(Scenario(
                    del_v=len(g.v_nt) - v_t, ins_vd=p - v_t, d_call=allele,
                    del_dl=dl, del_dr=dr, ins_dj=q - j_t, del_j=len(g.j_nt) - j_t))
    return out


def estimate(records, *, organism: str = "human", iterations: int = 5,
             max_del: int = 24, echo=None) -> SufficientStats:
    """EM over ``(junction_nt, v_call, j_call, weight)`` records; returns the final counts.

    Never: **weighted by the model, never by 1/n over the ambiguity set.** Uniform weighting is
    itself a strong and wrong prior -- it says a 12-nt insertion is as likely as a 2-nt one, and
    the ambiguity sets are much larger for long insertions, so it would pull every distribution
    toward them. Initialisation is the shipped prior where one exists (``insVD``, ``insDJ``,
    ``dlen``, ``d_given_j``); ``delV``/``delJ`` start uniform because nothing shipped them.

    Never: **each record contributes 1.0 in total, not one count per scenario.** A junction with
    400 admissible readings must not outvote one with 3.
    """
    records = list(records)
    model = _Model(organism, max_del=max_del)
    stats = SufficientStats(counts={})
    for it in range(max(1, iterations)):
        stats = SufficientStats(counts={})
        for junction, v_call, j_call, weight in records:
            if not accumulate(junction, v_call, j_call, model, stats, weight, organism):
                stats.skipped += 1
        if echo:
            echo(f"scenarios: iteration {it + 1}/{iterations} -- {stats.records} records, "
                 f"log-likelihood {stats.log_likelihood:.1f}")
        if it + 1 < iterations:
            model = _refit(model, stats, organism, max_del)
    return stats


def _refit(old: _Model, stats: SufficientStats, organism: str, max_del: int) -> _Model:
    """The M-step: the normalised expected counts become the next iteration's parameters."""
    m = _Model(organism, max_del=max_del, max_ins=old.max_ins)
    by_kind: dict[tuple[str, str], dict[str, float]] = {}
    for locus, kind, key, value in stats.rows():
        by_kind.setdefault((locus, kind), {})[key] = value
    for (locus, kind), table in by_kind.items():
        if kind in ("insVD", "insDJ"):
            dense = _dist({int(k): v for k, v in table.items()}, old.max_ins + 1)
            (m.ins_vd if kind == "insVD" else m.ins_dj)[locus] = dense
        elif kind == "dlen":
            for key, v in table.items():
                allele, _, n = key.rpartition(":")
                m.dlen.setdefault(locus, {}).setdefault(allele, {})[int(n)] = v
        elif kind == "d_given_j":
            for key, v in table.items():
                d_allele, _, j_allele = key.partition("|")
                m.d_given_j.setdefault(locus, {}).setdefault(j_allele, {})[d_allele] = v
        elif kind in ("delV", "delJ"):
            for key, v in table.items():
                allele, _, n = key.rpartition(":")
                dst = m.del_v if kind == "delV" else m.del_j
                dst.setdefault(locus, {}).setdefault(allele, {})[int(n)] = v
    # The per-allele tables were filled as sparse dicts; densify them to the shape the E-step
    # indexes, keeping the flat floor so a value this pass gave zero mass stays reachable.
    for locus, per in m.dlen.items():
        for allele, sparse in list(per.items()):
            per[allele] = _dist(sparse, max(sparse) + 1 if sparse else 1)
    for table in (m.del_v, m.del_j):
        for locus, per in table.items():
            for allele, sparse in list(per.items()):
                per[allele] = _dist(sparse, max(sparse) + 1 if sparse else 1)
    return m
