"""Which D, and where — from an amino-acid junction alone.

A VDJdb-style record has no nucleotides, and D segments are short and trimmed at both
ends, so the D is often invisible in the translated junction. Two independent sources of
information remain, and they are complementary:

**Where.** The junction's nucleotide length is known (3x its amino-acid length), so
``insVD + |D surviving| + insDJ`` is *pinned*. Marginalising the generative model's
insertion-length and D-trimming distributions therefore places the D even when the
sequence says nothing at all about it. Measured against OLGA ground truth, the MAP
``d_start`` is a median 1 nt off for mouse TRB, 2 nt for human TRB, and 3 nt for TRD and
IGH.

**Which.** The length constraint is nearly useless for identity — the D length
distributions overlap, so the posterior barely moves off the prior. What identity the
prior does carry is ``P(D | J)``, and for TRB that is mostly genomic order: TRBD2 lies 3'
of the whole TRBJ1 cluster, so a TRBJ1 junction can only have used TRBD1 (see
``_mask_forbidden``). What otherwise identifies a D is the amino-acid match, and only where
enough D survives: median surviving D is 17 nt for IGH (~5.7 aa) but 5 nt for human TRB
(~1.7 aa).

So neither source alone is enough, and which one dominates flips by locus:

    locus       prior only   aa only   combined     n     (held-out seed, generated)
    human IGH      15 %       81 %       82 %      345
    human TRB      76 %       70 %       82 %      595
    human TRD      86 %       88 %       87 %      699
    mouse TRB      76 %       83 %       85 %      699

"prior only" is ``beta = 0``; "aa only" is argmax of the match score under a uniform prior,
ties broken by marginal usage. Combining wins at IGH and both TRB; TRD is a wash, because
one D gene (TRDD3) accounts for 85 % of rearrangements and the aa match already finds it.

The combination is ``log P(D | M, J) + beta * s_D``: the length-and-J prior, tempered by
the best gapless local alignment score ``s_D`` of the D's three-frame translations against
the non-templated middle of the junction. ``beta`` is fitted per locus and shipped in
``database/vdj/<org>/d_prior.tsv`` with the distributions themselves, so nothing here needs
OLGA at runtime. It is flat above ~1.25 for TRB, so the shipped values are not delicate.

**Honesty about the numbers.** The table is measured on junctions drawn from the same
generative model that supplies the prior, so the prior's contribution is flattered. The
amino-acid contribution is not — it is germline matching. Rearrangements that genomic order
forbids are excluded from the truth: OLGA's human TRB model emits TRBD2 x TRBJ1 in 8.7 % of
draws, and scoring against those measures agreement with a model artifact.

Out of model, against nucleotide D calls (E <= 0.05) on the real GenBank fixtures: human
TRB 94 %, IGH 85 %, TRD 91 %, mouse TRB 85 %. On TRB, note that both this posterior and the
nucleotide caller enforce the same D2-x-J1 constraint, so their agreement on TRBJ1 records
is guaranteed rather than earned; the TRBJ2 rows, where both D genes stay possible, score
91 % (human) and 81 % (mouse).

Priors exist only for the (organism, locus) pairs with a published model: human IGH, TRB
and TRD, and mouse TRB. Everything else returns ``None`` rather than guessing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache

from . import _markup
from .annotate.reference import _load_d_germlines
from .cdr3fix import load_anchors, markup_cdr3, resolve_locus, resolve_species
from .paths import vdj_dir
from .refbuild.translate import translate

__all__ = ["DPosterior", "posterior_d", "load_d_prior"]


def _gene(allele: str) -> str:
    return allele.split("*")[0]


@dataclass
class DPrior:
    """The shipped generative-model summaries for one locus."""

    ins_vd: list[float]
    ins_dj: list[float]
    dlen: dict[str, list[float]]                 # allele -> P(surviving nt length)
    d_given_j: dict[str, dict[str, float]]       # j allele -> {d allele: P}
    d_marginal: dict[str, float]
    beta: float


@dataclass
class DPosterior:
    """Posterior over the D gene, and over where it sits in the junction."""

    locus: str
    d_call: str                                  # MAP D gene
    posterior: float                             # P(MAP gene)
    entropy: float                               # bits, over genes
    by_gene: dict[str, float] = field(default_factory=dict)
    support_aa: int = 0                          # best aa local-align score for the MAP gene
    d_start: int = -1                            # MAP nt offset of D start within the junction
    d_start_ci90: tuple[int, int] = (-1, -1)     # narrowest 90 % credible interval
    n_middle_nt: int = 0                         # insVD + |D| + insDJ, pinned by the length

    @property
    def confident(self) -> bool:
        """A hard call. 0.9 keeps ~the top decile of TRB records and most of IGH."""
        return self.posterior >= 0.9


@lru_cache(maxsize=8)
def load_d_prior(organism: str) -> dict[str, DPrior]:
    """``{locus: DPrior}``; empty when the organism has no shipped model."""
    path = vdj_dir(organism) / "d_prior.tsv"
    if not path.exists():
        return {}
    raw: dict[str, dict] = {}
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            locus, kind, key, value = parts
            raw.setdefault(locus, {}).setdefault(kind, {})[key] = float(value)

    out: dict[str, DPrior] = {}
    for locus, k in raw.items():
        ins_vd = _dense(k.get("insVD", {}))
        ins_dj = _dense(k.get("insDJ", {}))
        dlen: dict[str, dict[int, float]] = {}
        for key, p in k.get("dlen", {}).items():
            allele, _, length = key.rpartition(":")
            dlen.setdefault(allele, {})[int(length)] = p
        dgj: dict[str, dict[str, float]] = {}
        for key, p in k.get("d_given_j", {}).items():
            d_allele, _, j_allele = key.partition("|")
            dgj.setdefault(j_allele, {})[d_allele] = p
        out[locus] = DPrior(
            ins_vd=ins_vd, ins_dj=ins_dj,
            dlen={a: _dense(v) for a, v in dlen.items()},
            d_given_j=dgj, d_marginal=dict(k.get("d_marginal", {})),
            beta=float(next(iter(k.get("beta", {"": 1.0}).values()))),
        )
    return out


def _dense(sparse) -> list[float]:
    if not sparse:
        return []
    n = max(int(i) for i in sparse) + 1
    out = [0.0] * n
    for i, p in sparse.items():
        out[int(i)] = p
    return out


@lru_cache(maxsize=8)
def _d_aa_frames(organism: str) -> dict[str, tuple[str, str, str]]:
    """``{D allele: (aa frame 0, 1, 2)}`` from the shipped IMGT D germlines."""
    out: dict[str, tuple[str, str, str]] = {}
    for _locus, alleles in _load_d_germlines(vdj_dir(organism)).items():
        for allele, seq in alleles:
            out[allele] = tuple(translate(seq[f:], 0) for f in (0, 1, 2))
    return out


def _mask_forbidden(pd: dict[str, float], j_call: str) -> dict[str, float]:
    """Zero the D alleles lying 3' of the J, then renormalise.

    Mirrors :func:`arda.annotate.transfer._allowed_d`: TRBD2 sits downstream of the whole
    TRBJ1 cluster, so no deletional join can reach it. An ambiguous J spanning both clusters
    forbids nothing.
    """
    genes = {a.split("*")[0] for a in j_call.split(",") if a.strip()}
    if not genes or not all(g.startswith("TRBJ1-") for g in genes):
        return pd
    kept = {a: p for a, p in pd.items() if not a.startswith("TRBD2")}
    total = sum(kept.values())
    return {a: p / total for a, p in kept.items()} if total > 0 else pd


def _logsumexp(xs: list[float]) -> float:
    if not xs:
        return -math.inf
    m = max(xs)
    if m == -math.inf:
        return -math.inf
    return m + math.log(sum(math.exp(x - m) for x in xs))


def posterior_d(junction_aa: str, v_call: str, j_call: str,
                species: str = "human") -> DPosterior | None:
    """Posterior over the D gene (and its position) for an amino-acid junction.

    ``junction_aa`` is junction space (Cys104 .. Phe/Trp118, both included), as in
    :mod:`arda.cdr3fix`. Returns ``None`` when the locus has no D, no shipped model, or
    the junction cannot be marked up.
    """
    organism = resolve_species(species)
    locus = resolve_locus(v_call, j_call)
    prior = load_d_prior(organism).get(locus)
    if prior is None:
        return None

    mk = markup_cdr3(junction_aa, v_call, j_call, organism,
                     anchors=load_anchors(organism))
    if mk.v_end < 0 or mk.j_start < 0 or mk.j_start < mk.v_end:
        return None
    middle = mk.cdr3_repaired[mk.v_end : mk.j_start]
    n_middle = 3 * len(middle)                   # insVD + |D surviving| + insDJ, pinned

    frames = _d_aa_frames(organism)
    alleles = sorted(set(prior.dlen) | set(prior.d_marginal))
    if not alleles:
        return None

    # P(D allele | J): the load-bearing prior for TRB, where genomic order forbids TRBD2 x
    # TRBJ1 outright. Back off to the marginal when the J allele is outside the model -- but
    # the marginal pools both J clusters, so re-apply the same mask (the shipped human model
    # has no TRBJ1-6*01 row at all, and would otherwise let TRBD2 back in through the door).
    j_id = j_call.split(",")[0].strip()
    pd_given_j = prior.d_given_j.get(j_id) or _mask_forbidden(prior.d_marginal, j_call)

    # Joint over (D, insVD), marginalising the surviving-D length and insDJ.
    joint: dict[str, list[float]] = {}
    pa = [0.0] * max(1, len(prior.ins_vd))
    for allele in alleles:
        base = pd_given_j.get(allele, 0.0)
        if base <= 0:
            continue
        dl = prior.dlen.get(allele, [])
        row = [0.0] * len(prior.ins_vd)
        for a, p_a in enumerate(prior.ins_vd):
            if p_a <= 0:
                continue
            tot = 0.0
            for L, p_l in enumerate(dl):
                if p_l <= 0:
                    continue
                b = n_middle - a - L
                if 0 <= b < len(prior.ins_dj):
                    tot += p_l * prior.ins_dj[b]
            if tot > 0:
                row[a] = base * p_a * tot
                pa[a] += row[a]
        if any(row):
            joint[allele] = row
    total = sum(pa)
    if total <= 0:
        return None

    # Amino-acid evidence: best gapless local alignment of each D's three-frame
    # translations against the non-templated middle.
    score: dict[str, int] = {}
    for allele in joint:
        fr = frames.get(allele)
        score[allele] = 0 if not fr or not middle else max(
            (_markup.d_local_align(middle, f)[0] for f in fr if f), default=0)

    # Combine: log-prior + beta * aa score, marginalised over alleles within a gene.
    log_by_allele = {a: math.log(sum(joint[a])) + prior.beta * score[a] for a in joint}
    by_gene_log: dict[str, list[float]] = {}
    for allele, lp in log_by_allele.items():
        by_gene_log.setdefault(_gene(allele), []).append(lp)
    gene_log = {g: _logsumexp(v) for g, v in by_gene_log.items()}
    z = _logsumexp(list(gene_log.values()))
    by_gene = {g: math.exp(lp - z) for g, lp in gene_log.items()}
    best_gene = max(by_gene, key=by_gene.get)
    # `+ 0.0` so a degenerate posterior reports 0.0 rather than -0.0.
    entropy = -sum(p * math.log2(p) for p in by_gene.values() if p > 0) + 0.0

    # Where: d_start = (nt templated by V) + insVD, marginalising D and its trimming.
    pa_norm = [p / total for p in pa]
    order = sorted(range(len(pa_norm)), key=lambda i: pa_norm[i], reverse=True)
    cum, chosen = 0.0, []
    for i in order:
        chosen.append(i)
        cum += pa_norm[i]
        if cum >= 0.90:
            break
    v_nt = 3 * mk.v_end
    map_a = order[0] if order else -1
    ci = (v_nt + min(chosen), v_nt + max(chosen)) if chosen else (-1, -1)

    best_allele = max((a for a in log_by_allele if _gene(a) == best_gene),
                      key=lambda a: log_by_allele[a])
    return DPosterior(
        locus=locus, d_call=best_gene, posterior=by_gene[best_gene], entropy=entropy,
        by_gene=by_gene, support_aa=score.get(best_allele, 0),
        d_start=(v_nt + map_a) if map_a >= 0 else -1, d_start_ci90=ci,
        n_middle_nt=n_middle)
