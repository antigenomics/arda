"""The V -> N1 -> D -> N2 -> J semi-Markov model, read as inference.

:mod:`arda.scenarios` estimates this model's parameters by EM. This module is the *same
recursion* used the other way round: given parameters, score a junction and ask what it was made
of. :func:`arda.scenarios.lattice` is the forward-backward pass and there is one implementation
of it, because an E-step and a posterior differ only in what you do with the same term weights.

**Semi-Markov, not Markov**, because the durations are explicit and non-geometric: germline
deletions and insertion lengths each have their own tabulated distribution -- exactly the shape
``database/vdj/<org>/d_prior.tsv`` already ships. A plain HMM would impose geometric durations on
quantities that are measured not to be.

**Conditioned on the V and J call**, which mmseqs already made. That is what keeps the state
space to ``(delV, insVD, D, delDl, delDr, insDJ, delJ)`` and the recursion cheap -- per clonotype
after ``correct``, never per read.

How this relates to :mod:`arda.dpost`
-------------------------------------

It does not replace it. ``dpost`` answers the **amino-acid** question -- a VDJdb-style record with
no nucleotides, where the D is often invisible in the translated junction and the length
constraint plus an aa match is all there is. This module needs nucleotides. They are two
different inputs, and both ship.

``shm`` threads an :class:`arda.shmmodel.ShmModel` into the same recursion: without one the
templated V length is bounded by an exact common prefix, with one it is priced by the model, so a
hypermutated V tail stops being read as N-region. Default ``None``, which is what every shipped
caller passes and is byte-identical to the behaviour before it existed.

Never: **do not re-run the two negatives the roadmap already records.** (i) Re-ranking nucleotide
D candidates by a scenario likelihood changes nothing -- 98.9 -> 97.8 % gene accuracy on IGH,
94.2 -> 94.5 % on huTRB, flat elsewhere, because with 10-18 matched nt the alignment term is
11-20 nats and the prior moves it by 3. (ii) Replacing ``_D_MAX_EVALUE`` with a present/absent
Bayes factor buys IGH ~+3 pp recall at matched FP and nothing for TRD, but needs a *per-locus*
threshold -- which is a shipped constant of exactly the kind this project refuses. So this module
reports likelihoods and posteriors; it does not gate anything, and nothing in the annotation path
calls it.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log

from .scenarios import Scenario, _Model, lattice

__all__ = ["DPosterior", "model_for", "log_likelihood", "posterior_d", "best_scenario"]


@dataclass(frozen=True, slots=True)
class DPosterior:
    """``P(D | junction, V, J)`` and the evidence behind it."""

    locus: str
    probabilities: dict[str, float]
    log_likelihood: float
    scenarios: int

    @property
    def best(self) -> str:
        """The most probable D allele, ties broken lexicographically for determinism."""
        if not self.probabilities:
            return ""
        top = max(self.probabilities.values())
        return min(a for a, p in self.probabilities.items() if p == top)


def model_for(organism: str = "human", *, prior: str | None = None) -> _Model:
    """The parameter set to score against.

    ``prior`` names a TSV written by ``arda scenarios`` -- so a model estimated on one cohort can
    score another. ``None`` uses the shipped ``d_prior.tsv``, whose numbers are OLGA's.
    """
    model = _Model(organism)
    if prior is not None:
        _load_into(model, prior, organism)
    return model


def _load_into(model: _Model, path, organism: str) -> None:
    """Overwrite the shipped parameters with an ``arda scenarios`` table."""
    from .scenarios import _dist

    raw: dict[tuple[str, str], dict[str, float]] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4 or parts[0] == "locus":
                continue
            raw.setdefault((parts[0], parts[1]), {})[parts[2]] = float(parts[3])
    for (locus, kind), table in raw.items():
        if kind in ("insVD", "insDJ"):
            dense = _dist({int(k): v for k, v in table.items()}, model.max_ins + 1)
            (model.ins_vd if kind == "insVD" else model.ins_dj)[locus] = dense
        elif kind in ("dlen", "delV", "delJ"):
            dst = {"dlen": model.dlen, "delV": model.del_v, "delJ": model.del_j}[kind]
            per: dict[str, dict[int, float]] = {}
            for key, v in table.items():
                allele, _, n = key.rpartition(":")
                per.setdefault(allele, {})[int(n)] = v
            dst.setdefault(locus, {}).update(
                {a: _dist(sp, max(sp) + 1) for a, sp in per.items()})
        elif kind == "d_given_j":
            for key, v in table.items():
                d_allele, _, j_allele = key.partition("|")
                model.d_given_j.setdefault(locus, {}).setdefault(j_allele, {})[d_allele] = v


def log_likelihood(junction_nt: str, v_call: str, j_call: str, model: _Model,
                   species: str = "human", shm=None) -> float:
    """``log P(junction | V, J)``, marginalised over every scenario. ``-inf`` if unscoreable.

    Marginal, not the best path. A junction explainable many mediocre ways is more probable than
    one explainable a single slightly-better way, and a Viterbi score cannot say so.
    """
    built = lattice(junction_nt, v_call, j_call, model, species, shm)
    if built is None:
        return float("-inf")
    z = sum(t[0] for t in built[2])
    return log(z) if z > 0 else float("-inf")


def posterior_d(junction_nt: str, v_call: str, j_call: str, model: _Model,
                species: str = "human", shm=None) -> DPosterior:
    """``P(D | junction, V, J)``, summing over trimming, insertions and placement.

    Never: this sums over scenarios rather than taking the best one. A D that fits many
    placements passably is better supported than one that fits a single placement well, and the
    trimming/insertion nuisance parameters are exactly what has to be integrated out to see that.
    """
    built = lattice(junction_nt, v_call, j_call, model, species, shm)
    if built is None:
        return DPosterior("", {}, float("-inf"), 0)
    g, _L, terms, _left, _right = built
    z = sum(t[0] for t in terms)
    if z <= 0:
        return DPosterior(g.locus, {}, float("-inf"), 0)
    probs: dict[str, float] = {}
    for w, _p, _q, allele, _m in terms:
        if allele:
            probs[allele] = probs.get(allele, 0.0) + w / z
    return DPosterior(g.locus, probs, log(z), len(terms))


def best_scenario(junction_nt: str, v_call: str, j_call: str, model: _Model,
                  species: str = "human", shm=None) -> Scenario | None:
    """The single most probable scenario -- the Viterbi path, for reading, not for counting.

    Never: this is NOT what the estimator counts, and ``project/design-scenarios.md`` says why.
    One reading carries weight 1 and the rest 0, which biases trimming short and insertions
    shorter. It is here because a human inspecting one junction wants one answer.
    """
    built = lattice(junction_nt, v_call, j_call, model, species, shm)
    if built is None:
        return None
    g, L, terms, left, right = built
    if not terms:
        return None

    def flank(parts):
        # (germline_len, ins_len) of the heaviest split of this flank.
        return max(parts, key=lambda t: t[2])[:2] if parts else (0, 0)

    if not g.d_germlines:
        # VJ locus: terms are (w, v_templated, j_templated, "", -1).
        w, v_t, j_t, _a, _m = max(terms, key=lambda t: (t[0], -t[1], -t[2]))
        return Scenario(del_v=len(g.v_nt) - v_t, ins_vd=L - v_t - j_t,
                        del_j=len(g.j_nt) - j_t)
    # Deterministic: ties broken on (allele, start), never on term order.
    w, p, q, allele, m = max(terms, key=lambda t: (t[0], [-ord(c) for c in t[3]], -t[1]))
    v_t, ins1 = flank(left[p][1])
    j_t, ins2 = flank(right[q][1])
    dg = dict(g.d_germlines)[allele]
    # The germline offset is recovered from the placement: the surviving D sits at `p`, and among
    # equal-length matches `_d_placements` kept the first offset, which is what this reproduces.
    start = dg.find(junction_nt.strip().upper()[p:p + m]) if m else 0
    return Scenario(del_v=len(g.v_nt) - v_t, ins_vd=ins1, d_call=allele,
                    del_dl=max(start, 0), del_dr=len(dg) - max(start, 0) - m,
                    ins_dj=ins2, del_j=len(g.j_nt) - j_t)
