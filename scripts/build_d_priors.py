#!/usr/bin/env python
"""Derive the D-inference priors used by ``arda.dpost`` from generative V(D)J models.

# Reads OLGA's shipped models (and the vdjrearm human TRD/TRG models, which OLGA lacks)
# and writes a small per-organism table so that the runtime posterior needs no OLGA
# dependency and no model files.
#
#     python scripts/build_d_priors.py            # writes database/vdj/<org>/d_prior.tsv
#
# 2026-07-10

What ships, per (locus):

* ``insVD`` / ``insDJ`` — P(number of inserted nt) on each side of D.
* ``dlen``  — P(surviving D nt length | D allele), marginalised over delDl/delDr.
* ``d_given_j`` — P(D allele | J allele). Load-bearing for TRB, where genomic order forbids
  half the matrix outright (see ``_forbidden``).
* ``beta`` — weight of the amino-acid match score relative to the log-prior, fitted per
  locus on generated data (see ``arda.dpost``).

The IGoR models are fitted without the genomic-order constraint, and the human TRB model
visibly fails to recover it: TRBD2*02 absorbs 21-27 % of every TRBJ1 row (TRBD2*01 correctly
falls to ~1e-5), so at gene level the unmasked model claims P(TRBD2 | TRBJ1) ~ 0.23 for a
join that cannot physically occur. The mouse TRB model does learn it (residual ~0.003). We
mask and renormalise rather than trust either.

Together with the junction length these pin ``insVD + |D| + insDJ``, which is what lets a
D be *placed* (to a median of 1-4 nt) even when the amino-acid sequence says nothing about
*which* D it is.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from arda.paths import vdj_dir  # noqa: E402

# Fitted on 700 generated junctions per model (train seed 101, held out on seed 202).
# beta = 0 ignores the amino acids; beta -> inf ignores the prior. See `arda.dpost`.
BETA = {("human", "IGH"): 3.25, ("human", "TRB"): 1.00,
        ("human", "TRD"): 1.75, ("mouse", "TRB"): 1.50}

MODELS = {
    "human_T_beta": ("human", "TRB"),
    "human_B_heavy": ("human", "IGH"),
    "mouse_T_beta": ("mouse", "TRB"),
}


def _model_dirs() -> dict[str, tuple[str, str]]:
    import olga

    base = Path(olga.__file__).parent / "default_models"
    out = {str(base / d): v for d, v in MODELS.items() if (base / d).is_dir()}
    vr = os.environ.get("ARDA_VDJREARM")
    if vr:
        trd = Path(vr) / "model" / "Homo+sapiens" / "TRD"
        if (trd / "model_params.txt").exists():
            out[str(trd)] = ("human", "TRD")
    return out


def _load(model_dir: str):
    import olga.load_model as lm

    gd = lm.GenomicDataVDJ()
    gd.load_igor_genomic_data(f"{model_dir}/model_params.txt",
                              f"{model_dir}/V_gene_CDR3_anchors.csv",
                              f"{model_dir}/J_gene_CDR3_anchors.csv")
    gm = lm.GenerativeModelVDJ()
    gm.load_and_process_igor_model(f"{model_dir}/model_marginals.txt")
    return gd, gm


def _forbidden(d_name: str, j_name: str) -> bool:
    """True when the D lies 3' of the J, so deletional joining cannot produce the pair.

    TRB genomic order is TRBD1 - TRBJ1 cluster - TRBC1 - TRBD2 - TRBJ2 cluster - TRBC2, so
    TRBD2 x TRBJ1 is unproducible. IGH and TRD put every D 5' of every J: nothing forbidden.
    """
    return d_name.startswith("TRBD2") and j_name.startswith("TRBJ1-")


def rows_for(locus: str, gd, gm) -> list[tuple[str, str, str, float]]:
    out: list[tuple[str, str, str, float]] = []
    d_names = [rec[0].strip() for rec in gd.genD]   # genJ/genV carry a third field
    j_names = [rec[0].strip() for rec in gd.genJ]

    for i, p in enumerate(gm.PinsVD):
        if p > 0:
            out.append((locus, "insVD", str(i), float(p)))
    for i, p in enumerate(gm.PinsDJ):
        if p > 0:
            out.append((locus, "insDJ", str(i), float(p)))

    # P(surviving nt length | D allele): cutD already carries the max palindrome on both ends.
    for di, name in enumerate(d_names):
        full = len(gd.cutD_genomic_CDR3_segs[di])
        acc = np.zeros(full + 1)
        pd = gm.PdelDldelDr_given_D[:, :, di]
        for dl in range(pd.shape[0]):
            for dr in range(pd.shape[1]):
                surviving = full - dl - dr
                if surviving >= 0:
                    acc[surviving] += pd[dl, dr]
        if acc.sum() <= 0:
            continue
        acc /= acc.sum()
        for L, p in enumerate(acc):
            if p > 0:
                out.append((locus, "dlen", f"{name}:{L}", float(p)))

    # P(D | J), after masking the pairs genomic order forbids and renormalising each J column.
    pdj = np.asarray(gm.PDJ, dtype=float)
    for di, dn in enumerate(d_names):
        for ji, jn in enumerate(j_names):
            if _forbidden(dn, jn):
                pdj[di, ji] = 0.0
    for ji, jn in enumerate(j_names):
        col = pdj[:, ji]
        if col.sum() <= 0:
            continue
        col = col / col.sum()
        for di, dn in enumerate(d_names):
            if col[di] > 0:
                out.append((locus, "d_given_j", f"{dn}|{jn}", float(col[di])))

    # Marginal P(D), used when the J allele is unknown to the model.
    marg = pdj.sum(axis=1)
    marg = marg / marg.sum()
    for di, dn in enumerate(d_names):
        if marg[di] > 0:
            out.append((locus, "d_marginal", dn, float(marg[di])))
    return out


def main() -> None:
    dirs = _model_dirs()
    if not dirs:
        raise SystemExit("no models found; pip install olga and set $ARDA_VDJREARM")
    by_org: dict[str, list[tuple[str, str, str, float]]] = {}
    for model_dir, (org, locus) in sorted(dirs.items(), key=lambda kv: kv[1]):
        gd, gm = _load(model_dir)
        rows = rows_for(locus, gd, gm)
        beta = BETA.get((org, locus))
        if beta is None:
            raise SystemExit(f"no fitted beta for {org}/{locus}")
        rows.append((locus, "beta", "", beta))
        by_org.setdefault(org, []).extend(rows)
        print(f"{org}/{locus}: {len(rows)} rows  (D={len(gd.genD)}, J={len(gd.genJ)}, beta={beta})")

    for org, rows in by_org.items():
        out = vdj_dir(org) / "d_prior.tsv"
        with open(out, "w") as fh:
            fh.write("locus\tkind\tkey\tvalue\n")
            for locus, kind, key, value in rows:
                fh.write(f"{locus}\t{kind}\t{key}\t{value:.8g}\n")
        print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
