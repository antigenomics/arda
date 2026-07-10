"""Ground-truth junction generator built on OLGA's generative model.

``SequenceGeneration*.gen_rnd_prod_CDR3()`` returns only ``(nt, aa, V, J)`` and
throws away the recombination events, but ``choose_random_recomb_events()`` keeps
them. So we re-run OLGA's own accept/reject loop (OLGA is GPL-3.0, as is arda) and
retain ``delV``/``delJ``/``delDl``/``delDr``/``insVD``/``insDJ``. That gives the
exact number of germline nucleotides each segment contributed to the junction --
the only ground truth that can validate arda's error-vs-boundary heuristic.

Two OLGA details matter and are easy to get wrong:

* ``cutV`` is the germline CDR3 segment with ``max_delV_palindrome`` **palindromic**
  nucleotides appended on the right; ``cutJ`` has them prepended on the left. Those
  P-nucleotides are not germline, so the germline-derived length must be computed
  net of them.
* OLGA's germline sequences are **not** IMGT's (its ``TRBV3-1*01`` starts three
  nucleotides later, which is why its CDR3 anchor is 267 where IMGT's is 270).
  ``Truth.germline_matches_imgt`` says whether the two agree for this record; a
  ground-truth comparison is only meaningful where they do.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# Model dirs shipped by OLGA, mapped to (arda organism, locus, is_vdj).
OLGA_MODELS = {
    "human_T_alpha": ("human", "TRA", False),
    "human_T_beta": ("human", "TRB", True),
    "human_B_heavy": ("human", "IGH", True),
    "human_B_kappa": ("human", "IGK", False),
    "human_B_lambda": ("human", "IGL", False),
    "mouse_T_alpha": ("mouse", "TRA", False),
    "mouse_T_beta": ("mouse", "TRB", True),
}

# OLGA ships no TRG/TRD. `$ARDA_VDJREARM` supplies them (human only); TRD carries a D.
VDJREARM_MODELS = {
    "TRD": ("human", "TRD", True),
    "TRG": ("human", "TRG", False),
}


@dataclass
class Truth:
    """One generated junction plus what each germline actually contributed."""

    cdr3_nt: str
    cdr3_aa: str
    v_call: str
    j_call: str
    v_germline_nt: int          # germline (non-palindromic) nt the V contributed
    j_germline_nt: int          # germline (non-palindromic) nt the J contributed
    d_call: str = ""
    d_start: int = -1           # 0-based, inclusive, within cdr3_nt
    d_end: int = -1
    germline_matches_imgt: bool = True
    # Full event record: needed to reason about where D *must* be, from the
    # insertion-length model alone (junction length pins insVD + |D| + insDJ).
    v_nt: int = 0               # nt the V contributed, palindromes included
    j_nt: int = 0
    ins_vd: int = 0
    ins_dj: int = 0
    d_nt: int = 0               # surviving D nt (0 = trimmed away entirely)
    v_idx: int = -1
    j_idx: int = -1
    d_idx: int = -1

    @property
    def v_end(self) -> int:
        """Junction residues wholly templated by the V germline."""
        return self.v_germline_nt // 3

    @property
    def j_start(self) -> int:
        """Index of the first junction residue wholly templated by the J germline."""
        return len(self.cdr3_aa) - (self.j_germline_nt // 3)


def olga_model_dirs() -> dict[str, tuple[str, str, bool]]:
    """``{dir: (organism, locus, is_vdj)}`` for every model we can actually load."""
    import olga

    base = Path(olga.__file__).parent / "default_models"
    out = {str(base / d): v for d, v in OLGA_MODELS.items() if (base / d).is_dir()}
    vr = os.environ.get("ARDA_VDJREARM")
    if vr:
        root = Path(vr) / "model" / "Homo+sapiens"
        for d, v in VDJREARM_MODELS.items():
            if (root / d / "model_params.txt").exists():
                out[str(root / d)] = v
    return out


def load_generator(model_dir: str, is_vdj: bool):
    import olga.load_model as lm
    import olga.sequence_generation as sg

    gd = lm.GenomicDataVDJ() if is_vdj else lm.GenomicDataVJ()
    gd.load_igor_genomic_data(f"{model_dir}/model_params.txt",
                              f"{model_dir}/V_gene_CDR3_anchors.csv",
                              f"{model_dir}/J_gene_CDR3_anchors.csv")
    gm = lm.GenerativeModelVDJ() if is_vdj else lm.GenerativeModelVJ()
    gm.load_and_process_igor_model(f"{model_dir}/model_marginals.txt")
    gen = sg.SequenceGenerationVDJ(gm, gd) if is_vdj else sg.SequenceGenerationVJ(gm, gd)
    return gd, gen


def generate(model_dir: str, is_vdj: bool, n: int, anchors: dict, seed: int = 0,
             max_tries: int = 200_000) -> list[Truth]:
    """Generate ``n`` productive junctions with their germline contributions."""
    import numpy as np
    import olga.sequence_generation as sg
    from olga.utils import nt2aa

    np.random.seed(seed)
    gd, gen = load_generator(model_dir, is_vdj)
    out: list[Truth] = []
    for _ in range(max_tries):
        if len(out) >= n:
            break
        e = gen.choose_random_recomb_events()
        v_seg, j_seg = gd.genV[e["V"]][1], gd.genJ[e["J"]][1]
        cut_v, cut_j = gd.cutV_genomic_CDR3_segs[e["V"]], gd.cutJ_genomic_CDR3_segs[e["J"]]
        if len(cut_v) <= max(e["delV"], 0) or len(cut_j) < e["delJ"]:
            continue
        v_nt = cut_v[: len(cut_v) - e["delV"]]
        j_nt = cut_j[e["delJ"]:]

        if is_vdj:
            cut_d = gd.cutD_genomic_CDR3_segs[e["D"]]
            if len(cut_d) < e["delDl"] + e["delDr"]:
                continue
            d_nt = cut_d[e["delDl"]: len(cut_d) - e["delDr"]]
            if (len(v_nt) + len(d_nt) + len(j_nt) + e["insVD"] + e["insDJ"]) % 3:
                continue
            ins_vd = sg.rnd_ins_seq(e["insVD"], gen.C_Rvd, gen.C_first_nt_bias_insVD)
            ins_dj = sg.rnd_ins_seq(e["insDJ"], gen.C_Rdj, gen.C_first_nt_bias_insDJ)[::-1]
            nt = v_nt + ins_vd + d_nt + ins_dj + j_nt
            d_start = len(v_nt) + len(ins_vd)
            d_call, d_span = gd.genD[e["D"]][0].strip(), (d_start, d_start + len(d_nt) - 1)
            n_ins_vd, n_ins_dj, n_d = e["insVD"], e["insDJ"], len(d_nt)
        else:
            if (len(v_nt) + len(j_nt) + e["insVJ"]) % 3:
                continue
            ins_vj = sg.rnd_ins_seq(e["insVJ"], gen.C_Rvj, gen.C_first_nt_bias_insVJ)
            nt = v_nt + ins_vj + j_nt
            d_call, d_span = "", (-1, -1)
            n_ins_vd, n_ins_dj, n_d = e["insVJ"], 0, 0

        aa = nt2aa(nt)
        if "*" in aa or not aa or aa[0] != "C" or aa[-1] not in "FVW":
            continue

        # Strip the palindromic nucleotides: they are not germline-templated.
        v_germ = min(len(v_nt), len(v_seg))
        j_germ = len(j_nt) - max(0, gd.max_delJ_palindrome - e["delJ"])
        v_name, j_name = gd.genV[e["V"]][0], gd.genJ[e["J"]][0]
        va, ja = anchors.get(("V", v_name)), anchors.get(("J", j_name))
        agrees = bool(
            va and ja and va.status == "ok" and ja.status == "ok"
            and va.germline_nt.upper() == v_seg.upper()
            and ja.germline_nt.upper().endswith(j_seg.upper()))
        out.append(Truth(cdr3_nt=nt, cdr3_aa=aa, v_call=v_name, j_call=j_name,
                         v_germline_nt=max(0, v_germ), j_germline_nt=max(0, j_germ),
                         d_call=d_call, d_start=d_span[0], d_end=d_span[1],
                         germline_matches_imgt=agrees,
                         v_nt=len(v_nt), j_nt=len(j_nt), ins_vd=n_ins_vd, ins_dj=n_ins_dj,
                         d_nt=n_d, v_idx=e["V"], j_idx=e["J"],
                         d_idx=e.get("D", -1)))
    return out
