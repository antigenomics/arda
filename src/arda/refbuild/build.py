"""Orchestrate the per-species reference database build.

For each locus: enumerate deduplicated V-J scaffolds, annotate them with IgBLAST,
keep those with complete FR1-FR4 + CDR1-3 markup, translate to protein, and
derive protein markup. Writes the committed artifacts under
``database/vdj/<organism>/`` plus a comprehensive ``build.log``.
"""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import polars as pl

from ..paths import data_dir, vdj_dir
from ..igblast import SUPPORTED_ORGANISMS
from .loci import LOCI, IMGT_SPECIES_DIR
from . import imgt, combinations, airr_extract
from .translate import translate, aa_coords_from_nt
from .airr_extract import REGION_NAMES

__all__ = ["build", "build_species"]


def _scaffold_fasta_path(species_dir: str, locus_name: str) -> Path:
    d = data_dir() / "scaffolds" / species_dir
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{locus_name}.fasta"


def _setup_logger(out_dir: Path) -> logging.Logger:
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(f"arda.build.{out_dir.name}")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(out_dir / "build.log", mode="w")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(fh)
    logger.addHandler(logging.StreamHandler())
    return logger


def _process_locus(organism, species_dir, locus, j_frames, logger):
    """Return (nt_rows, aa_rows, combo_rows, fasta_nt, fasta_aa) lists for a locus."""
    from ..igblast import has_internal_annotation
    if not has_internal_annotation(organism, locus.group):
        logger.info("%s: no IgBLAST %s internal annotation for %s — skipped",
                    locus.name, locus.group, organism)
        return [], [], [], [], []

    v = imgt.load_functional_alleles(species_dir, locus.group, locus.v)
    j = imgt.load_functional_alleles(species_dir, locus.group, locus.j)
    if not v or not j:
        logger.warning("%s: missing V (%d) or J (%d) alleles — skipped",
                       locus.name, len(v), len(j))
        return [], [], [], [], []

    scaffolds = combinations.build_locus_scaffolds(locus, v, j, j_frames)
    if not scaffolds:
        logger.warning("%s: no scaffolds produced — skipped", locus.name)
        return [], [], [], [], []
    by_id = {s.scaffold_id: s for s in scaffolds}

    fa = _scaffold_fasta_path(species_dir, locus.name)
    fa.write_text("".join(f">{s.scaffold_id}\n{s.sequence}\n" for s in scaffolds))

    df = airr_extract.annotate_scaffolds(
        fa, organism, species_dir, locus, num_threads=max(1, (os.cpu_count() or 2)))
    raw = len(v) * len(j)
    logger.info("%s: V=%d J=%d raw_combos=%d unique_scaffolds=%d annotated=%d",
                locus.name, len(v), len(j), raw, len(scaffolds), df.height)

    nt_rows, aa_rows, combo_rows, fasta_nt, fasta_aa = [], [], [], [], []
    incomplete = 0
    coord_cols = [f"{r}_start" for r in REGION_NAMES] + [f"{r}_end" for r in REGION_NAMES]

    for rec in df.iter_rows(named=True):
        sid = rec["sequence_id"]
        sc = by_id.get(sid)
        if sc is None:
            continue
        # Require complete markup (all region coordinates present).
        if any(rec.get(c) in (None, "", "NA") for c in coord_cols):
            incomplete += 1
            continue

        coding_start = int(rec["fwr1_start"])  # 1-based; scaffold reads frame 0
        protein = translate(sc.sequence, coding_start - 1)
        v_call = ",".join(sc.v_calls)
        j_call = ",".join(sc.j_calls)

        combo_rows.append({
            "scaffold_id": sid, "locus": locus.name,
            "v_calls": v_call, "j_calls": j_call, "n_pad": sc.n_pad,
        })
        fasta_nt.append((sid, sc.sequence))
        fasta_aa.append((sid, protein))

        nt_row = {"scaffold_id": sid, "locus": locus.name,
                  "v_call": v_call, "j_call": j_call,
                  "productive": rec.get("productive"),
                  # Extended markup: scaffold nt positions of the V germline end and
                  # J germline start, transferred to queries to locate the V/J split
                  # inside the junction (and to bridge frame for out-of-frame calls).
                  "v_sequence_end": rec.get("v_sequence_end") or "",
                  "j_sequence_start": rec.get("j_sequence_start") or "",
                  "junction": rec.get("junction"), "junction_aa": rec.get("junction_aa")}
        aa_row = {"scaffold_id": sid, "locus": locus.name,
                  "v_call": v_call, "j_call": j_call, "coding_start": coding_start,
                  "junction_aa": rec.get("junction_aa")}
        for r in REGION_NAMES:
            ns, ne = int(rec[f"{r}_start"]), int(rec[f"{r}_end"])
            nt_row[f"{r}_start"], nt_row[f"{r}_end"] = ns, ne
            nt_row[r] = rec.get(r)
            a_s, a_e = aa_coords_from_nt(ns, ne, coding_start)
            aa_row[f"{r}_start"], aa_row[f"{r}_end"] = a_s, a_e
            # Slice our own protein so aa coords round-trip exactly (igblast's
            # independently-translated *_aa can differ at ragged boundaries).
            aa_row[r] = protein[a_s - 1 : a_e]
        nt_rows.append(nt_row)
        aa_rows.append(aa_row)

    if incomplete:
        logger.info("%s: dropped %d scaffolds with incomplete markup", locus.name, incomplete)
    return nt_rows, aa_rows, combo_rows, fasta_nt, fasta_aa


def build_species(organism: str) -> Path:
    """Build the reference DB for one organism. Returns the output directory."""
    if organism not in IMGT_SPECIES_DIR:
        raise ValueError(f"Unknown organism {organism!r}; one of {list(IMGT_SPECIES_DIR)}")
    species_dir = IMGT_SPECIES_DIR[organism]
    out_dir = vdj_dir(organism)
    logger = _setup_logger(out_dir)
    t0 = time.perf_counter()
    logger.info("=== arda reference build: %s (%s) ===", organism, species_dir)

    imgt.download_reference()
    j_frames = combinations.load_j_frames(organism)

    nt_all, aa_all, combo_all, fa_nt, fa_aa = [], [], [], [], []
    for locus in LOCI:
        try:
            nt, aa, combo, fnt, faa = _process_locus(
                organism, species_dir, locus, j_frames, logger)
        except Exception as exc:  # noqa: BLE001 — one bad locus must not kill the species
            logger.warning("%s: failed (%s) — skipped", locus.name, exc)
            continue
        nt_all += nt; aa_all += aa; combo_all += combo; fa_nt += fnt; fa_aa += faa

    # Write artifacts.
    (out_dir / "alleles.fasta").write_text("".join(f">{i}\n{s}\n" for i, s in fa_nt))
    (out_dir / "alleles.aa.fasta").write_text("".join(f">{i}\n{s}\n" for i, s in fa_aa))
    pl.DataFrame(nt_all).write_csv(out_dir / "markup.tsv", separator="\t")
    pl.DataFrame(aa_all).write_csv(out_dir / "markup.aa.tsv", separator="\t")
    pl.DataFrame(combo_all).write_csv(out_dir / "combinations.tsv", separator="\t")

    dt = time.perf_counter() - t0
    logger.info("TOTAL: %d scaffolds across %d loci in %.1fs -> %s",
                len(fa_nt), len({r['locus'] for r in combo_all}), dt, out_dir)
    return out_dir


def build(organism: str = "all") -> None:
    """Build one organism or ``"all"`` supported organisms."""
    organisms = SUPPORTED_ORGANISMS if organism == "all" else (organism,)
    for org in organisms:
        build_species(org)
