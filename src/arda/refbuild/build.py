"""Orchestrate the per-species reference database build.

For each locus: enumerate deduplicated V-J scaffolds, annotate them with IgBLAST,
keep those with complete FR1-FR4 + CDR1-3 markup, translate to protein, and
derive protein markup. Writes the committed artifacts under
``database/vdj/<organism>/`` plus a comprehensive ``build.log``.
"""

from __future__ import annotations

import logging
import os
import re
import time
from collections import Counter
from pathlib import Path

import polars as pl

from ..paths import data_dir, vdj_dir
from ..igblast import SUPPORTED_ORGANISMS
from .loci import LOCI, VDJ_LOCI, IMGT_SPECIES_DIR
from . import imgt, combinations, airr_extract, constant
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


def _process_locus(organism, species_dir, locus, j_frames, logger,
                   fwr3_stops=None, fr4_offsets=None, one_allele_per_gene=False):
    """Return (nt_rows, aa_rows, combo_rows, fasta_nt, fasta_aa) lists for a locus."""
    from ..igblast import has_internal_annotation
    if not has_internal_annotation(organism, locus.group):
        logger.info("%s: no IgBLAST %s internal annotation for %s — skipped",
                    locus.name, locus.group, organism)
        return [], [], [], [], []

    fwr3_stops = combinations.load_v_fwr3_stops(organism) if fwr3_stops is None else fwr3_stops
    fr4_offsets = combinations.load_j_fr4_offsets(organism) if fr4_offsets is None else fr4_offsets

    v = imgt.load_functional_alleles(species_dir, locus.group, locus.v)
    if locus.v_shared:
        # TRD shares V genes with TRA (chr14 interleaved locus): the shared genes are filed under TRAV
        # with a ".../DV..." designation. Add them so TRAV/DV δ rearrangements match a TRD scaffold
        # instead of being miscalled TRA. Safe: the J region (TRDJ vs TRAJ) disambiguates at runtime.
        stem, needle = locus.v_shared
        shared = imgt.load_functional_alleles(species_dir, locus.group, stem)
        added = {a: s for a, s in shared.items() if needle in a}
        v.update(added)
        logger.info("%s: +%d shared V alleles from %s (%s)", locus.name, len(added), stem, needle)
    j = imgt.load_functional_alleles(species_dir, locus.group, locus.j)
    if not v or not j:
        logger.warning("%s: missing V (%d) or J (%d) alleles — skipped",
                       locus.name, len(v), len(j))
        return [], [], [], [], []

    # Drop germlines that cannot produce a trustworthy junction. Must run after the v_shared merge
    # (the TRAV/DV alleles are subject to both rules too) and before the scaffold cross-product.
    v_tails = {a: v_cys_tail(s, fwr3_stops.get(a)) for a, s in v.items()}
    j_tails = {a: j_junction_tail(s, fr4_offsets.get(a)) for a, s in j.items()}
    #   (a) anchor unfindable, but a sibling allele proves what the junction should be
    v = drop_unanchorable(v, v_tails, logger=logger, locus=locus.name, segment="V")
    j = drop_unanchorable(j, j_tails, logger=logger, locus=locus.name, segment="J")
    #   (b) anchor found, but the germline stops short of its own gene's best allele
    v = drop_truncated(v, v_tails, logger=logger, locus=locus.name, segment="V")
    j = drop_truncated(j, j_tails, logger=logger, locus=locus.name, segment="J")
    if one_allele_per_gene:
        v = select_one_allele_per_gene(v, logger=logger, locus=locus.name, segment="V")
        j = select_one_allele_per_gene(j, logger=logger, locus=locus.name, segment="J")

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
    incomplete = unanchored = 0
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

        junction = rec.get("junction") or ""
        junction_aa = rec.get("junction_aa") or ""
        productive = rec.get("productive")
        if not junction_aa.startswith("C"):
            # A junction that does not begin at Cys104 is not a junction. It happens when the V
            # germline has no findable anchor: IgBLAST ships no `ndm` entry for the allele and the
            # FR3-motif fallback fails (9% of functional mouse V alleles, 5% of human). Such a
            # scaffold used to emit a garbage junction AND claim productive=T.
            #
            # Keep the scaffold -- the read still MAPS, so filter recall is untouched -- but emit no
            # junction and do not claim it is productive, so `correct --complete-only` drops it
            # instead of minting a bogus clonotype. Dropping the allele outright would delete 41
            # *functional* mouse IGHV/IGLV genes.
            junction, junction_aa, productive = "", "", "F"
            unanchored += 1

        combo_rows.append({
            "scaffold_id": sid, "locus": locus.name,
            "v_calls": v_call, "j_calls": j_call, "n_pad": sc.n_pad,
        })
        fasta_nt.append((sid, sc.sequence))
        fasta_aa.append((sid, protein))

        nt_row = {"scaffold_id": sid, "locus": locus.name,
                  "v_call": v_call, "j_call": j_call,
                  "productive": productive,
                  # Extended markup: scaffold nt positions of the V germline end and
                  # J germline start, transferred to queries to locate the V/J split
                  # inside the junction (and to bridge frame for out-of-frame calls).
                  "v_sequence_end": rec.get("v_sequence_end") or "",
                  "j_sequence_start": rec.get("j_sequence_start") or "",
                  "junction": junction, "junction_aa": junction_aa}
        aa_row = {"scaffold_id": sid, "locus": locus.name,
                  "v_call": v_call, "j_call": j_call, "coding_start": coding_start,
                  "junction_aa": junction_aa}
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
    if unanchored:
        logger.info("%s: %d scaffolds have no Cys104 anchor — kept for mapping, junction suppressed",
                    locus.name, unanchored)
    return nt_rows, aa_rows, combo_rows, fasta_nt, fasta_aa


def _collect_d_germlines(species_dir: str, logger) -> list[tuple[str, str, str]]:
    """Return ``(locus, allele, ungapped_seq)`` for every D allele of each VDJ locus.

    D segments map at runtime via gapless local alignment (not the V·J scaffold
    DB), so we ship the raw germline sequences alongside the scaffolds. Loci with
    no IMGT D file for this species (e.g. rat TRD) are skipped.
    """
    out: list[tuple[str, str, str]] = []
    for locus in VDJ_LOCI:
        try:
            path = imgt.ungap_gene(species_dir, locus.group, locus.d)
        except FileNotFoundError:
            logger.info("%s: no IMGT D file (%s) for %s — no D germlines",
                        locus.name, locus.d, species_dir)
            continue
        for header, seq in imgt.read_fasta(path):
            allele = header.split("|")[0].strip().split()[0]
            s = seq.upper().replace(".", "")
            if allele and s:
                out.append((locus.name, allele, s))
    return out


# --------------------------------------------------------------------------
# CDR3 junction anchors (per V/J allele).
#
# These let us mark up a bare (junction_aa, V, J) record -- the VDJdb case --
# without any alignment: we ship the residues each germline templates into the
# junction. Coordinates are JUNCTION space (Cys104..Phe/Trp118 inclusive), which
# is what VDJdb's `cdr3` column actually holds -- NOT the IMGT CDR3 (105..117).
#
# V anchor: IgBLAST's `.ndm.imgt` FWR3 stop (FR3-IMGT ends at 104 = the 2nd-CYS),
# so the Cys104 codon starts at `fwr3_stop - 3`. IgBLAST ships only a subset of
# IMGT alleles, so the rest fall back to the conserved FR3 motif (below).
#
# NOTE: the tempting shortcut "IMGT position 104 == gapped nt 310..312" is WRONG.
# It holds for human and rabbit but not for mouse/rhesus, whose V-QUEST gapped
# FASTAs carry insertion positions that break the linear slot arithmetic (rhesus
# IGHV1-111*01 has its Cys at gapped slot 106). It silently produced 671 wrong
# rhesus anchors before being caught by the Cys check. Do not reintroduce it.
#
# J anchor: IgBLAST's aux `cdr3_stop` is the last nt of position 117, so the
# [FW]118 codon starts at cdr3_stop + 1.
# --------------------------------------------------------------------------

# Conserved FR3 3' motif: the 2nd-CYS is preceded two residues back by an aromatic
# ("YYC", "YFC", "YHC"...). Used when IgBLAST has no entry for an allele. The
# 5'-most such Cys is taken because the V's CDR3 tail can hold a *second* Cys
# (`YYC-AC-DT` in TRDV2*01, `YYCC...` in IGLV2-11*01). Verified against every
# allele IgBLAST does annotate: 1183 agree, 0 disagree.
_FR3_AROMATIC = frozenset("YFHW")
_V_TAIL_MIN, _V_TAIL_MAX = 1, 14   # residues the V templates into the junction after Cys104

# Conserved J FR4 motif [FW]-G-X-G at IMGT 118..121. Used only when the aux file
# has no entry for an allele (rat/rhesus ship no TR aux). Verified against the
# human aux: 117 alleles anchored, 0 disagreements. It does not fire for the
# non-canonical J alleles (TRAJ35*01 has Cys118, TRBJ2-7*02 has Val118,
# TRAJ16*01 is FARG, TRAJ61*01 is FGAN) -- those are flagged `no_anchor` rather
# than guessed, because a wrong anchor silently corrupts every emitted coordinate.
_J_FR4_MOTIF = re.compile(r"[FW]G.G")

_ANCHOR_COLUMNS = ["locus", "segment", "allele", "functionality", "anchor_nt",
                   "partial_nt", "templated_aa", "germline_nt", "status", "source"]


def _gapped_alleles(species_dir: str, group: str, stem: str) -> list[tuple[str, str, str]]:
    """``(allele, functionality, gapped_seq)`` from an IMGT V-QUEST gene file.

    Keeps every allele regardless of functionality: VDJdb cites pseudogene and
    ORF alleles, and ``markup.tsv`` (which is scaffold-gated and deduped) is not
    a valid substitute.
    """
    out: list[tuple[str, str, str]] = []
    for header, seq in imgt.read_fasta(imgt.gene_fasta_path(species_dir, group, stem)):
        fields = header.split("|")
        if len(fields) < 4 or not seq:
            continue
        allele = fields[1].strip()
        func = fields[3].strip().strip("()[]").split("/")[0].split()[0] or "?"
        out.append((allele, func, seq.upper()))
    return out


def _v_anchor(ungapped: str, fwr3_stop: int | None) -> tuple[int, str]:
    """``(anchor_nt, source)`` of the Cys104 codon in an ungapped V germline.

    ``source`` is ``ndm`` (IgBLAST), ``motif`` (conserved FR3 aromatic), or
    ``no_anchor``. Every returned anchor is verified to be a Cys codon -- an
    anchor that is off by one codon silently corrupts every coordinate we emit,
    so we would rather refuse than guess.
    """
    if fwr3_stop:
        a = fwr3_stop - 3
        if a >= 0 and translate(ungapped[a : a + 3], 0) == "C":
            return a, "ndm"
    prot = translate(ungapped, 0)
    for i, res in enumerate(prot):       # 5'-most qualifying Cys
        if res != "C" or i < 60 or i < 2:
            continue
        if not (_V_TAIL_MIN <= len(prot) - 1 - i <= _V_TAIL_MAX):
            continue
        if prot[i - 2] in _FR3_AROMATIC:
            return 3 * i, "motif"
    return -1, "no_anchor"


def _j_anchor_from_motif(seq: str) -> int:
    """0-based nt offset of the [FW]118 codon via the FR4 motif, or ``-1``.

    Picks the frame with the fewest stop codons; only the first motif hit in a
    frame is considered.
    """
    best: tuple[int, int, int] | None = None
    for frame in (0, 1, 2):
        aa = translate(seq[frame:], 0)
        hit = _J_FR4_MOTIF.search(aa)
        if hit is None:
            continue
        cand = (aa.count("*"), frame, frame + 3 * hit.start())
        if best is None or cand < best:
            best = cand
    return best[2] if best else -1


_TAIL_SLACK_NT = 3   # one codon


def v_cys_tail(ungapped: str, fwr3_stop: int | None) -> int:
    """nt this V templates INTO the junction (Cys104 -> its 3' end); ``-1`` if the anchor is unfindable."""
    anchor, _ = _v_anchor(ungapped, fwr3_stop)
    return len(ungapped) - anchor if anchor >= 0 else -1


def j_junction_tail(seq: str, fr4_offset=None) -> int:
    """nt this J templates INTO the junction (its 5' end -> [FW]118); ``-1`` if unfindable.

    Mirror of :func:`v_cys_tail`: the V reaches into the junction from the left, the J from the right.
    """
    anchor = fr4_offset[0] + 1 if fr4_offset else _j_anchor_from_motif(seq)
    return anchor + 3 if 0 <= anchor and anchor + 3 <= len(seq) else -1


def drop_unanchorable(alleles: dict[str, str], tails: dict[str, int], *,
                      logger, locus: str, segment: str) -> dict[str, str]:
    """Drop alleles whose junction anchor is unfindable **when a sibling allele proves the truth**.

    ``tails[a] < 0`` means arda cannot locate Cys104 in this germline. The scaffold still gets a CDR3
    -- IgBLAST annotates one -- but at a position arda has no way to verify, and it is demonstrably
    the WRONG one: ``TRAV23/DV6*04`` yields a scaffold junction ``CTTSGTYKYIF`` (11 aa) while every
    other allele of its gene templates ``CAAS`` (13 nt, 14 aa junctions). Worse, those scaffolds are
    read magnets -- ``TRAV23/DV6*04`` took 55 % of all TRA reads in a real tumour, every one with a
    junction 3 aa short. Blanking the scaffold's ``junction_aa`` is not enough: the runtime projects
    ``cdr3_start``/``cdr3_end`` onto the read and re-derives a junction from those.

    Only drop when the gene keeps at least one anchored allele. Then the drop is FREE (no gene is
    lost, in any organism) and the surviving sibling is what proves the dropped one wrong. A gene
    whose alleles are ALL unanchorable is a different case -- nothing contradicts its annotation, and
    dropping it would delete 41 *functional* mouse IGHV/IGLV genes -- so those are left alone.
    """
    by_gene: dict[str, list[str]] = {}
    for a in alleles:
        by_gene.setdefault(a.split("*")[0], []).append(a)
    drop = {a for al in by_gene.values() if any(tails.get(x, -1) >= 0 for x in al)
            for a in al if tails.get(a, -1) < 0}
    if drop:
        logger.info("%s: dropped %d unanchorable %s allele(s) with an anchored sibling: %s%s",
                    locus, len(drop), segment, ", ".join(sorted(drop)[:6]),
                    " …" if len(drop) > 6 else "")
    return {a: s for a, s in alleles.items() if a not in drop}


def drop_truncated(alleles: dict[str, str], tails: dict[str, int], *,
                   logger, locus: str, segment: str) -> dict[str, str]:
    """Drop alleles whose germline is truncated INTO the junction, relative to their own gene.

    A 3'-truncated V (or 5'-truncated J) builds a scaffold that **lacks nucleotides a real
    rearrangement has**, so every read projected onto it yields a junction short by exactly that
    much -- and the junction still starts with Cys and ends with [FW], so it looks canonical and
    nothing downstream catches it. ``TRAV20*03`` templates 3 nt where its siblings template 13:
    every clonotype built on it came out 3 aa short, which is why arda's TRA CDR3 length peaked
    3 aa below every other tool's.

    Orphan-safe by construction: the threshold is relative to the gene's own longest allele, so that
    allele always survives -- no gene can lose all of its alleles.

    Alleles with no anchor at all (``tail < 0``) are deliberately NOT dropped here. They keep their
    scaffold so the read still maps; :func:`_process_locus` suppresses their junction instead.
    Dropping them would delete 41 *functional* mouse IGHV/IGLV genes.
    """
    by_gene: dict[str, list[str]] = {}
    for a, t in tails.items():
        if t >= 0:
            by_gene.setdefault(a.split("*")[0], []).append(a)
    drop = {a for al in by_gene.values()
            for a in al if max(tails[x] for x in al) - tails[a] >= _TAIL_SLACK_NT}
    if drop:
        shown = ", ".join(f"{a}({tails[a]}nt)" for a in sorted(drop)[:6])
        logger.info("%s: dropped %d truncated %s allele(s): %s%s", locus, len(drop), segment,
                    shown, " …" if len(drop) > 6 else "")
    return {a: s for a, s in alleles.items() if a not in drop}


def select_one_allele_per_gene(alleles: dict[str, str], *, logger, locus: str,
                               segment: str) -> dict[str, str]:
    """Keep one representative allele per gene: ``*01`` where it exists, else the lowest-numbered.

    Shrinks the scaffold set ~4x and removes allele-level ambiguity. Deliberately NOT a literal
    ``*01`` filter -- 19 human V genes (IGHV2-70D, IGHV3-25, IGHV3-43D, …) have no ``*01`` record
    and would silently vanish.
    """
    by_gene: dict[str, list[str]] = {}
    for a in alleles:
        by_gene.setdefault(a.split("*")[0], []).append(a)
    keep: dict[str, str] = {}
    no01: list[str] = []
    for al in by_gene.values():
        pick = next((a for a in al if a.endswith("*01")), None)
        if pick is None:
            pick = min(al, key=lambda a: a.rsplit("*", 1)[-1])
            no01.append(pick)
        keep[pick] = alleles[pick]
    if no01:
        logger.info("%s: %d %s gene(s) have no *01 — kept the lowest allele: %s",
                    locus, len(no01), segment, ", ".join(sorted(no01)[:6]))
    logger.info("%s: one-allele-per-gene → %d/%d %s alleles",
                locus, len(keep), len(alleles), segment)
    return keep


def _collect_cdr3_anchors(organism: str, species_dir: str, logger) -> list[dict]:
    """Per-allele junction anchors for every V and J of every locus.

    Independent of the scaffold build: a locus with no scaffolds still gets
    anchors (they are all a bare-junction markup needs). Loci with no IMGT file
    for this species (e.g. rat TR) are skipped.
    """
    fr4_offsets = combinations.load_j_fr4_offsets(organism)
    fwr3_stops = combinations.load_v_fwr3_stops(organism)
    if not fr4_offsets:
        logger.warning("no IgBLAST aux for %s — J anchors fall back to the FR4 motif", organism)
    if not fwr3_stops:
        logger.warning("no IgBLAST ndm for %s — V anchors fall back to the FR3 motif", organism)

    rows: list[dict] = []
    counts: Counter = Counter()

    for locus in LOCI:
        # ---- V: Cys104 from IgBLAST's FWR3 stop, else the conserved FR3 motif.
        try:
            v_alleles = _gapped_alleles(species_dir, locus.group, locus.v)
        except FileNotFoundError:
            logger.info("%s: no IMGT V file (%s) for %s — no anchors",
                        locus.name, locus.v, species_dir)
            v_alleles = []
        for allele, func, gapped in v_alleles:
            ungapped = gapped.replace(".", "")
            anchor, source = _v_anchor(ungapped, fwr3_stops.get(allele))
            status = "ok" if anchor >= 0 else "no_anchor"
            tail = ungapped[anchor:] if anchor >= 0 else ""
            whole = len(tail) // 3 * 3
            rows.append({
                "locus": locus.name, "segment": "V", "allele": allele,
                "functionality": func, "anchor_nt": anchor,
                "partial_nt": len(tail) - whole,          # V ends mid-codon
                "templated_aa": translate(tail[:whole], 0),
                "germline_nt": tail, "status": status, "source": source,
            })
            counts[f"V:{status}"] += 1

        # ---- J: [FW]118 from the aux `cdr3_stop`, else the FR4 motif.
        try:
            j_alleles = _gapped_alleles(species_dir, locus.group, locus.j)
        except FileNotFoundError:
            logger.info("%s: no IMGT J file (%s) for %s — no anchors",
                        locus.name, locus.j, species_dir)
            j_alleles = []
        for allele, func, gapped in j_alleles:
            seq = gapped.replace(".", "")
            if allele in fr4_offsets:
                anchor, source = fr4_offsets[allele][0] + 1, "aux"
            else:
                anchor = _j_anchor_from_motif(seq)
                source = "motif" if anchor >= 0 else "no_anchor"
            if anchor < 0 or anchor + 3 > len(seq):
                anchor, source = -1, "no_anchor"
            # The anchor pins the [FW]118 codon, so every junction codon sits at
            # anchor - 3k and the J's 5' partial run is exactly `anchor % 3`. Do NOT
            # use the aux frame column: for TRAJ31*01 it disagrees with its own
            # cdr3_stop ((anchor - frame) % 3 == 1) and yields a garbage translation.
            frame = anchor % 3 if anchor >= 0 else 0
            status = "ok" if anchor >= 0 else "no_anchor"
            slice_nt = seq[: anchor + 3] if anchor >= 0 else ""
            rows.append({
                "locus": locus.name, "segment": "J", "allele": allele,
                "functionality": func, "anchor_nt": anchor,
                "partial_nt": frame,                       # J starts mid-codon
                "templated_aa": translate(slice_nt[frame:], 0) if slice_nt else "",
                "germline_nt": slice_nt, "status": status, "source": source,
            })
            counts[f"J:{status}"] += 1

    # `truncated`: the germline DOES reach into the junction, but stops at least a codon short of the
    # longest allele of its own gene. `docs/reference_build.rst` has documented this status from the
    # start; it was never implemented, so a truncated allele was indistinguishable from a complete one
    # and its scaffolds silently emitted short junctions (TRAV20*03 templates 3 nt where its siblings
    # template 13 -> every clonotype on it came out 3 aa short). `_process_locus` drops these from the
    # SCAFFOLD set; they stay in this file because VDJdb cites them and bare-junction markup needs
    # their anchors.
    by_gene: dict[tuple[str, str], list[dict]] = {}
    for r in rows:
        if r["status"] == "ok":
            by_gene.setdefault((r["segment"], r["allele"].split("*")[0]), []).append(r)
    for group in by_gene.values():
        best = max(len(r["germline_nt"]) for r in group)
        for r in group:
            if best - len(r["germline_nt"]) >= _TAIL_SLACK_NT:
                r["status"] = "truncated"
                counts[f"{r['segment']}:ok"] -= 1
                counts[f"{r['segment']}:truncated"] += 1

    bad = {k: v for k, v in counts.items() if not k.endswith(":ok")}
    if bad:
        logger.warning("cdr3 anchors: %d allele(s) not clean %s", sum(bad.values()), dict(bad))
    logger.info("cdr3 anchors: %d alleles (%s)", len(rows), dict(counts))
    return rows


def _locus_manifest(nt_all: list[dict], d_germ: list[tuple[str, str, str]]) -> list[dict]:
    """Per-locus reference coverage for EVERY defined locus (0 where the build produced nothing).

    A locus with no scaffolds cannot be annotated for this organism -- reads from it fall through
    silently. A locus with D germlines but no scaffolds ships dead weight: runtime D lookup is keyed
    by a hit scaffold's locus, and that scaffold does not exist. Both are invisible today; this turns
    each into a manifest row and (in ``build_species``) a warning. ``nt_all`` mixes V-J rows (``c_call``
    empty) and J+C rows (``c_call`` set); count them apart."""
    vj = Counter(r["locus"] for r in nt_all if not r.get("c_call"))
    jc = Counter(r["locus"] for r in nt_all if r.get("c_call"))
    dc = Counter(loc for loc, _, _ in d_germ)
    rows = []
    for loc in LOCI:
        n_scaf = vj[loc.name] + jc[loc.name]
        rows.append({
            "locus": loc.name, "group": loc.group,
            "n_vj_scaffolds": vj[loc.name], "n_jc_scaffolds": jc[loc.name],
            "n_d_germlines": dc[loc.name],
            "unreachable_d_germlines": dc[loc.name] if n_scaf == 0 else 0,
            "status": "ok" if n_scaf else "EMPTY",
        })
    return rows


def build_species(organism: str, *, one_allele_per_gene: bool = False) -> Path:
    """Build the reference DB for one organism. Returns the output directory.

    ``one_allele_per_gene`` builds scaffolds from a single representative allele per gene (``*01``
    where it exists, else the lowest-numbered) -- roughly a 4x smaller reference with no
    allele-level ambiguity. Off by default.
    """
    if organism not in IMGT_SPECIES_DIR:
        raise ValueError(f"Unknown organism {organism!r}; one of {list(IMGT_SPECIES_DIR)}")
    species_dir = IMGT_SPECIES_DIR[organism]
    out_dir = vdj_dir(organism)
    logger = _setup_logger(out_dir)
    t0 = time.perf_counter()
    logger.info("=== arda reference build: %s (%s)%s ===", organism, species_dir,
                " [one allele per gene]" if one_allele_per_gene else "")

    imgt.download_reference()
    j_frames = combinations.load_j_frames(organism)
    fwr3_stops = combinations.load_v_fwr3_stops(organism)
    fr4_offsets = combinations.load_j_fr4_offsets(organism)

    nt_all, aa_all, combo_all, fa_nt, fa_aa = [], [], [], [], []
    for locus in LOCI:
        try:
            nt, aa, combo, fnt, faa = _process_locus(
                organism, species_dir, locus, j_frames, logger,
                fwr3_stops=fwr3_stops, fr4_offsets=fr4_offsets,
                one_allele_per_gene=one_allele_per_gene)
        except Exception as exc:  # noqa: BLE001 — one bad locus must not kill the species
            logger.warning("%s: failed (%s) — skipped", locus.name, exc)
            continue
        nt_all += nt
        aa_all += aa
        combo_all += combo
        fa_nt += fnt
        fa_aa += faa

    # A V-J scaffold is V-J all the way to its 3' end and carries no constant region.
    lengths = {i: len(s) for i, s in fa_nt}
    for r in nt_all:
        r["c_call"] = ""
        r["vj_end"] = lengths[r["scaffold_id"]]

    # `J + C` scaffolds, appended -- deliberately NOT routed through IgBLAST, which cannot annotate a
    # V-less sequence; their missing region coordinates would trip the completeness gate in
    # `_process_locus` and every one of them would be dropped, silently. Their geometry is known by
    # construction: J occupies [1, j_len], C occupies [j_len+1, len]. Region coords are -1, which is
    # what `transfer_regions` already emits for a region the query does not reach.
    #
    # FR4 is the exception, and it is not optional. It lies wholly inside J, so a `J + C` scaffold
    # contains all of it, and IgBLAST's own aux file says exactly where (`load_j_fr4_offsets`).
    # Leaving it at -1 cost every J->C read the only markup it can carry -- 19.2 % of mapped reads on
    # SRR5233639 -- and with it `fwr4_aa`, the frame check that separates a real J->C read from a
    # chance hit inside the constant region. Every other region needs the V's conserved Cys104.
    n_vj = len(fa_nt)
    jc = constant.build_jc_scaffolds(organism, species_dir, log=logger)
    for s in jc:
        fa_nt.append((s.scaffold_id, s.sequence))
        row = {"scaffold_id": s.scaffold_id, "locus": s.locus,
               "v_call": "", "j_call": s.j_call, "c_call": s.c_call, "productive": "",
               "v_sequence_end": "", "j_sequence_start": 1,
               "vj_end": s.j_len, "junction": "", "junction_aa": ""}
        for r in REGION_NAMES:
            row[f"{r}_start"], row[f"{r}_end"], row[r] = -1, -1, ""
        if s.fwr4_start > 0:
            row["fwr4_start"], row["fwr4_end"] = s.fwr4_start, s.fwr4_end
            row["fwr4"] = s.sequence[s.fwr4_start - 1 : s.fwr4_end]
        nt_all.append(row)
    logger.info("constant region: %d J+C scaffolds across %d loci (+%.1f%% over %d V-J scaffolds)",
                len(jc), len({s.locus for s in jc}), 100.0 * len(jc) / max(n_vj, 1), n_vj)

    # Write artifacts.
    (out_dir / "alleles.fasta").write_text("".join(f">{i}\n{s}\n" for i, s in fa_nt))
    (out_dir / "alleles.aa.fasta").write_text("".join(f">{i}\n{s}\n" for i, s in fa_aa))
    pl.DataFrame(nt_all, infer_schema_length=None).write_csv(out_dir / "markup.tsv", separator="\t")
    pl.DataFrame(aa_all).write_csv(out_dir / "markup.aa.tsv", separator="\t")
    pl.DataFrame(combo_all).write_csv(out_dir / "combinations.tsv", separator="\t")

    # D germlines for runtime D-segment mapping (VDJ loci only).
    d_germ = _collect_d_germlines(species_dir, logger)
    (out_dir / "d_germlines.fasta").write_text(
        "".join(f">{loc}|{al}\n{s}\n" for loc, al, s in d_germ))
    logger.info("D germlines: %d alleles across %d loci",
                len(d_germ), len({loc for loc, _, _ in d_germ}))

    # Per-allele junction anchors for bare (junction_aa, V, J) markup — see `arda.cdr3fix`.
    anchors = _collect_cdr3_anchors(organism, species_dir, logger)
    pl.DataFrame(anchors, schema={c: pl.Utf8 if c not in ("anchor_nt", "partial_nt") else pl.Int64
                                  for c in _ANCHOR_COLUMNS}).write_csv(
        out_dir / "cdr3_anchors.tsv", separator="\t")

    # Locus coverage manifest — makes a silently-absent locus visible (build.py used to just
    # `continue` past a locus with no IMGT V/J, so rat/rabbit/rhesus shipped no TCR reference and
    # nothing said so). A committed, checkable artifact plus an end-of-build warning.
    manifest = _locus_manifest(nt_all, d_germ)
    pl.DataFrame(manifest).write_csv(out_dir / "loci_manifest.tsv", separator="\t")
    empty = [m["locus"] for m in manifest if m["status"] == "EMPTY"]
    if empty:
        logger.warning("%s: %d/%d loci have NO reference scaffolds — reads from %s will not be "
                       "annotated (see loci_manifest.tsv)", organism, len(empty), len(LOCI),
                       ", ".join(empty))
    for m in manifest:
        if m["unreachable_d_germlines"]:
            logger.warning("%s: %s ships %d D germlines but no scaffolds — unreachable, dead weight",
                           organism, m["locus"], m["unreachable_d_germlines"])

    dt = time.perf_counter() - t0
    logger.info("TOTAL: %d scaffolds across %d loci in %.1fs -> %s",
                len(fa_nt), len({r['locus'] for r in combo_all}), dt, out_dir)
    return out_dir


def build(organism: str = "all", *, one_allele_per_gene: bool = False) -> None:
    """Build one organism or ``"all"`` supported organisms."""
    organisms = SUPPORTED_ORGANISMS if organism == "all" else (organism,)
    for org in organisms:
        build_species(org, one_allele_per_gene=one_allele_per_gene)
