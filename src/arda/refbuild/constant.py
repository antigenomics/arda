"""Constant-region (``J + C``) scaffolds — the reads a ``V + pad + J`` reference cannot represent.

arda's V-J scaffold ends at the J segment's 3' terminus, so a read running from the J across the
splice junction into the constant region has nowhere legal to align. That class is **22 % of all real
receptor fragments** in bulk RNA-seq, and arda maps only about half of them.

The fix is *not* to append C to every V-J scaffold. **A J->C read has no V** — that is the definition
of the class — so C never has to sit downstream of a V. Appending ``C[:150]`` to all 17,244 human V-J
scaffolds multiplies them by the number of distinct C stubs per locus (11 for IGH) and yields 73,360
scaffolds, 37 Mnt: mmseqs' per-query cost rises 1.71x and its peak RSS by 367 MB.

Adding 345 ``J + C`` scaffolds *beside* the 17,244 V-J ones (human) is accuracy-identical and costs
+2.0 % scaffolds and +9 MB RSS. A read spanning V->J->C still maps to its V-J scaffold and soft-clips the C tail
(V-supported recall is already 100 %); a V-less J->C read finally has somewhere legal to end.

These scaffolds never see IgBLAST. It cannot annotate a sequence with no V, and ``build.build_species``
drops any scaffold with a missing region coordinate — so routing them through the annotator would
silently delete every one of them. Their geometry is known by construction instead: ``J`` occupies
``[1, j_len]`` and ``C`` occupies ``[j_len + 1, len]``.

``vj_end`` is the length of the V-J part of a scaffold (``j_len`` here, the full length for a V-J
scaffold). A hit with ``tstart >= vj_end`` lies wholly inside the constant region: real receptor mRNA,
but carrying no V(D)J and therefore no clonotype. Callers drop those.

Source of the C sequences: ``database/c_genes/<organism>.fasta`` — the CH1 exon of each constant-region
gene, lifted from the reference genome assembly by exon coordinates. CH1 is the exon that splices
directly onto J, and it begins **mid-codon**: the codon straddles the J-C splice. So ``J + CH1``
reconstructs the mRNA and translates contiguously, which is the bundle's acceptance test
(``IGHJ4 + IGHG1`` -> ``WGQGTLVTVSS|ASTKGPSVFP``). See ``database/c_genes/README.md``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from ..paths import database_dir
from . import imgt
from .loci import LOCI

__all__ = ["JCScaffold", "C_STUB_NT", "c_genes_path", "isotype_class", "build_jc_scaffolds"]

logger = logging.getLogger(__name__)

# A read cannot see further into C than its own length. 150 nt covers every Illumina read in the
# benchmark set (51-152 bp) with room for the J-side overhang. Longer stubs only inflate the DB:
# IGHG1-4 are ~95 % identical, so their 5' stubs collapse to fewer distinct sequences as it shrinks.
C_STUB_NT = 150


@dataclass(frozen=True)
class JCScaffold:
    scaffold_id: str
    locus: str
    j_call: str      # comma-joined: distinct J alleles sharing this exact sequence
    c_call: str      # comma-joined: distinct C genes sharing this exact stub
    j_len: int       # nt length of the J part == `vj_end`
    sequence: str


def c_genes_path(organism: str) -> Path:
    return database_dir() / "c_genes" / f"{organism}.fasta"


def germline_supplement_path(organism: str, gene_stem: str) -> Path:
    """Germline for a locus the IMGT V-QUEST reference directory omits entirely.

    V-QUEST ships **no TR directory at all** for rat, so arda's rat ``TRAC1``/``TRBC1``/``TRBC2`` CH1
    exons had no J to splice onto and rat got zero ``J + C`` scaffolds.

    Consulted **only when IMGT has no file for that stem** — a supplement that shadowed an IMGT gene
    would be far worse than a missing one. It cannot rescue a V-J scaffold: those need IgBLAST, whose
    ``internal_data`` carries TR annotation for human and mouse alone, so ``_process_locus`` skips rat
    TR regardless of germline. A ``J + C`` scaffold never sees IgBLAST.
    """
    return database_dir() / "germline" / organism / f"{gene_stem}.fasta"


def _load_j_alleles(organism: str, species_dir: str, locus, log=logger) -> dict[str, str]:
    """IMGT first; the shipped supplement only if IMGT has nothing for this stem."""
    try:
        alleles = imgt.load_functional_alleles(species_dir, locus.group, locus.j)
        if alleles:
            return alleles
    except FileNotFoundError:
        pass
    path = germline_supplement_path(organism, locus.j)
    if not path.exists():
        log.info("%s: no IMGT J file (%s) for %s and no supplement — no J+C scaffolds",
                 locus.name, locus.j, species_dir)
        return {}
    alleles = {hdr.split()[0]: seq.upper() for hdr, seq in imgt.read_fasta(path)}
    log.info("%s: J germline from SUPPLEMENT (%d alleles) — IMGT ships no %s for %s",
             locus.name, len(alleles), locus.j, species_dir)
    return alleles


def isotype_class(c_call: str) -> str:
    """``IGHG1`` -> ``IGHG``; ``IGHG1,IGHG3`` -> ``IGHG``; ``IGKC`` -> ``IGKC``; ``IGHG1,IGHM`` -> ``IGHC``.

    **Report the class, never the subclass.** IGHG1-4 are ~95 % identical over CH1, so the top-scoring
    *gene* is a coin-flip between them: on real data the best gene ties 26.7 % of the time, while the
    best *class* is unique on every read. A subclass call from a 100 bp read is noise dressed as data.

    When the calls straddle more than one class the answer is the locus-level constant gene
    (``IGHC``, ``IGLC``, ``TRBC`` ...), not the empty string: "some IGH constant region" is a true and
    useful statement, and an empty ``c_class`` would be indistinguishable from "no C hit at all".
    """
    if not c_call:
        return ""
    classes = set()
    for g in c_call.split(","):
        g = g.split("*")[0]
        # IGHG1..4 -> IGHG, IGHA1/2 -> IGHA; IGHM/IGHD/IGHE and all TR/IGK/IGL genes are their own class
        classes.add(g[:4] if g.startswith(("IGHG", "IGHA")) else g)
    if len(classes) == 1:
        return classes.pop()
    return _locus_of(c_call.split(",")[0]) + "C"      # e.g. IGHG1,IGHM -> IGHC


def _locus_of(c_gene: str) -> str:
    """``IGHG1`` -> ``IGH``, ``IGKC`` -> ``IGK``, ``TRBC2`` -> ``TRB``."""
    return c_gene.split("*")[0][:3]


def build_jc_scaffolds(organism: str, species_dir: str, c_len: int = C_STUB_NT,
                       log=logger) -> list[JCScaffold]:
    """One scaffold per (distinct J allele sequence x distinct C stub) within each locus.

    ``log`` is the per-species build logger, so a supplement fallback is recorded in that species'
    ``build.log`` rather than vanishing into a module logger nobody reads.
    """
    cpath = c_genes_path(organism)
    if not cpath.exists():
        return []

    # locus -> {stub_seq: [c_gene, ...]}   (IGHG1-4 collapse; IGLC1-7 mostly do not)
    by_locus: dict[str, dict[str, list[str]]] = {}
    for header, seq in imgt.read_fasta(cpath):
        gene = header.split()[0]
        stub = seq.upper()[:c_len]
        if len(stub) < 20:
            continue
        by_locus.setdefault(_locus_of(gene), {}).setdefault(stub, []).append(gene)

    out: list[JCScaffold] = []
    for locus in LOCI:
        stubs = by_locus.get(locus.name)
        if not stubs:
            continue
        # A locus without an IMGT J file must not kill the species build (`_process_locus` guards its
        # own per-locus work the same way). Rat has TRAC/TRBC genes but IMGT ships no rat TR directory,
        # so its J comes from the shipped supplement. A J+C scaffold needs no IgBLAST and no V, so it is
        # worth building even for a locus whose V-J scaffolds were skipped for want of TR annotation.
        j_alleles = _load_j_alleles(organism, species_dir, locus, log)
        if not j_alleles:
            continue
        # collapse J alleles that share an identical sequence, exactly as the V-J builder does
        j_by_seq: dict[str, list[str]] = {}
        for allele, s in j_alleles.items():
            j_by_seq.setdefault(s.upper(), []).append(allele)

        idx = 0
        for jseq, jcalls in sorted(j_by_seq.items()):
            for stub, cgenes in sorted(stubs.items()):
                out.append(JCScaffold(
                    scaffold_id=f"{locus.name}_JC_{idx}", locus=locus.name,
                    j_call=",".join(sorted(jcalls)), c_call=",".join(sorted(cgenes)),
                    j_len=len(jseq), sequence=jseq + stub))
                idx += 1
    return out
