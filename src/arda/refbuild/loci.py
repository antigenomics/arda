"""Locus / chain definitions and species name mappings.

Antigen-receptor loci split into VJ (no D segment) and VDJ (with D). Each locus
maps to the IMGT gene-type FASTA file names and to the IgBLAST ``-ig_seqtype``.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "Locus",
    "LOCI",
    "loci_for",
    "VJ_LOCI",
    "VDJ_LOCI",
    "IMGT_SPECIES_DIR",
    "VDJDB_SPECIES",
    "RECEPTOR_GROUP",
]

# IgBLAST organism name -> IMGT V-QUEST reference directory name.
IMGT_SPECIES_DIR = {
    "human": "Homo_sapiens",
    "mouse": "Mus_musculus",
    "rat": "Rattus_norvegicus",
    "rabbit": "Oryctolagus_cuniculus",
    "rhesus_monkey": "Macaca_mulatta",
}

# VDJdb `species` column -> arda organism (lowercased key; VDJdb writes "HomoSapiens").
VDJDB_SPECIES = {
    "homosapiens": "human",
    "musmusculus": "mouse",
    "macacamulatta": "rhesus_monkey",
    "rattusnorvegicus": "rat",
}

# IMGT splits files under IG/ and TR/ subfolders.
RECEPTOR_GROUP = {"IG": "IG", "TR": "TR"}


@dataclass(frozen=True)
class Locus:
    """A single antigen-receptor locus.

    Attributes:
        name: Locus symbol (e.g. ``"TRB"``, ``"IGH"``).
        group: IMGT receptor group, ``"IG"`` or ``"TR"``.
        v: V gene-type file stem (e.g. ``"TRBV"``).
        j: J gene-type file stem.
        d: D gene-type file stem, or ``None`` for VJ loci.
        ig_seqtype: Value for IgBLAST ``-ig_seqtype`` (``"Ig"`` or ``"TCR"``).
        v_shared: Optional ``(gene_stem, name_substring)`` — also pull V alleles from another gene
            file whose name contains ``name_substring``; an **empty** substring takes the whole
            stem. TRA and TRD share V genes (the loci are interleaved on chr14): IMGT files the
            shared ones under ``TRAV`` with a ``.../DV...`` designation (e.g. ``TRAV14/DV4``).
            Without this, a δ rearrangement on such a V gene has no TRD scaffold to match and is
            miscalled TRA (the locus is set by J/D/C, never the V).

            TRDV-named genes are distinct from the dual-designation TRAV/DV genes. IMGT
            documents TRDV1 rearranged to TRAJ and expressed with TRAC; its V name alone does
            not make that a TRD chain. The optional TRA cross-product includes the TRDV stem.
    """

    name: str
    group: str
    v: str
    j: str
    d: str | None
    ig_seqtype: str
    v_shared: tuple[str, str] | None = None

    @property
    def has_d(self) -> bool:
        return self.d is not None

    @property
    def gene_files(self) -> tuple[str, ...]:
        stems = (self.v, self.d, self.j) if self.has_d else (self.v, self.j)
        return tuple(s for s in stems if s)


LOCI: tuple[Locus, ...] = (
    # VJ loci
    # Keep the shipped TRAV-based TRA vocabulary until the broader cross-product is promoted.
    # This is a reference-scope policy, not a claim that TRDV–TRAJ is impossible: IMGT documents
    # TRDV1–TRAJ/TRAC (https://www.imgt.org/IMGTrepertoire/Regulation/cDNA/human/TRD/TRDV/Hu_TRDVcdna.html).
    # Dual-designation TRAV/DV genes are already covered in both TRA and TRD; J/D/C define locus.
    Locus("TRA", "TR", "TRAV", "TRAJ", None, "TCR"),
    Locus("TRG", "TR", "TRGV", "TRGJ", None, "TCR"),
    Locus("IGK", "IG", "IGKV", "IGKJ", None, "Ig"),
    Locus("IGL", "IG", "IGLV", "IGLJ", None, "Ig"),
    # VDJ loci
    Locus("IGH", "IG", "IGHV", "IGHJ", "IGHD", "Ig"),
    Locus("TRB", "TR", "TRBV", "TRBJ", "TRBD", "TCR"),
    Locus("TRD", "TR", "TRDV", "TRDJ", "TRDD", "TCR", v_shared=("TRAV", "/DV")),
)

VJ_LOCI = tuple(loc for loc in LOCI if not loc.has_d)
VDJ_LOCI = tuple(loc for loc in LOCI if loc.has_d)


def loci_for(*, allow_chimeras: bool = False) -> tuple[Locus, ...]:
    """:data:`LOCI`, optionally with the ``TRDV × TRAJ`` scaffolds the default build refuses.

    The legacy flag name does not establish that a pairing is chimeric. IMGT explicitly
    documents TRDV1*01–TRAJ54*01 and another TRDV1–TRAJ example. That supports admitting the
    observed pairing; it does not by itself validate every member of the broader cross-product.

    Historical real-read evidence (arda-benchmark ``results/round18`` §5c): on 48,030 TRA
    amplicon reads, IgBLAST calls ``TRDV1`` + a ``TRAJ`` on **530 reads (1.10 %)** and MiXCR
    independently agrees. All have ``v_score ≥ 70`` (median 93.8) and an IgBLAST junction.
    Arda calls their J identically while emitting no V/junction: **83 %** of its remaining
    V-gene gap (.9867 vs MiXCR's .9973; .9978 excluding these reads).

    ``False`` preserves the shipped reference scope. Ordinary TRG/TRD annotation already exists;
    this is not a gamma-delta mode. The opt-in build requires complete expected scaffold coverage
    and appends targets without renumbering or reframing existing ones. An earlier build added
    only 7 of the expected 483 scaffolds because of missing markup; assert coverage, not the flag.
    """
    if not allow_chimeras:
        return LOCI
    return tuple(
        Locus(loc.name, loc.group, loc.v, loc.j, loc.d, loc.ig_seqtype,
              v_shared=("TRDV", ""))            # empty needle == the whole TRDV stem
        if loc.name == "TRA" else loc
        for loc in LOCI
    )
