"""Locus / chain definitions and species name mappings.

Antigen-receptor loci split into VJ (no D segment) and VDJ (with D). Each locus
maps to the IMGT gene-type FASTA file names and to the IgBLAST ``-ig_seqtype``.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "Locus",
    "LOCI",
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

            ⛔ The sharing runs **both ways**, and only one direction was wired. TRDV1/2/3 are
            dedicated δ V genes but lie *inside* the TRA locus, between the TRAV genes and the TRAJ
            cluster, so an α rearrangement can join one to a TRAJ. Without a ``TRDV × TRAJ``
            scaffold such a read gets its J called and **no V at all**, hence no junction.
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
    # ⛔ TRA does NOT pull the TRDV stem, and a previous attempt to make it do so was biologically
    # wrong. The sharing is not symmetric:
    #
    #   * **TRDV1/2/3 are dedicated delta V genes.** They rearrange to TRDJ. A `TRDV1 + TRAJ`
    #     scaffold is not a rearrangement that happens, so building one invites reads onto a
    #     chimera the biology does not contain.
    #   * **TRAV/DV genes pair with EITHER**, and which J they took is what defines the locus.
    #     Those are already covered from both sides: `TRAV/DV x TRAJ` comes free with the TRA
    #     build because IMGT files them under TRAV, and `TRAV/DV x TRDJ` is exactly what TRD's
    #     `v_shared=("TRAV", "/DV")` below exists to add.
    #
    # So the asymmetry in this table is the correct encoding of the biology, not an oversight.
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
