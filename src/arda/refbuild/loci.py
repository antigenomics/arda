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
    """

    name: str
    group: str
    v: str
    j: str
    d: str | None
    ig_seqtype: str

    @property
    def has_d(self) -> bool:
        return self.d is not None

    @property
    def gene_files(self) -> tuple[str, ...]:
        stems = (self.v, self.d, self.j) if self.has_d else (self.v, self.j)
        return tuple(s for s in stems if s)


LOCI: tuple[Locus, ...] = (
    # VJ loci
    Locus("TRA", "TR", "TRAV", "TRAJ", None, "TCR"),
    Locus("TRG", "TR", "TRGV", "TRGJ", None, "TCR"),
    Locus("IGK", "IG", "IGKV", "IGKJ", None, "Ig"),
    Locus("IGL", "IG", "IGLV", "IGLJ", None, "Ig"),
    # VDJ loci
    Locus("IGH", "IG", "IGHV", "IGHJ", "IGHD", "Ig"),
    Locus("TRB", "TR", "TRBV", "TRBJ", "TRBD", "TCR"),
    Locus("TRD", "TR", "TRDV", "TRDJ", "TRDD", "TCR"),
)

VJ_LOCI = tuple(l for l in LOCI if not l.has_d)
VDJ_LOCI = tuple(l for l in LOCI if l.has_d)
