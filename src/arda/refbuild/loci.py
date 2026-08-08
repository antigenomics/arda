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


def loci_for(*, allow_chimeras: bool = False) -> tuple[Locus, ...]:
    """:data:`LOCI`, optionally with the ``TRDV × TRAJ`` scaffolds the default build refuses.

    The refusal above is a *biological* claim: TRDV1/2/3 are dedicated δ V genes, so a
    ``TRDV × TRAJ`` scaffold is a chimera. The claim is testable and the test disagrees with it on
    the only external evidence available. Measured on 48,030 TRA amplicon reads against an IgBLAST
    truth (arda-benchmark ``results/round18`` §5c): **530 reads — 1.10 % of the library — are
    called ``TRDV1`` + a ``TRAJ`` by IgBLAST, and MiXCR independently agrees**; all sit at IgBLAST
    ``v_score ≥ 70`` (median 93.8), all carry an IgBLAST junction, and **arda calls their J
    identically to both tools** while emitting no ``v_call`` and no junction at all. They are
    **83 % of arda's entire remaining ``v_gene`` gap** on that library (.9867 vs MiXCR's .9973;
    .9978 excluding them).

    So the two readings are:

    * the pairing is real — arda silently drops 1.1 % of a TRA repertoire, or
    * it is a chimera — arda is right and IgBLAST and MiXCR both report it because they call V and
      J independently, with no combination constraint at all.

    This flag exists because that is a domain judgement, not a code decision, and the default must
    not quietly encode either answer as if it were settled. ``False`` keeps the shipped biology.

    ⛔ It is **not** free to turn on. An earlier attempt at the same scaffolds measured **+7**
    scaffolds, not the ~483 the cross-product predicts, because 62 of 69 dropped for incomplete
    markup — so a reference built with this on may not contain what you expect. Assert the
    scaffold count, never the flag (CLAUDE.md: a reference swap can silently be a no-op).
    """
    if not allow_chimeras:
        return LOCI
    return tuple(
        Locus(loc.name, loc.group, loc.v, loc.j, loc.d, loc.ig_seqtype,
              v_shared=("TRDV", ""))            # empty needle == the whole TRDV stem
        if loc.name == "TRA" else loc
        for loc in LOCI
    )
