"""Enumerate in-frame V-J reference scaffolds.

For markup transfer the FR/CDR region *coordinates* are fully determined by the
V gene (FR1-3, CDR1-2, CDR3 start at the conserved Cys104) and the J gene (CDR3
end, FR4). The D segment lies inside the hypervariable CDR3 — query-specific at
runtime — so we enumerate **V×J** scaffolds for every locus and, for VDJ loci,
insert a short frame-neutral N spacer where D would sit so IgBLAST still
annotates a plausible CDR3 + FR4.

Each scaffold is ``V + N*pad + J`` where ``pad`` keeps the J coding frame aligned
to V's reading frame (``jframe`` from the IgBLAST aux file). Byte-identical
scaffolds are **deduplicated**: one DB entry, with all contributing (V,J) allele
pairs recorded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..paths import bin_dir
from .loci import Locus
from .translate import detect_coding_frame

__all__ = [
    "Scaffold",
    "DEFAULT_D_SPACER_NT",
    "load_j_frames",
    "load_j_fr4_offsets",
    "load_v_fwr3_stops",
    "build_locus_scaffolds",
]

# Frame-neutral N run (3 codons) standing in for a D segment in VDJ loci, so the
# CDR3 is a realistic length. Must be a multiple of 3 to preserve frame.
DEFAULT_D_SPACER_NT = 9


@dataclass
class Scaffold:
    """A deduplicated V-J reference scaffold.

    Fields: ``scaffold_id`` (stable ``"{locus}_{index}"``), ``locus``,
    ``sequence`` (assembled ``V + N*pad + J``), ``v_calls`` / ``j_calls`` (all
    alleles producing this scaffold), and ``n_pad`` (N nucleotides between V and J).
    """

    scaffold_id: str
    locus: str
    sequence: str
    v_calls: list[str] = field(default_factory=list)
    j_calls: list[str] = field(default_factory=list)
    n_pad: int = 0


def load_j_frames(organism: str) -> dict[str, int]:
    """Parse ``bin/optional_file/<organism>_gl.aux`` -> {J allele: frame}.

    Frame is the 0-based "first coding frame start position" (column 2).
    """
    aux = bin_dir() / "optional_file" / f"{organism}_gl.aux"
    frames: dict[str, int] = {}
    if not aux.exists():
        return frames
    with open(aux) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    frames[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return frames


def load_v_fwr3_stops(organism: str) -> dict[str, int]:
    """Parse ``bin/internal_data/<organism>/<organism>.ndm.imgt`` -> {V allele: FWR3 stop}.

    IgBLAST's own IMGT annotation of its V germlines. Column 11 is the 1-based FWR3
    stop, and FR3-IMGT ends at position 104 -- the conserved 2nd-CYS -- so the Cys104
    codon starts at ``fwr3_stop - 3`` (0-based). This is the authoritative V junction
    anchor, and it is the same IgBLAST metadata the rest of the build already trusts.

    It does NOT cover every IMGT allele (IgBLAST ships a subset), hence the motif
    fallback in ``refbuild.build._v_anchor``.
    """
    ndm = bin_dir() / "internal_data" / organism / f"{organism}.ndm.imgt"
    out: dict[str, int] = {}
    if not ndm.exists():
        return out
    with open(ndm) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                out[parts[0]] = int(parts[10])
            except ValueError:
                continue
    return out


def load_j_fr4_offsets(organism: str) -> dict[str, tuple[int, int]]:
    """Parse the same aux file -> ``{J allele: (cdr3_stop, extra_bp)}``, both 0-based nt counts.

    IgBLAST's aux carries five columns: allele, coding-frame start, chain type, **CDR3 stop**, and
    **extra bp beyond the J coding end**. arda only ever read column 2. Columns 4 and 5 pin FR4
    inside the J exactly::

        fwr4 = j_seq[cdr3_stop + 1 : len(j_seq) - extra_bp]

    That is how a ``J + C`` scaffold gets an FR4 at all: ``igblastn`` cannot annotate a V-less
    sequence, so those scaffolds are not routed through it (see ``refbuild.build``).

    Verified against every V-J scaffold arda builds: the string this yields is byte-identical to
    IgBLAST's own ``fwr4`` on all 125 human J alleles where both exist -- including IgBLAST's own
    non-multiple-of-3 cases (``IGHJ6*02`` has ``extra_bp = 0`` and a 34 nt FR4; 726 V-J scaffolds
    already carry one). Reproducing IgBLAST, quirks included, is the requirement here: the two
    scaffold kinds must agree, or a J->C read and a V-J read of the same clone disagree on FR4.

    Pseudogene J entries carry only three columns and are skipped -- they have no FR4 to report.
    """
    aux = bin_dir() / "optional_file" / f"{organism}_gl.aux"
    out: dict[str, tuple[int, int]] = {}
    if not aux.exists():
        return out
    with open(aux) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    out[parts[0]] = (int(parts[3]), int(parts[4]))
                except ValueError:
                    continue
    return out


def build_locus_scaffolds(
    locus: Locus,
    v_alleles: dict[str, str],
    j_alleles: dict[str, str],
    j_frames: dict[str, int],
    *,
    d_spacer: int | None = None,
) -> list[Scaffold]:
    """Build deduplicated V×J scaffolds for one locus.

    Args:
        locus: The locus definition.
        v_alleles: ``{allele: ungapped V sequence}``.
        j_alleles: ``{allele: ungapped J sequence}``.
        j_frames: ``{J allele: 0-based coding frame}`` from the aux file.
        d_spacer: N spacer length for VDJ loci (default ``DEFAULT_D_SPACER_NT``);
            forced to 0 for VJ loci.

    Returns:
        Scaffolds, one per unique assembled sequence.
    """
    spacer = (DEFAULT_D_SPACER_NT if d_spacer is None else d_spacer) if locus.has_d else 0

    # Normalize each V to its coding frame: trim leading nt so the V (and thus
    # the whole scaffold) reads in frame 0. Partial-5' alleles otherwise start
    # mid-codon and frameshift the J/FR4 translation.
    v_inframe = {a: s[detect_coding_frame(s):] for a, s in v_alleles.items()}

    # seq -> (set of v_alleles, set of j_alleles, n_pad)
    seen: dict[str, tuple[set[str], set[str], int]] = {}
    for v_allele, v_seq in v_inframe.items():
        lv = len(v_seq)
        for j_allele, j_seq in j_alleles.items():
            jframe = j_frames.get(j_allele, 0)
            base_pad = (3 - ((lv + jframe) % 3)) % 3
            n_pad = base_pad + spacer
            seq = v_seq + ("N" * n_pad) + j_seq
            if seq in seen:
                vs, js, _ = seen[seq]
                vs.add(v_allele)
                js.add(j_allele)
            else:
                seen[seq] = ({v_allele}, {j_allele}, n_pad)

    scaffolds: list[Scaffold] = []
    for idx, (seq, (vs, js, n_pad)) in enumerate(sorted(seen.items())):
        scaffolds.append(
            Scaffold(
                scaffold_id=f"{locus.name}_{idx}",
                locus=locus.name,
                sequence=seq,
                v_calls=sorted(vs),
                j_calls=sorted(js),
                n_pad=n_pad,
            )
        )
    return scaffolds
