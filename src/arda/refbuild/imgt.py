"""IMGT/V-QUEST germline reference download, parsing, and ungapping.

The IMGT V-QUEST reference directory ships gapped germline FASTAs laid out as
``<Species>/<IG|TR>/<GENE>.fasta`` (e.g. ``Homo_sapiens/IG/IGHV.fasta``).
Sequences carry IMGT-numbering gap dots; IgBLAST's ``edit_imgt_file.pl`` ungaps
them and rewrites headers to bare allele names (what ``makeblastdb`` wants).

This module:

* downloads & extracts the reference zip into ``data/imgt`` (gitignored),
* parses the original gapped FASTA headers for per-allele *functionality*,
* ungaps a gene file via ``edit_imgt_file.pl`` into ``data/imgt/ungapped``.
"""

from __future__ import annotations

import zipfile
from dataclasses import dataclass
from pathlib import Path

import requests

from ..paths import data_dir
from .. import igblast

__all__ = [
    "REFERENCE_URL",
    "ImgtAllele",
    "download_reference",
    "gene_fasta_path",
    "parse_functionality",
    "ungap_gene",
    "read_fasta",
]

REFERENCE_URL = (
    "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory.zip"
)

# Functionality codes we treat as usable germline (functional + ORF). IMGT may
# wrap codes in parentheses/brackets for inference, e.g. "(F)" or "[F]".
_FUNCTIONAL = {"F", "ORF"}


@dataclass(frozen=True)
class ImgtAllele:
    """A germline allele parsed from an IMGT FASTA header + sequence."""

    allele: str          # e.g. "IGHV4-59*01"
    functionality: str   # normalized: "F", "ORF", "P", ...
    sequence: str        # as stored in the source file (gapped or ungapped)

    @property
    def is_functional(self) -> bool:
        return self.functionality in _FUNCTIONAL


def _imgt_root() -> Path:
    return data_dir() / "imgt"


def reference_dir() -> Path:
    """Directory the reference zip is extracted into."""
    return _imgt_root() / "reference"


def download_reference(*, force: bool = False) -> Path:
    """Download and extract the IMGT V-QUEST reference directory.

    Returns the extraction root (containing the per-species directories).
    Idempotent unless ``force``.
    """
    root = reference_dir()
    if root.is_dir() and any(root.iterdir()) and not force:
        return root
    root.mkdir(parents=True, exist_ok=True)
    zip_path = _imgt_root() / "IMGT_V-QUEST_reference_directory.zip"
    if force or not zip_path.exists():
        with requests.get(REFERENCE_URL, stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with open(zip_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(root)
    return root


def gene_fasta_path(species_dir: str, group: str, gene_stem: str) -> Path:
    """Path to a gene-type FASTA, e.g. ``Homo_sapiens/IG/IGHV.fasta``.

    Handles the occasional top-level wrapper directory inside the zip.
    """
    root = reference_dir()
    rel = Path(species_dir) / group / f"{gene_stem}.fasta"
    direct = root / rel
    if direct.exists():
        return direct
    # The zip sometimes nests everything under a single wrapper folder.
    matches = list(root.glob(f"*/{rel}")) + list(root.glob(str(rel)))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"IMGT gene file not found: {rel} under {root}")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA file into ``(header, sequence)`` pairs (sequence joined)."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            elif line:
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def parse_functionality(path: Path) -> dict[str, str]:
    """Map allele name -> normalized functionality from gapped IMGT headers.

    IMGT header: ``accession|allele|species|functionality|region|...``. The
    functionality field may be wrapped, e.g. ``(F)`` / ``[F]`` for inferred.
    """
    out: dict[str, str] = {}
    for header, _seq in read_fasta(path):
        fields = header.split("|")
        if len(fields) < 4:
            continue
        allele = fields[1].strip()
        func = fields[3].strip().strip("()[]")
        # Take the first token if multiple (e.g. "ORF/F").
        func = func.split("/")[0].split()[0] if func else ""
        out[allele] = func
    return out


def ungap_gene(species_dir: str, group: str, gene_stem: str) -> Path:
    """Ungap a gene file with ``edit_imgt_file.pl``; return the ungapped path."""
    src = gene_fasta_path(species_dir, group, gene_stem)
    out_dir = _imgt_root() / "ungapped" / species_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{gene_stem}.fasta"
    igblast.edit_imgt_file(src, out)
    return out


def load_functional_alleles(
    species_dir: str, group: str, gene_stem: str
) -> dict[str, str]:
    """Return ``{allele: ungapped_sequence}`` for functional/ORF alleles only."""
    functionality = parse_functionality(gene_fasta_path(species_dir, group, gene_stem))
    ungapped = ungap_gene(species_dir, group, gene_stem)
    out: dict[str, str] = {}
    for header, seq in read_fasta(ungapped):
        # edit_imgt_file rewrites the header to the bare allele name.
        allele = header.split("|")[0].strip().split()[0]
        func = functionality.get(allele, "")
        if func in _FUNCTIONAL and seq:
            out[allele] = seq.upper().replace(".", "")
    return out
