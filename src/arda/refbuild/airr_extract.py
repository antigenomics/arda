"""Annotate V-J scaffolds with IgBLAST and extract AIRR region markup.

For each locus we build germline BLAST databases from the ungapped IMGT files,
run ``igblastn -outfmt 19`` on the assembled scaffolds, and read the AIRR TSV
with polars. The scaffold sequences contain verbatim germline V and J, so
IgBLAST finds exact matches and reports precise FR/CDR coordinates.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from ..paths import data_dir, bin_dir
from .. import igblast
from .imgt import ungap_gene
from .loci import Locus

__all__ = ["REGION_NAMES", "AIRR_MARKUP_COLUMNS", "build_germline_dbs", "annotate_scaffolds"]

REGION_NAMES = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")

# Columns we pull from the AIRR output (1-based, closed coords for *_start/_end).
AIRR_MARKUP_COLUMNS = (
    ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call",
     "v_sequence_start", "v_sequence_end", "j_sequence_start", "j_sequence_end",
     "junction", "junction_aa"]
    + [f"{r}" for r in REGION_NAMES]
    + [f"{r}_aa" for r in REGION_NAMES]
    + [f"{r}_start" for r in REGION_NAMES]
    + [f"{r}_end" for r in REGION_NAMES]
)


def _blastdb_dir(species_dir: str) -> Path:
    d = data_dir() / "blastdb" / species_dir
    d.mkdir(parents=True, exist_ok=True)
    return d


def _dummy_d_db(species_dir: str) -> Path:
    """A placeholder D database for VJ loci (IgBLAST requires -germline_db_D)."""
    prefix = _blastdb_dir(species_dir) / "_dummyD"
    if not Path(str(prefix) + ".nin").exists():
        fa = _blastdb_dir(species_dir) / "_dummyD.fasta"
        fa.write_text(">dummyD\nGGGGGGGGGGGGGGGGGGGG\n")
        igblast.makeblastdb(fa, prefix, dbtype="nucl")
    return prefix


def build_germline_dbs(species_dir: str, locus: Locus) -> dict[str, Path]:
    """Ungap each gene file and build a germline BLAST DB; return {role: prefix}."""
    out: dict[str, Path] = {}
    roles = {"V": locus.v, "J": locus.j}
    if locus.has_d:
        roles["D"] = locus.d  # type: ignore[assignment]
    for role, stem in roles.items():
        ungapped = ungap_gene(species_dir, locus.group, stem)
        prefix = _blastdb_dir(species_dir) / stem
        igblast.makeblastdb(ungapped, prefix, dbtype="nucl")
        out[role] = prefix
    return out


def annotate_scaffolds(
    scaffold_fasta: Path,
    organism: str,
    species_dir: str,
    locus: Locus,
    *,
    num_threads: int = 1,
) -> pl.DataFrame:
    """Run IgBLAST on a scaffold FASTA and return the markup columns as polars."""
    dbs = build_germline_dbs(species_dir, locus)
    aux = bin_dir() / "optional_file" / f"{organism}_gl.aux"
    out_tsv = scaffold_fasta.with_suffix(".airr.tsv")
    igblast.igblastn_airr(
        scaffold_fasta,
        out_tsv,
        organism=organism,
        germline_db_v=dbs["V"],
        germline_db_j=dbs["J"],
        germline_db_d=dbs.get("D") or _dummy_d_db(species_dir),
        auxiliary_data=aux if aux.exists() else None,
        ig_seqtype=locus.ig_seqtype,
        num_threads=num_threads,
    )
    df = pl.read_csv(out_tsv, separator="\t", infer_schema_length=0)
    keep = [c for c in AIRR_MARKUP_COLUMNS if c in df.columns]
    return df.select(keep)
