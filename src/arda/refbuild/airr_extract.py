"""Annotate V-J scaffolds with IgBLAST and extract AIRR region markup.

For each locus we build germline BLAST databases from the ungapped IMGT files,
run ``igblastn -outfmt 19`` on the assembled scaffolds, and read the AIRR TSV
with polars. The scaffold sequences contain verbatim germline V and J, so
IgBLAST finds exact matches and reports precise FR/CDR coordinates.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from ..paths import data_dir
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


def _records(path: Path) -> list[tuple[str, str]]:
    """``[(seq_id, ">header\nseq...")]`` from a FASTA, headers kept verbatim."""
    text = path.read_text()
    return [(b.split("\n", 1)[0].split()[0], ">" + b) for b in text.split(">") if b.strip()]


def _shared_v_fasta(species_dir: str, locus: Locus) -> tuple[Path, str]:
    """Ungapped V FASTA covering ``locus.v`` **and** ``locus.v_shared``; ``(path, db_name)``.

    ``build._process_locus`` merges the shared stem's dual-use genes into the allele set the
    scaffolds are built from, so a V database built from ``locus.v`` alone does not contain the
    germline half of the scaffolds it is about to mark up. IgBLAST then calls the nearest
    same-stem gene instead and every region coordinate collapses, which is silent: the scaffold
    is simply dropped for incomplete markup. Measured on the chimera-enabled TRA locus
    (benchmark round 27), complete markup went **7/483 -> 49/483** on this one line.

    ⚠ Dedupe by seq id, first wins. IMGT files the dual-use ``TRAV*/DV*`` genes under BOTH the
    TRAV and TRDV stems and ``makeblastdb`` dies hard on
    ``Duplicate seq_ids are found: LCL|TRAV14/DV4*01``.
    """
    stem, needle = locus.v_shared          # type: ignore[misc]
    base = ungap_gene(species_dir, locus.group, locus.v)
    shared = ungap_gene(species_dir, locus.group, stem)
    seen: set[str] = set()
    keep: list[str] = []
    for path, want in ((base, None), (shared, needle)):
        for sid, rec in _records(path):
            if (want is not None and want not in sid) or sid in seen:
                continue
            seen.add(sid)
            keep.append(rec)
    name = f"{locus.v}_{stem}"
    merged = _blastdb_dir(species_dir) / f"{name}.fasta"
    merged.write_text("".join(keep))
    return merged, name


def build_germline_dbs(species_dir: str, locus: Locus) -> dict[str, Path]:
    """Ungap each gene file and build a germline BLAST DB; return {role: prefix}."""
    out: dict[str, Path] = {}
    roles = {"V": locus.v, "J": locus.j}
    if locus.has_d:
        roles["D"] = locus.d  # type: ignore[assignment]
    for role, stem in roles.items():
        ungapped, name = ungap_gene(species_dir, locus.group, stem), stem
        if role == "V" and locus.v_shared:
            ungapped, name = _shared_v_fasta(species_dir, locus)
        prefix = _blastdb_dir(species_dir) / name
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
    # Resolved through `igblast.auxiliary_data`, which raises when it is missing. It used to be
    # looked up under `bin_dir()` and passed as None when absent -- the same directory in a
    # source checkout, a different one on every auto-fetched install, and the silent result is
    # an IgBLAST run with NO junction on any read.
    aux = igblast.auxiliary_data(organism)
    out_tsv = scaffold_fasta.with_suffix(".airr.tsv")
    igblast.igblastn_airr(
        scaffold_fasta,
        out_tsv,
        organism=organism,
        germline_db_v=dbs["V"],
        germline_db_j=dbs["J"],
        germline_db_d=dbs.get("D") or _dummy_d_db(species_dir),
        auxiliary_data=aux,
        ig_seqtype=locus.ig_seqtype,
        num_threads=num_threads,
    )
    df = pl.read_csv(out_tsv, separator="\t", infer_schema_length=0)
    keep = [c for c in AIRR_MARKUP_COLUMNS if c in df.columns]
    return df.select(keep)
