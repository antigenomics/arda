"""Gold-standard read annotation with IgBLAST (AIRR ``-outfmt 19``), all loci.

IgBLAST is arda's offline gold standard. This runs it on arbitrary reads the way a
repertoire study does: one combined germline BLAST DB per receptor type (all TR loci
in one V/D/J DB, all IG loci in another), ``igblastn`` once per type the organism has
internal annotation for, then the per-type AIRR TSVs merged keeping the best-scoring
hit per read. The IMGT reference is fetched once if missing; combined DBs are cached
under ``data/blastdb``. Exposed as ``arda igblast`` and used by the benchmark.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import polars as pl

from ..paths import data_dir, bin_dir
from ..annotate import io as seqio
from .. import igblast
from . import imgt
from .airr_extract import _dummy_d_db
from .loci import LOCI, IMGT_SPECIES_DIR

__all__ = ["igblast_reads"]

_IG_SEQTYPE = {"TR": "TCR", "IG": "Ig"}


def _blastdb_dir(species_dir: str) -> Path:
    d = data_dir() / "blastdb" / species_dir
    d.mkdir(parents=True, exist_ok=True)
    return d


def _combined_db(species_dir: str, group: str, role: str, stems: list[str]) -> Path:
    """Concatenate the ungapped gene files for ``stems`` into one germline BLAST DB.

    Cached: reused across runs once ``.nin`` exists. This is the standard IgBLAST
    setup — a single V (or D/J) database spanning every locus of one receptor type.
    """
    prefix = _blastdb_dir(species_dir) / f"_all_{group}_{role}"
    if Path(str(prefix) + ".nin").exists():
        return prefix
    fa = _blastdb_dir(species_dir) / f"_all_{group}_{role}.fasta"
    # Dedup by allele: dual-designation genes (e.g. TRAV14/DV4*01) live in both the
    # TRAV and TRDV files, and makeblastdb -parse_seqids rejects duplicate seq_ids.
    seen: set[str] = set()
    with open(fa, "w") as out:
        for stem in stems:
            ungapped = Path(imgt.ungap_gene(species_dir, group, stem))
            for header, seq in imgt.read_fasta(ungapped):
                allele = header.split("|")[0].strip().split()[0]
                if not seq or allele in seen:
                    continue
                seen.add(allele)
                out.write(f">{allele}\n{seq}\n")
    igblast.makeblastdb(fa, prefix, dbtype="nucl")
    return prefix


def _run_group(query_fa: Path, organism: str, species_dir: str, group: str,
               num_threads: int, out_tsv: Path) -> Path:
    loci = [l for l in LOCI if l.group == group]
    v = _combined_db(species_dir, group, "V", [l.v for l in loci])
    j = _combined_db(species_dir, group, "J", [l.j for l in loci])
    d_stems = [l.d for l in loci if l.has_d]
    d = (_combined_db(species_dir, group, "D", d_stems)  # type: ignore[arg-type]
         if d_stems else _dummy_d_db(species_dir))
    aux = bin_dir() / "optional_file" / f"{organism}_gl.aux"
    igblast.igblastn_airr(
        query_fa, out_tsv, organism=organism,
        germline_db_v=v, germline_db_j=j, germline_db_d=d,
        auxiliary_data=aux if aux.exists() else None,
        ig_seqtype=_IG_SEQTYPE[group], num_threads=num_threads)
    return out_tsv


def igblast_reads(
    query: str | Path,
    out_tsv: str | Path,
    *,
    organism: str = "human",
    num_threads: int = 1,
    groups: tuple[str, ...] | None = None,
) -> Path:
    """Run IgBLAST on reads across all annotatable loci; write a merged AIRR TSV.

    Args:
        query: FASTA/FASTQ (gzip by ``.gz``); converted to a plain FASTA for IgBLAST.
        groups: restrict to receptor types (``"TR"``/``"IG"``); default both.

    Returns:
        ``out_tsv`` — one AIRR row per read, keeping the highest ``v_score`` hit when
        a read aligns under more than one receptor type.
    """
    out_tsv = Path(out_tsv)
    imgt.download_reference()  # idempotent
    species_dir = IMGT_SPECIES_DIR[organism]
    wanted = groups or ("TR", "IG")
    run_groups = [g for g in wanted if igblast.has_internal_annotation(organism, g)]
    if not run_groups:
        raise igblast.IgBlastError(
            f"IgBLAST ships no internal annotation for organism {organism!r}")

    with tempfile.TemporaryDirectory(prefix="arda_igblast_") as td:
        tmp = Path(td)
        qfa = seqio.write_fasta(seqio.read_sequences(query), tmp / "query.fasta")
        frames = []
        for g in run_groups:
            gtsv = _run_group(qfa, organism, species_dir, g, num_threads,
                              tmp / f"{g}.tsv")
            frames.append(pl.read_csv(gtsv, separator="\t", infer_schema_length=0))

    df = pl.concat(frames, how="vertical_relaxed") if len(frames) > 1 else frames[0]
    vs = (pl.col("v_score").cast(pl.Float64, strict=False) if "v_score" in df.columns
          else pl.lit(None, dtype=pl.Float64))
    df = (df.with_columns(vs.alias("_vs"))
            .sort("_vs", descending=True, nulls_last=True)
            .unique(subset="sequence_id", keep="first")
            .drop("_vs"))
    df.write_csv(out_tsv, separator="\t")
    return out_tsv
