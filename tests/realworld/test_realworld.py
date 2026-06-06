"""Real-world concordance test: arda vs IgBLAST on GenBank mRNA.

Fetches human immunoglobulin heavy-chain mRNA records from NCBI, annotates them
both with IgBLAST (the gold standard) and with arda, and checks that the FR/CDR
amino-acid region calls agree. Guarded by ``ARDA_REALWORLD=1`` (needs network +
the human reference DB + the per-locus germline BLAST DBs from the build).

    env ARDA_REALWORLD=1 ARDA_MMSEQS=$(which mmseqs) \\
        python -m pytest tests/realworld -s
"""

from __future__ import annotations

import os
import urllib.parse
import urllib.request
from pathlib import Path

import pytest
import polars as pl

from arda.paths import data_dir, vdj_dir
from arda import igblast
from arda.annotate.mapper import annotate_records
from tests.conftest import requires_mmseqs, requires_human_db

pytestmark = [
    pytest.mark.skipif(not os.getenv("ARDA_REALWORLD"), reason="set ARDA_REALWORLD=1"),
    requires_mmseqs,
    requires_human_db,
]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = '"immunoglobulin heavy chain"[Title] AND Homo sapiens[Organism] AND mRNA[Filter]'
N_RECORDS = 50
REGIONS = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3")


def _fetch_mrna(n: int) -> Path:
    """Fetch n IGH mRNA records as FASTA (cached under data/realworld)."""
    out_dir = data_dir() / "realworld"
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / "igh_mrna.fasta"
    if fasta.exists() and fasta.stat().st_size > 0:
        return fasta
    es = f"{EUTILS}/esearch.fcgi?db=nuccore&retmax={n}&term=" + urllib.parse.quote(QUERY)
    with urllib.request.urlopen(es, timeout=120) as r:
        xml = r.read().decode()
    ids = [s.split("</Id>")[0] for s in xml.split("<Id>")[1:]]
    if not ids:
        pytest.skip("NCBI returned no records")
    ef = f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ids)
    with urllib.request.urlopen(ef, timeout=300) as r:
        fasta.write_bytes(r.read())
    return fasta


def _igblast_airr(fasta: Path) -> pl.DataFrame:
    db = data_dir() / "blastdb" / "Homo_sapiens"
    aux = vdj_dir().parent.parent / "bin" / "optional_file" / "human_gl.aux"
    out = fasta.with_suffix(".igblast.airr.tsv")
    igblast.igblastn_airr(
        fasta, out, organism="human",
        germline_db_v=db / "IGHV", germline_db_j=db / "IGHJ", germline_db_d=db / "IGHD",
        auxiliary_data=aux if Path(aux).exists() else None,
        ig_seqtype="Ig", num_threads=os.cpu_count() or 1,
    )
    return pl.read_csv(out, separator="\t", infer_schema_length=0)


def test_arda_matches_igblast_regions():
    from arda.annotate.io import read_sequences

    fasta = _fetch_mrna(N_RECORDS)
    queries = list(read_sequences(fasta))
    igb = {r["sequence_id"]: r for r in _igblast_airr(fasta).iter_rows(named=True)}
    arda = {r["sequence_id"]: r for r in annotate_records(queries, "human", "nt", threads=4)}

    compared = agree = 0
    for sid, ig in igb.items():
        ar = arda.get(sid)
        if not ar:
            continue
        for r in REGIONS:
            ig_aa = (ig.get(f"{r}_aa") or "").strip()
            ar_aa = (ar.get(f"{r}_aa") or "").strip()
            if not ig_aa or not ar_aa:
                continue
            compared += 1
            # Terminal regions (FR1 start, CDR3 ends) legitimately differ by one
            # boundary residue between tools, so accept substring containment there;
            # internal regions must match exactly.
            if r in ("fwr1", "cdr3"):
                agree += ig_aa in ar_aa or ar_aa in ig_aa
            else:
                agree += ig_aa == ar_aa
    assert compared > 0, "no comparable regions"
    frac = agree / compared
    print(f"\n[realworld] region concordance: {agree}/{compared} = {frac:.1%}")
    assert frac >= 0.90
