"""Fetch GenBank mRNA records + IgBLAST reference output as committed test fixtures.

Stores under tests/data/realworld/ (committed) so realworld tests run offline and
reproducibly. Covers a BCR locus (IGH — has D and long CDR3) and a TCR locus (TRB).

    env ARDA_MMSEQS=$(which mmseqs) python scripts/build_test_fixtures.py
"""
from __future__ import annotations

import os
import urllib.parse
import urllib.request
from pathlib import Path

from arda.paths import data_dir, bin_dir
from arda import igblast

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "realworld"

# locus -> (query, igblast V/D/J db stems, ig_seqtype, N)
SETS = {
    "IGH": ('"immunoglobulin heavy chain"[Title] AND Homo sapiens[Organism] AND mRNA[Filter]',
            ("IGHV", "IGHD", "IGHJ"), "Ig", 150),
    "TRB": ('"T cell receptor beta"[Title] AND Homo sapiens[Organism] AND mRNA[Filter]',
            ("TRBV", "TRBD", "TRBJ"), "TCR", 100),
}


def fetch(query: str, n: int) -> str:
    es = f"{EUTILS}/esearch.fcgi?db=nuccore&retmax={n}&term=" + urllib.parse.quote(query)
    ids = [s.split("</Id>")[0] for s in urllib.request.urlopen(es, timeout=120).read().decode().split("<Id>")[1:]]
    if not ids:
        raise SystemExit(f"no records for {query!r}")
    ef = f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ids)
    return urllib.request.urlopen(ef, timeout=300).read().decode()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    threads = os.cpu_count() or 1
    db = data_dir() / "blastdb" / "Homo_sapiens"
    aux = bin_dir() / "optional_file" / "human_gl.aux"
    for locus, (query, (v, d, j), seqtype, n) in SETS.items():
        fa = OUT / f"{locus.lower()}_mrna.fasta"
        fa.write_text(fetch(query, n))
        nrec = fa.read_text().count(">")
        out = OUT / f"{locus.lower()}_mrna.igblast.airr.tsv"
        igblast.igblastn_airr(fa, out, organism="human",
                              germline_db_v=db / v, germline_db_j=db / j, germline_db_d=db / d,
                              auxiliary_data=aux if aux.exists() else None,
                              ig_seqtype=seqtype, num_threads=threads)
        print(f"{locus}: {nrec} records -> {fa.name} + {out.name}")


if __name__ == "__main__":
    main()
