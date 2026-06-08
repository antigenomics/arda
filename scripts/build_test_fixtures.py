"""Fetch a balanced ~10k GenBank mRNA set across all 5 species + loci, with an
IgBLAST reference, stored gzipped as committed test fixtures.

Per species we write ``tests/assets/realworld/<organism>.fasta.gz`` (all loci) and
``<organism>.igblast.airr.tsv.gz`` (IgBLAST AIRR, concatenated across loci). Only
loci IgBLAST can annotate are fetched: IG for all species, TR only for human/mouse
(no TR internal annotation ships for rat/rabbit/rhesus). Runs offline at test time.

    env ARDA_MMSEQS=$(which mmseqs) python scripts/build_test_fixtures.py
"""
from __future__ import annotations

import gzip
import os
import time
import urllib.parse
import urllib.request
from pathlib import Path

from arda.paths import data_dir, bin_dir
from arda.igblast import SUPPORTED_ORGANISMS, igblastn_airr, has_internal_annotation
from arda.refbuild.loci import LOCI, IMGT_SPECIES_DIR

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = Path(__file__).resolve().parents[1] / "tests" / "assets" / "realworld"
PER_GROUP = 500          # records per (species, locus)
SPECIES_LATIN = {
    "human": "Homo sapiens", "mouse": "Mus musculus", "rat": "Rattus norvegicus",
    "rabbit": "Oryctolagus cuniculus", "rhesus_monkey": "Macaca mulatta",
}
LOCUS_TITLE = {
    "IGH": "immunoglobulin heavy chain", "IGK": "immunoglobulin kappa",
    "IGL": "immunoglobulin lambda", "TRA": "T cell receptor alpha",
    "TRB": "T cell receptor beta", "TRG": "T cell receptor gamma",
    "TRD": "T cell receptor delta",
}
LOCUS = {l.name: l for l in LOCI}


def esearch(term: str, n: int) -> list[str]:
    url = f"{EUTILS}/esearch.fcgi?db=nuccore&retmax={n}&term=" + urllib.parse.quote(term)
    xml = urllib.request.urlopen(url, timeout=120).read().decode()
    return [s.split("</Id>")[0] for s in xml.split("<Id>")[1:]]


def efetch(ids: list[str]) -> str:
    out = []
    for i in range(0, len(ids), 200):
        url = f"{EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ids[i:i+200])
        out.append(urllib.request.urlopen(url, timeout=300).read().decode())
        time.sleep(0.4)  # be polite to NCBI
    return "".join(out)


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    threads = os.cpu_count() or 1
    aux = bin_dir() / "optional_file"
    total = 0
    for org in SUPPORTED_ORGANISMS:
        sp = IMGT_SPECIES_DIR[org]
        db = data_dir() / "blastdb" / sp
        fasta_all, airr_all, header = [], [], None
        for loc in LOCI:
            if not has_internal_annotation(org, loc.group):
                continue
            term = f'"{LOCUS_TITLE[loc.name]}"[Title] AND {SPECIES_LATIN[org]}[Organism] AND mRNA[Filter]'
            try:
                ids = esearch(term, PER_GROUP)
            except Exception as e:  # noqa: BLE001
                print(f"  {org}/{loc.name}: esearch failed ({e})"); continue
            if not ids:
                continue
            fa = efetch(ids)
            tmp = OUT / f"_{org}_{loc.name}.fasta"
            tmp.write_text(fa)
            d_db = db / loc.d if loc.has_d else db / "_dummyD"
            out_tsv = OUT / f"_{org}_{loc.name}.airr.tsv"
            try:
                igblastn_airr(tmp, out_tsv, organism=org,
                              germline_db_v=db / loc.v, germline_db_j=db / loc.j,
                              germline_db_d=d_db,
                              auxiliary_data=(aux / f"{org}_gl.aux"),
                              ig_seqtype=loc.ig_seqtype, num_threads=threads)
            except Exception as e:  # noqa: BLE001
                print(f"  {org}/{loc.name}: igblast failed ({e})"); tmp.unlink(); continue
            fasta_all.append(fa)
            lines = out_tsv.read_text().splitlines()
            header = header or lines[0]
            airr_all.extend(lines[1:])
            tmp.unlink(); out_tsv.unlink()
            print(f"  {org}/{loc.name}: {fa.count('>')} records")
        if not fasta_all:
            continue
        with gzip.open(OUT / f"{org}.fasta.gz", "wt") as fh:
            fh.write("".join(fasta_all))
        with gzip.open(OUT / f"{org}.igblast.airr.tsv.gz", "wt") as fh:
            fh.write(header + "\n" + "\n".join(airr_all) + "\n")
        nrec = sum(f.count(">") for f in fasta_all)
        total += nrec
        print(f"{org}: {nrec} records -> {org}.fasta.gz")
    print(f"TOTAL: {total} records")


if __name__ == "__main__":
    build()
