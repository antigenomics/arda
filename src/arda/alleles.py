"""Shared allele-name ordering; callers supply biologically eligible names."""
from __future__ import annotations

import csv


def select_major_alleles(alleles):
    """Return gene -> major allele, preferring ``*01`` then the lowest numeric allele.

    Names without a star are valid (notably imported D gene references). This
    function imposes no functionality/completeness policy on the supplied names.
    """
    selected = {}
    def order(name):
        suffix = name.rsplit("*", 1)[1] if "*" in name else ""
        return suffix != "01", int(suffix) if suffix.isdecimal() else 0, name
    for allele in alleles:
        gene = allele.split("*", 1)[0]
        if gene not in selected or order(allele) < order(selected[gene]):
            selected[gene] = allele
    return selected


def major_allele_map(organism="human", locus=None):
    """Return usable runtime gene -> major allele for an organism or one locus.

    Prefer ok F/ORF anchors among available V/J targets. Genes with only
    incomplete targets retain one preliminary alignment template, with its
    original anchor status; this never establishes complete emission eligibility.
    D names come from exported germlines. No raw IMGT installation is needed.
    """
    from .paths import vdj_dir
    base = vdj_dir(organism)
    with (base / "markup.tsv").open() as handle:
        available = {allele for row in csv.DictReader(handle, delimiter="\t")
                     if locus is None or row["locus"] == locus
                     for key in ("v_call", "j_call") for allele in row.get(key, "").split(",") if allele}
    with (base / "cdr3_anchors.tsv").open() as handle:
        eligible = {row["allele"] for row in csv.DictReader(handle, delimiter="\t")
                    if row["allele"] in available and row["status"] == "ok"
                    and row["functionality"] in ("F", "ORF")
                    and (locus is None or row["locus"] == locus)}
    d_names = set()
    with (base / "d_germlines.fasta").open() as handle:
        for line in handle:
            if line.startswith(">"):
                loc, allele = line[1:].strip().split("|", 1)
                if locus is None or loc == locus:
                    d_names.add(allele)
    available.update(d_names)
    d_metadata = base / "d_germlines.tsv"
    if d_metadata.exists():
        with d_metadata.open() as handle:
            eligible.update(row["allele"] for row in csv.DictReader(handle, delimiter="\t")
                            if row["allele"] in d_names and row["status"] == "ok")
    else:
        eligible.update(d_names)
    selected = select_major_alleles(available)
    selected.update(select_major_alleles(eligible))
    return selected


def d_germline_metadata(organism="human"):
    """Return exact D-template metadata keyed by (allele, nucleotide sequence).

    Preliminary alignment targets may be partial; complete emission consumers
    must reject status=partial. Functionality policy belongs to the consumer.
    """
    from .paths import vdj_dir
    with (vdj_dir(organism) / "d_germlines.tsv").open() as handle:
        return {(row["allele"], row["sequence"]): row for row in csv.DictReader(handle, delimiter="\t")}
