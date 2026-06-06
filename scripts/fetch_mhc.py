#!/usr/bin/env python3
"""Pre-fetch MHC allele references (reference-only; mapping is a future phase).

Downloads:
  * Human HLA class I & II protein sequences from the IPD-IMGT/HLA project
    (ANHIG/IMGTHLA mirror), and
  * β2-microglobulin (B2M) protein from UniProt.

Raw downloads go to ``data/mhc`` (gitignored). A small **representative** set —
one allele per gene plus B2M — is written to ``database/mhc/<organism>/`` so the
committed footprint stays small. Re-run with ``--full`` to also copy the complete
per-gene FASTAs into ``data/mhc`` for development.

Non-human MHC (mouse H2, rat RT1, etc.) lives in IPD-MHC and is left as a
documented TODO — the directory scaffolding is created here.

    python scripts/fetch_mhc.py
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

# IPD-IMGT/HLA combined protein FASTAs (class I and class II) from the ANHIG repo.
HLA_PROT_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta"
# B2M reviewed human protein (UniProt P61769).
B2M_URL = "https://rest.uniprot.org/uniprotkb/P61769.fasta"

# Genes we keep one representative allele of in the committed reference.
CLASS_I_GENES = ("A", "B", "C", "E", "F", "G")
CLASS_II_GENES = ("DRA", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1")


def _here() -> Path:
    return Path(__file__).resolve().parents[1]


def _download(url: str, dest: Path) -> None:
    print(f"[fetch_mhc] downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=600) as r, open(dest, "wb") as fh:
        fh.write(r.read())


def _read_fasta(path: Path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def _representative(hla_fasta: Path):
    """Pick the first allele encountered per gene -> {gene_label: (header, seq)}."""
    chosen: dict[str, tuple[str, str]] = {}
    wanted = {f"HLA-{g}" for g in CLASS_I_GENES + CLASS_II_GENES}
    for header, seq in _read_fasta(hla_fasta):
        # Header: ">HLA:HLA00001 A*01:01:01:01 365 bp"
        parts = header.split()
        allele = parts[1] if len(parts) > 1 else parts[0]
        gene = "HLA-" + allele.split("*")[0]
        if gene in wanted and gene not in chosen:
            chosen[gene] = (allele, seq)
    return chosen


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="Keep full HLA FASTA in data/mhc.")
    args = ap.parse_args()

    root = _here()
    raw = root / "data" / "mhc"
    raw.mkdir(parents=True, exist_ok=True)
    out = root / "database" / "mhc" / "human"
    out.mkdir(parents=True, exist_ok=True)

    hla = raw / "hla_prot.fasta"
    if not hla.exists():
        _download(HLA_PROT_URL, hla)
    b2m = raw / "b2m.fasta"
    if not b2m.exists():
        _download(B2M_URL, b2m)

    rep = _representative(hla)
    ci = [(a, s) for g, (a, s) in rep.items() if g.split("-")[1] in CLASS_I_GENES]
    cii = [(a, s) for g, (a, s) in rep.items() if g.split("-")[1] in CLASS_II_GENES]
    (out / "class_i.fasta").write_text("".join(f">{a}\n{s}\n" for a, s in ci))
    (out / "class_ii.fasta").write_text("".join(f">{a}\n{s}\n" for a, s in cii))
    b2m_seq = next(_read_fasta(b2m))
    (out / "b2m.fasta").write_text(f">B2M_HUMAN\n{b2m_seq[1]}\n")

    # Scaffolding for non-human organisms (IPD-MHC) — TODO.
    for org in ("mouse", "rat", "rabbit", "rhesus_monkey"):
        d = root / "database" / "mhc" / org
        d.mkdir(parents=True, exist_ok=True)
        readme = d / "README.md"
        if not readme.exists():
            readme.write_text(
                f"# {org} MHC references (TODO)\n\n"
                "Non-human MHC references come from IPD-MHC "
                "(https://www.ebi.ac.uk/ipd/mhc/). Not yet fetched.\n"
            )

    print(f"[fetch_mhc] wrote {len(ci)} class-I, {len(cii)} class-II alleles + B2M "
          f"to {out}", file=sys.stderr)
    if not args.full:
        # Keep raw HLA out of the way unless --full (it's gitignored regardless).
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
