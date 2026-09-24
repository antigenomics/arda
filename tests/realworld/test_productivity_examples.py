"""The worked examples in ``docs/productivity.rst`` and the README, pinned.

Those two pages teach ``productive`` / ``stop_codon`` / ``vj_in_frame`` with named GenBank
reads and quote their translations. A documentation example nobody runs is how the page rots,
so every claim the pages make about these reads is asserted here.

The load-bearing one is ``JX277388.1``: it carries two nucleotides ``TRBJ2-5*01`` does not,
so the J's 5' germline residues and its own FR4 cannot both sit in the germline frame. arda
lands on the FR4 side, which is the side the constant exon splices to -- splice ``TRBC2*01``
on and the whole thing translates without a stop. That read is the reason the page exists:
matching a junction's tail against germline is NOT a frame test, and reading it as one
produced a wrong verdict on this exact row once already.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from arda.annotate.mapper import annotate_records
from tests.conftest import requires_mmseqs

pytestmark = [requires_mmseqs]

ROOT = Path(__file__).resolve().parents[2]

# (id, organism, sequence, j_call, junction_aa, fwr4_aa)  -- all documented on the page.
CASES = [
    ("JX277388.1", "mouse",
     "AACATGAGCCAGGGCAGAACCTTGTACTGCACCTGCAGTGCAGGGGGCCCAAGACACCCAGTACTTTTTGGGCCAGGCACTCGG",
     "TRBJ2-5*01", "CTCSAGGPRHPVLF", "FGPGTR"),
    ("JX277345.1", "mouse",
     "AAAACAAACCAGACATCTGTGTACTTCTGTGCTAGCAGTTTATTGAGGCAGGGGGAGGTGGAAATACGCTCTATTTTTGGAGAAGG"
     "AAGCCGG",
     "TRBJ1-3*01", "CASSLLRQGEVEIRSIF", "FGEGSR"),
    ("MH918759.1", "human",
     "GAAGATGGAAGGTTTACAGCACAGCTCAATAGAGCCAGCCAGTATATTTCCCTGCTCATCAGAGACTCCAAGCTCAGTGATTCAGC"
     "CACCTACCTCTGTGTGGTGAACATTCGGGGAAATGCCAGACTCATGTTTGGAGATGGAACTCAGCTGGTGGTGAAGCCCAATATCC"
     "AGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTC",
     "TRAJ31*01", "CVVNIRGNARLMF", "FGDGTQLVVKP"),
]

CODONS = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        [(a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}


def _translate(nt: str) -> str:
    return "".join(CODONS.get(nt[i:i + 3].upper(), "X") for i in range(0, len(nt) - 2, 3))


def _fasta(path: Path) -> dict[str, str]:
    out, key = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            out[key] = ""
        elif key:
            out[key] += line.strip()
    return out


def _germline_j_tail(org: str, j_call: str) -> str:
    """Germline nucleotides from the [FW]118 codon to the J's 3' end, from COMMITTED artifacts.

    Never ``database/vdj/<org>/segments.fasta``: it is GENERATED from the reference, is not
    committed, and the fetched reference tarball does not carry it -- so a test reading it passes
    on a built checkout and fails everywhere else. ``alleles.fasta`` and ``markup.tsv`` are both
    committed, and ``fwr4_start .. vj_end`` on any scaffold carrying this J is exactly the span
    wanted: FR4 **plus** the trailing partial codon the C exon completes. ``markup.tsv``'s ``fwr4``
    column alone stops at the last whole codon and would shift the splice by one base.
    """
    base = ROOT / "database" / "vdj" / org
    scaffolds = _fasta(base / "alleles.fasta")
    rows = (base / "markup.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    sid, fwr4_start, vj_end = (header.index(c) for c in ("scaffold_id", "fwr4_start", "vj_end"))
    for row in rows[1:]:
        f = row.split("\t")
        if f[header.index("j_call")] == j_call and f[fwr4_start] and f[vj_end]:
            return scaffolds[f[sid]][int(f[fwr4_start]) - 1:int(f[vj_end])]
    raise AssertionError(f"{j_call} carries no markup row in {org}")


@pytest.mark.parametrize("sid,org,seq,j_call,junction_aa,fwr4_aa", CASES, ids=[c[0] for c in CASES])
def test_documented_call(sid, org, seq, j_call, junction_aa, fwr4_aa):
    """The page quotes these columns verbatim; they must still be what arda writes."""
    rec = annotate_records([(sid, seq)], organism=org)[0]
    assert rec["j_call"].split(",")[0] == j_call
    assert rec["junction_aa"] == junction_aa
    assert rec["fwr4_aa"] == fwr4_aa
    assert rec["vj_in_frame"] == "T"
    assert rec["stop_codon"] == "F"
    assert rec["productive"] == "T"


@pytest.mark.parametrize("sid,org,seq,j_call,_ja,_f4", CASES, ids=[c[0] for c in CASES])
def test_constant_exon_translates(sid, org, seq, j_call, _ja, _f4):
    """``vj_in_frame`` means the constant exon spliced onto that J translates.

    That is the only sound test -- comparing the junction's tail to the J's germline residues
    measures hypermutation, not frame. Complete the J from germline, append C, translate from
    Cys104 in arda's frame: there must be no stop codon.
    """
    rec = annotate_records([(sid, seq)], organism=org)[0]
    query = rec["sequence"]
    if rec["rev_comp"] == "T":
        query = query.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

    germline_j = _germline_j_tail(org, j_call)

    # A read that already runs past the J carries its own constant region (MH918759.1 has 66 nt
    # of TRAC) and needs no reconstruction. A read that stops inside the J is completed from
    # germline and spliced to the C exon -- both preserve arda's frame, which is the point.
    tail = ""
    for length in range(min(len(query), len(germline_j)), 5, -1):
        at = germline_j.find(query[-length:])
        if at >= 0:
            c_gene = _fasta(ROOT / "database" / "c_genes" / f"{org}.fasta")
            # TRBC1*01 and TRBC2*01 encode the same protein; TRBC1's stored exon carries one
            # extra leading nucleotide, a property of the reference rather than of the read.
            tail = germline_j[at + length:] + c_gene[
                "TRBC2*01" if rec["locus"] == "TRB" else "TRAC*01"]
            break
    else:
        assert germline_j[:12] in query, f"{sid}: read neither ends in nor spans the germline J"

    cys104 = int(rec["cdr3_start"]) - 3
    protein = _translate(query[cys104 - 1:] + tail)
    assert "*" not in protein, f"{sid}: frameshift into C -- {protein}"
    assert protein.startswith(rec["junction_aa"])


def test_junction_tail_vs_germline_is_not_a_frame_test():
    """JX277388.1 is the counterexample the page is built on.

    Its junction ends ``HPVLF``, nothing like ``TRBJ2-5*01``'s germline ``NQDTQYF`` -- yet the
    read is in frame and the receptor translates. Anything that infers non-productivity from
    that mismatch is wrong on this row.
    """
    sid, org, seq, *_ = CASES[0]
    rec = annotate_records([(sid, seq)], organism=org)[0]
    anchors = (ROOT / "database" / "vdj" / org / "cdr3_anchors.tsv").read_text().splitlines()
    header = anchors[0].split("\t")
    row = next(r.split("\t") for r in anchors[1:] if r.split("\t")[2] == "TRBJ2-5*01")
    templated = row[header.index("templated_aa")]

    assert templated == "NQDTQYF"
    assert not rec["junction_aa"].endswith(templated[:-1])   # the tails do NOT match
    assert rec["vj_in_frame"] == "T"                          # and it is still in frame
