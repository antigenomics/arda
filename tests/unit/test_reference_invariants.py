"""Invariants the committed reference must satisfy.

These read ``database/vdj/<org>/`` straight from the repo -- no build, no mmseqs, no IgBLAST -- so
they run in CI. That matters: the two reference bugs found so far (shared TRAV/DV missing from TRD,
and 3'-truncated V alleles building scaffolds) were both invisible to CI because the only tests
covering them needed a built DB and were skipped.

Skipped cleanly when the reference is absent (a bare PyPI install with no source tree).
"""

import csv

import pytest

from arda.paths import vdj_dir

ORGANISMS = ["human", "mouse", "rat", "rabbit", "rhesus_monkey"]
_TAIL_SLACK_NT = 3


def _rows(org, name):
    p = vdj_dir(org) / name
    if not p.is_file():
        pytest.skip(f"no committed reference for {org}")
    with p.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.parametrize("org", ORGANISMS)
def test_every_scaffold_junction_starts_at_cys104(org):
    """A junction that does not begin with Cys is not a junction.

    Before the fix, 570 human V-J scaffolds (3.2%) emitted a junction starting on some other
    residue, and 549 of those were flagged ``productive=T`` -- i.e. arda was minting clonotypes
    from germlines whose Cys104 it could not even locate. Such scaffolds are still built (the read
    must still map) but must carry no junction.
    """
    bad = [r for r in _rows(org, "markup.tsv")
           if (r.get("junction_aa") or "") and not r["junction_aa"].startswith("C")]
    assert not bad, (
        f"{org}: {len(bad)} scaffolds emit a junction not starting at Cys104, "
        f"e.g. {[(r['scaffold_id'], r['v_call'], r['junction_aa']) for r in bad[:3]]}"
    )


@pytest.mark.parametrize("org", ORGANISMS)
def test_no_scaffold_is_built_from_a_truncated_allele(org):
    """Every V allele that reaches a scaffold templates as far into the junction as its gene's best.

    This is the TRAV20*03 guard: its germline stopped 10 nt short of TRAV20*01/*02, so every
    clonotype built on it came out 3 aa short -- while still looking canonical (starts C, ends
    [FW]), so nothing downstream could catch it.
    """
    anchors = _rows(org, "cdr3_anchors.tsv")
    tails = {r["allele"]: len(r["germline_nt"]) for r in anchors
             if r["segment"] == "V" and r["status"] != "no_anchor"}
    best: dict[str, int] = {}
    for a, t in tails.items():
        g = a.split("*")[0]
        best[g] = max(best.get(g, 0), t)

    used = {a for r in _rows(org, "markup.tsv") for a in (r["v_call"] or "").split(",") if a}
    bad = [(a, tails[a], best[a.split("*")[0]]) for a in used
           if a in tails and best[a.split("*")[0]] - tails[a] >= _TAIL_SLACK_NT]
    assert not bad, f"{org}: {len(bad)} truncated V alleles still build scaffolds, e.g. {bad[:3]}"


@pytest.mark.parametrize("org", ORGANISMS)
def test_no_gene_lost_every_allele(org):
    """Dropping truncated alleles must never delete a whole gene.

    Safe by construction -- the predicate is relative to the gene's own longest allele -- but assert
    it, because a future absolute threshold would quietly break this.
    """
    # Three things must be excluded before this comparison is meaningful:
    #  - pseudogenes: cdr3_anchors.tsv keeps them (VDJdb cites them) but only F/ORF alleles ever
    #    reach the scaffold builder (`imgt.load_functional_alleles`);
    #  - unanchored alleles: they never had a scaffold to begin with;
    #  - whole loci with no scaffolds: rat/rabbit/rhesus are IG-only (IgBLAST ships no TR internal
    #    annotation for them), so their TR genes legitimately have anchors and no scaffolds.
    markup = _rows(org, "markup.tsv")
    # A locus counts as "built" only if it has V-J scaffolds. J+C scaffolds carry no V and bypass
    # IgBLAST entirely, so rat/rhesus/rabbit have TRA/TRB J+C rows despite having no TR V-J
    # scaffolds at all (IgBLAST ships them no TR internal annotation).
    built_loci = {r["locus"] for r in markup if r["v_call"]}
    eligible = {r["allele"].split("*")[0] for r in _rows(org, "cdr3_anchors.tsv")
                if r["segment"] == "V" and r["status"] != "no_anchor"
                and r["functionality"] in ("F", "ORF") and r["locus"] in built_loci}
    used = {a.split("*")[0] for r in markup for a in (r["v_call"] or "").split(",") if a}
    lost = sorted(eligible - used)
    assert not lost, f"{org}: {len(lost)} V genes lost every allele: {lost[:8]}"


def test_human_trd_carries_the_shared_trav_dv_genes():
    """Regression guard for the 2.5.2 fix: δ rearrangements on a shared TRAV*/DV* V gene need a TRD
    scaffold, or they are miscalled TRA (the locus follows J/D/C, never the shared V)."""
    rows = _rows("human", "markup.tsv")
    trd_dv = [r for r in rows if r["locus"] == "TRD" and "/DV" in (r["v_call"] or "")]
    assert trd_dv, "TRD has no scaffolds on a shared TRAV*/DV* gene — the 2.5.2 fix regressed"


def test_human_has_the_jc_constant_scaffolds():
    """A read spanning the J→C splice carries no V; it needs J+C scaffolds or it cannot map."""
    jc = [r for r in _rows("human", "markup.tsv") if not r["v_call"] and r.get("c_call")]
    assert len(jc) == 345, f"expected 345 J+C scaffolds, got {len(jc)}"
