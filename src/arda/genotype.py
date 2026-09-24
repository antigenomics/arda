"""The V alleles one donor actually carries, and restricting calls to them.

A germline reference is a catalogue of every allele anyone has, and no donor has all of them. Two
alleles of a gene at most -- so a call naming a third is wrong before any sequence is looked at,
and a read that could not choose between two alleles the donor does not carry was never ambiguous
in the first place. Immcantation's TIgGER measured the size of this on full-length BCR: restricting
V calls to an inferred genotype took ambiguous assignments from 11.2 % to 1.5 %.

What this module is, and is not
-------------------------------
**Never: this is a claim about the REFERENCE, never about the repertoire.** It answers "which
germline alleles does this donor carry", which is annotation -- the same kind of question as "which
J does this read use". It does not answer anything about clonal composition, and it must not grow
into one; that line is ``vdjtools``'.

It is also not :data:`arda.stats.ALLELE_MIN_FRAC`'s ``allele_candidate`` scope. That is a shortlist
of recurrent high-quality V mutations to *look at*, deliberately never a call, and it stays that
way. This is a separate, narrower object with its own evidence and its own file.

**Never: a restriction NEVER re-aligns and never rebuilds a reference.** TIgGER's ``reassignAlleles``
re-assigns rather than re-aligning, and so does this, for reasons specific to arda: scaffold ids are
positional, so changing the allele set renumbers every scaffold in the locus; ``build-db`` needs
IgBLAST and IMGT network access that a ``pip install`` does not have; and the mmseqs freshness
contract is an mtime with no allele-set identity recorded anywhere, so a donor-specific FASTA
dropped beside the shared one would silently invalidate the index for every concurrent process.
Given the span a read already aligned over, the restriction is a set intersection.

The file
--------
One row per carried allele -- ``locus, gene, allele, votes, gene_votes, explained, note`` -- because
a cell holding ``TRBV19*01,TRBV19*03`` cannot be sorted, joined or summed. A gene that could not be
called gets a row with an empty ``allele`` and a ``note`` saying why: **omitted with a reason, never
silently absent.**

Only ``allele`` is required when reading, and the gene is derived from it, so the simplest possible
user-supplied genotype is a one-column TSV of allele names -- which is the point. An OGRDB set, a
library inferred by MiXCR, or a list typed by hand are all first-class inputs, and none of them
needs arda to have inferred anything.

How the call is made
--------------------
**The junction is de facto a UMI.** V(D)J junctional diversity makes a nucleotide junction
essentially unique to one rearrangement, so grouping reads by ``(locus, gene, J gene, junction)``
groups them by molecule -- and each distinct junction is one independent draw from the donor's two
chromosomes. Two things follow, and the inference rests on both:

* **The clonotype, not the read, is the unit of observation.** One expanded clone is one
  chromosome's worth of evidence however many reads it has. Read depth is reported beside the
  clonotype count rather than substituted for it.
* **Disagreement WITHIN a junction is error, not allele.** Every read of one rearrangement carries
  the same V allele by construction, so a read naming a different one is hypermutation or
  sequencing error. That is where :func:`error_rate` comes from -- the miscall rate is measured on
  the library itself instead of being a constant somebody picked.

The call is then a **likelihood ratio between diploid genotypes**: every single allele and every
pair is scored by its multinomial likelihood under that error rate, and the winner is called only
if it beats the runner-up by ``--min-log10-bf``. Homozygous and heterozygous are the same formula
at genotype size 1 and 2, so neither is special-cased.

⚠ **Why it is a likelihood and not a coverage rule.** The first version of this was TIgGER's
frequency rule -- the fewest alleles explaining 7/8 of the calls -- which has no error model and so
cannot tell overwhelming evidence from none. Measured on a real TRB amplicon it called
``TRBV11-2`` off **754 of 757** clonotypes that singled the allele out and ``TRBV20-1`` off
**1 of 2,544**, reported both as ``explained = 1.0000, ok``, and returned 43 of 53 genes with not
one heterozygous. A ratio of counts is not a confidence.

⚠ **Allele-level genotyping is a read-length feature, and most libraries do not have it.**
Separating a gene's alleles needs a median of 150 nt of TRBV measured from the 3' end (175 TRAV,
230 IGHV). The TRB amplicon above is 151 nt paired and its reads cover a median of **72 nt** of V
germline -- **56 nt** after clipping at the Cys104 anchor, and **not one read** reaches 150. So on
that library most genes genuinely cannot be resolved, and the right output is a refusal. Genes the
data cannot separate are reported ``undetermined`` with no alleles: omitted with a reason, never
guessed.
"""

from __future__ import annotations

import math
import typing
from collections import Counter
from pathlib import Path

import polars as pl

from ._log import logger
from .germline import segment_germlines

__all__ = ["GENOTYPE_COLUMNS", "MIN_LOG10_BF", "MIN_CLONOTYPES", "MIN_ERROR_RATE", "PLOIDY",
           "FRACTION_TO_EXPLAIN", "Clonotype", "gene_of", "error_rate", "read_genotype",
           "write_genotype", "restrict", "restrict_airr", "infer_genotype"]

#: Written in this order. ``allele`` is the only one a reader requires.
GENOTYPE_COLUMNS = ("locus", "gene", "allele", "clonotypes", "reads", "gene_clonotypes",
                    "log10_bf", "note")


def gene_of(allele: str) -> str:
    """``TRBV19*01`` -> ``TRBV19``. Same rule as :func:`arda.dpost._gene`."""
    return allele.split("*")[0]


def read_genotype(path, *, organism: str | None = "human") -> dict[str, tuple[str, ...]]:
    """``{gene: carried alleles}`` from a genotype TSV.

    Rows with an empty ``allele`` are the "could not call this gene" rows and are skipped, which is
    what makes a gene *absent* from the returned mapping mean "not genotyped" rather than "carries
    nothing" -- :func:`restrict` leaves such a gene's calls alone.

    ``organism`` validates every named allele against the reference and **raises** on one the
    reference cannot call. Never: this is the user-supplied-list path, and a typo or an allele arda
    has no scaffold for must stop the run rather than silently restrict nothing -- 884 human V
    alleles are functional in IMGT and 801 reach a scaffold, so naming one of the others is an easy
    and completely invisible mistake. Pass ``organism=None`` to skip the check.
    """
    df = pl.read_csv(path, separator="\t", quote_char=None, infer_schema_length=0,
                     comment_prefix="#")
    if "allele" not in df.columns:
        raise ValueError(f"{path}: a genotype needs an `allele` column; got {list(df.columns)}")

    known = set(segment_germlines(organism, "v")) if organism else None
    out: dict[str, list[str]] = {}
    unknown: list[str] = []
    for allele in df["allele"].to_list():
        allele = (allele or "").strip()
        if not allele:
            continue
        if known is not None and allele not in known:
            unknown.append(allele)
            continue
        out.setdefault(gene_of(allele), []).append(allele)
    if unknown:
        raise ValueError(
            f"{path}: {len(unknown)} allele(s) are not in the {organism} reference and could "
            f"never be called: {', '.join(sorted(set(unknown))[:6])}"
            + (" ..." if len(set(unknown)) > 6 else ""))
    # Sorted, so a restriction's output does not depend on the order rows were written in.
    return {gene: tuple(sorted(set(alleles))) for gene, alleles in out.items()}


def write_genotype(rows, output, *, params: dict | None = None) -> int:
    """Write genotype rows as a TSV, with the parameters that produced them as ``#`` comments.

    The provenance header carries the allele universe the call was made over, because a genotype is
    only meaningful against the reference it was inferred from.
    """
    output = Path(output)
    frame = pl.DataFrame(
        list(rows),
        schema={c: pl.String for c in GENOTYPE_COLUMNS},
        orient="row",
    ) if not isinstance(rows, pl.DataFrame) else rows
    header = "".join(f"# {k}: {v}\n" for k, v in (params or {}).items())
    with open(output, "w") as fh:
        fh.write(header)
        frame.write_csv(fh, separator="\t", quote_style="never")
    return frame.height


def restrict(candidates: tuple[str, ...] | None, call: str,
             genotype: dict[str, tuple[str, ...]]) -> str:
    """The subset of ``candidates`` the donor carries, as a call string.

    Three outcomes, and they are deliberately distinguishable in one column:

    ``call`` unchanged
        ``candidates`` is ``None`` -- the tie machinery declined (span below ``MIN_SPAN``, allele
        absent, tie list past ``max_ties``) -- or no candidate's gene was genotyped at all.
        **Never: a refusal to answer is not a contradiction.** Treating the two alike would delete
        the call of every short read in the library.
    a narrower call
        the usual case, and the whole point.
    empty
        every candidate belongs to a genotyped gene and none is carried. This read contradicts the
        genotype. It is a COLUMN value and not a dropped row so a consumer can answer it per row,
        and because the per-gene rate of it is the sharpest available check that the genotype is
        right: a genotype that lost a carried allele makes one gene's empty rate spike while its
        neighbours sit at the error floor.
    """
    if candidates is None:
        return call
    # Never: an allele survives unless its OWN gene was genotyped and did not name it. A tie list
    # routinely spans genes -- the whole reason `resolve-ties` exists is that IGLV2-14 and
    # IGLV2-23 are indistinguishable over 70 nt -- so asking "is this allele carried" gene-blind
    # deletes every candidate belonging to a gene the genotype says nothing about. Measured on the
    # 453-read fixture with a one-gene genotype: 79 rows came back empty, and all of them were
    # reads of other genes whose tie list merely brushed the genotyped one.
    kept = [a for a in candidates if gene_of(a) not in genotype or a in genotype[gene_of(a)]]
    if kept:
        return ",".join(kept)
    return ""


def restrict_airr(path, out, genotype: dict[str, tuple[str, ...]], *, organism: str = "human",
                  column: str = "v_call_genotyped", echo=None) -> dict:
    """Add ``column`` to an AIRR TSV: ``v_call`` restricted to ``genotype``.

    ``v_call`` itself is left byte-identical. A genotype is an inference and this keeps the
    evidence it was made against in the same file, which is also what lets a reader re-restrict at
    a different genotype without re-running anything.
    """
    from .annotate.airr_out import read_airr
    from .annotate.ties import TieResolver

    df = read_airr(path)
    for col in ("v_call", "v_germline_start", "v_germline_end"):
        if col not in df.columns:
            raise ValueError(f"{path}: restriction needs {col}; got {len(df.columns)} columns")

    res = TieResolver(segment_germlines(organism, "v"))
    calls = df["v_call"].to_list()
    starts, ends = df["v_germline_start"].to_list(), df["v_germline_end"].to_list()
    values = [restrict(res.candidates(c or "", a, b), c or "", genotype)
              for c, a, b in zip(calls, starts, ends)]

    # Never: a row that never had a `v_call` has no input to a restriction, so it is counted in
    # `no_call` and NOT in `contradicted` -- "a metric with no input is OMITTED, never 0". On the
    # 453-read fixture 74 rows are J-only; scoring them as contradictions made a correct
    # restriction look like it had rejected a sixth of the library.
    assessed = [(a or "", b) for a, b in zip(calls, values) if (a or "")]
    report = {
        "rows": df.height,
        "no_call": df.height - len(assessed),
        "narrowed": sum(1 for a, b in assessed if b and b != a),
        "contradicted": sum(1 for _a, b in assessed if not b),
        "unchanged": sum(1 for a, b in assessed if b == a),
        "genes": len(genotype),
    }
    df.with_columns(pl.Series(column, values)).write_csv(out, separator="\t", quote_style="never")
    (echo or logger.info)(
        f"{column}: {report['narrowed']} narrowed, {report['contradicted']} contradicted, "
        f"{report['unchanged']} unchanged over {report['genes']} genotyped genes"
        + (f"; {report['no_call']} rows had no v_call" if report["no_call"] else ""))
    return report


# ── inference ─────────────────────────────────────────────────────────────────────────────────
#
# The call is a LIKELIHOOD RATIO between diploid genotypes, not a coverage rule.
#
# The first version of this was TIgGER's frequency rule -- the fewest alleles explaining 7/8 of the
# calls -- and it has no error model at all, so it cannot tell overwhelming evidence from none.
# Measured on a real TRB amplicon it called `TRBV11-2` off **754 of 757** clonotypes that singled
# the allele out and `TRBV20-1` off **1 of 2,544**, reported both as `explained = 1.0000, ok`, and
# called 43 of 53 genes with not one heterozygous. A ratio of counts is not a confidence.

#: Smallest ``log10`` Bayes factor between the best genotype and the runner-up before a gene is
#: called. 1.0 is "ten times more likely", the conventional "substantial" rung. Below it the gene
#: is reported ``undetermined`` -- which on a short-read library is the common and correct answer.
#: Written into the output file, and the per-gene factor is emitted beside every call, so the
#: number can be re-derived at another value without re-running anything.
MIN_LOG10_BF = 1.0

#: Below this many assigned clonotypes a gene is reported ``low_support`` with NO alleles.
#: Never: reported, not dropped -- a gene missing from the table is indistinguishable from one that
#: was never in the reference.
MIN_CLONOTYPES = 10

#: Two alleles per gene. The real constraint, not a proxy for it.
PLOIDY = 2

#: Floor on the estimated per-clonotype allele-miscall rate. A library whose clonotypes never
#: disagree internally would otherwise give an error rate of 0, which makes every homozygous
#: likelihood infinite and every call certain.
MIN_ERROR_RATE = 1e-4

#: Fraction of a gene's clonotype vote SETS the called genotype must still intersect. This is not
#: the call -- the Bayes factor is -- it is the guard that keeps an ambiguous-but-real allele from
#: being excluded. See :func:`_call_gene`.
FRACTION_TO_EXPLAIN = 0.875


class Clonotype(typing.NamedTuple):
    """One rearrangement's worth of evidence about which V allele it carries.

    Never: a clonotype is ONE rearrangement, so all of its reads carry the SAME V allele by
    construction. Reads that disagree are somatic hypermutation or sequencing error, never
    allelic variation -- which is what makes ``discordant`` an error measurement rather than a
    signal, and what makes the clonotype (not the read) the independent observation for a diploid
    call. In a UMI library the clonotype is the molecule.
    """

    locus: str
    gene: str
    allele: str | None          # the majority unambiguous assignment, or None
    compatible: frozenset       # every allele the reads could not rule out
    reads: int                  # reads behind it -- depth, reported alongside the clonotype count
    discordant: int             # reads disagreeing with the majority -> the error term
    observed: int               # reads that named ONE allele, mutated or not: the error term's n


def _clonotypes(df, res, *, unmutated_only: bool, scope: str) -> list[Clonotype]:
    """Collapse an AIRR into one record per clonotype.

    Grouping is on the Stage-1 AIRR itself -- ``(locus, gene, j gene, junction)``, ``correct``'s
    clonotype key with the V taken at GENE level so the allele being decided is not part of the
    key. No ``--read-map`` join and no second input file.
    """
    cols = df.columns
    has_anchor, has_mut = "v_anchor_nt" in cols, "v_mutations" in cols
    acc: dict[tuple, dict] = {}
    for row in df.iter_rows(named=True):
        call = (row.get("v_call") or "").strip()
        if not call:
            continue
        # Never: `candidates` slices the GERMLINE of the called allele -- the read's own bases never
        # enter it -- so a read carrying an error at a discriminating base expands exactly like one
        # that matches perfectly and votes just as confidently for the wrong allele. This filter is
        # what ties a vote to the read, and for TCR (no SHM) the rule is simply: germline-exact.
        mutated = bool(has_mut and (row.get("v_mutations") or "").strip())
        if unmutated_only and mutated:
            # Still counted for the error estimate below -- see `slot["seen"]`.
            pass
        end = row.get("v_germline_end")
        if scope == "framework" and has_anchor:
            try:
                anchor = int(row.get("v_anchor_nt") or 0)
                if anchor > 0:
                    end = min(int(end), anchor)
            except (TypeError, ValueError):
                pass
        cand = res.candidates(call, row.get("v_germline_start"), end)
        if cand is None:
            continue
        gene = gene_of(call.split(",")[0].strip())
        within = frozenset(a for a in cand if gene_of(a) == gene)
        if not within:
            continue
        key = (row.get("locus") or "", gene,
               gene_of((row.get("j_call") or "").split(",")[0].strip()),
               row.get("junction") or "")
        slot = acc.setdefault(key, {"reads": 0, "single": Counter(), "sets": [],
                                    "seen": Counter()})
        # Never: the ERROR estimate reads every read, the ASSIGNMENT reads only germline-exact
        # ones. Estimating the miscall rate from the unmutated subset is circular -- those reads
        # were selected for carrying no mismatch, so they essentially never disagree and the
        # estimate collapses onto its own floor, silently. Measured on a real TRB amplicon: 0
        # discordant reads out of 77,345, i.e. no estimate at all dressed up as a small one.
        if len(within) == 1:
            slot["seen"][next(iter(within))] += 1
        if unmutated_only and mutated:
            continue
        slot["reads"] += 1
        slot["sets"].append(within)
        if len(within) == 1:
            slot["single"][next(iter(within))] += 1

    out = []
    for (locus, gene, _j, _junction), slot in acc.items():
        single, seen = slot["single"], slot["seen"]
        if not slot["sets"]:
            continue
        allele = None
        if single:
            # Majority, ties broken by name -- a total order, so two runs cannot disagree.
            allele = min(single, key=lambda a: (-single[a], a))
        # Within one junction every read is the same rearrangement and so the same allele; a read
        # naming another is hypermutation or sequencing error by construction. Counted over ALL
        # reads of the clonotype, not just the germline-exact ones.
        discordant = (sum(seen.values()) - max(seen.values())) if seen else 0
        # Reads of one clonotype tile different windows of the same V, so their compatible sets
        # INTERSECT: two half-informative reads can pin an allele neither could alone. An empty
        # intersection means the reads disagree, and falls back to the union rather than claiming
        # more than the evidence supports.
        compatible = frozenset.intersection(*slot["sets"]) or frozenset.union(*slot["sets"])
        out.append(Clonotype(locus, gene, allele, compatible, slot["reads"], discordant,
                             sum(seen.values())))
    return out


def error_rate(clonotypes) -> float:
    """Per-observation allele-miscall rate, estimated from WITHIN-clonotype disagreement.

    Never: within one clonotype every read comes from the same rearrangement and therefore the same
    allele, so a read that names a different one is hypermutation or sequencing error by
    construction. That makes this an error rate measured on the library itself rather than a
    constant anyone had to pick -- the same reason arda ships no QC threshold.
    """
    discordant = sum(c.discordant for c in clonotypes)
    observed = sum(c.observed for c in clonotypes)
    return max(discordant / observed, MIN_ERROR_RATE) if observed else MIN_ERROR_RATE


def _log_likelihood(counts: Counter, genotype: tuple[str, ...], alleles: tuple[str, ...],
                    eps: float) -> float:
    """Multinomial log-likelihood of the observed per-allele clonotype counts under ``genotype``.

    An allele in the genotype carries ``(1 - eps) / |genotype|`` of the mass; every other allele of
    the gene shares ``eps``. Homozygous and heterozygous are the same formula at ``|genotype|`` 1
    and 2, so nothing has to special-case which is being tested.
    """
    others = len(alleles) - len(genotype)
    inside = (1.0 - eps) / len(genotype)
    outside = (eps / others) if others > 0 else 0.0
    total = 0.0
    for allele in alleles:
        n = counts.get(allele, 0)
        if not n:
            continue
        p = inside if allele in genotype else outside
        if p <= 0.0:
            return float("-inf")
        total += n * math.log(p)
    return total


def _call_gene(alleles: tuple[str, ...], clonotypes: list[Clonotype], eps: float, *,
               min_log10_bf: float, min_clonotypes: int, ploidy: int,
               fraction_to_explain: float):
    """``(chosen, log10_bf, note)`` for one gene.

    Every diploid genotype -- each single allele, then each pair -- is scored by its multinomial
    likelihood under an error rate measured from the library, and the winner is called only if it
    beats the runner-up by ``min_log10_bf``. **A Bayes factor is what a coverage ratio was
    pretending to be:** 754 of 757 clonotypes naming one allele is overwhelming, 1 of 2,544 is
    nothing, and only a likelihood tells them apart.

    Never: the ambiguous clonotypes still constrain, and dropping them is not the conservative
    choice. If a donor is ``*01/*02`` and only ``*02`` is separable at the read length, every
    unambiguous observation is ``*02``, the likelihood says homozygous, and the restriction then
    empties every ``*01`` read -- systematically, on every sample with that genotype. So a genotype
    that fails to intersect ``fraction_to_explain`` of the compatible SETS is rejected as
    ``unexplained`` however well it fits the counts.
    """
    from itertools import combinations

    if len(alleles) == 1:
        return alleles, float("inf"), "single_allele"
    assigned = [c for c in clonotypes if c.allele is not None]
    if len(assigned) < min_clonotypes:
        return (), 0.0, "low_support"

    counts = Counter(c.allele for c in assigned)
    sets = [c.compatible for c in clonotypes]
    scored = []
    for size in range(1, min(ploidy, len(alleles)) + 1):
        for combo in combinations(alleles, size):
            ll = _log_likelihood(counts, combo, alleles, eps)
            # A total order on ties: likelihood, then fewer alleles, then the names themselves.
            scored.append((-ll, len(combo), combo))
    scored.sort()
    best, runner = scored[0], scored[1]
    log10_bf = (runner[0] - best[0]) / math.log(10)      # -(-ll_best) - -(-ll_runner), in log10
    if log10_bf < min_log10_bf:
        return (), log10_bf, "undetermined"
    chosen = best[2]
    explained = sum(1 for s in sets if s & set(chosen)) / len(sets) if sets else 0.0
    if explained < fraction_to_explain:
        return (), log10_bf, "unexplained"
    return chosen, log10_bf, "ok"


def infer_genotype(airr, *, organism: str = "human", loci=None,
                   min_log10_bf: float = MIN_LOG10_BF,
                   min_clonotypes: int = MIN_CLONOTYPES, ploidy: int = PLOIDY,
                   fraction_to_explain: float = FRACTION_TO_EXPLAIN,
                   scope: str = "framework", unmutated_only: bool = True):
    """``(rows, report)`` -- the V alleles this donor carries, one row per carried allele.

    ``airr`` is a Stage-1 AIRR TSV: the only artifact carrying ``v_germline_start``/``_end``, which
    the tie test needs. The clonotype table has no germline coordinates and cannot be used.

    ``scope="framework"`` clips each read's germline span at ``v_anchor_nt``. The V germline runs a
    median of 13 nt (range 8-41 over 792 human V alleles) past that anchor, into the junction --
    where chew-back and N/P addition live and the read is not expected to match germline at all.
    Including it makes tie sets artificially narrow, which manufactures confident observations.
    """
    from .annotate.airr_out import read_airr
    from .annotate.ties import TieResolver

    df = read_airr(airr)
    for col in ("v_call", "v_germline_start", "v_germline_end"):
        if col not in df.columns:
            raise ValueError(f"{airr}: a genotype needs {col}")
    if unmutated_only and "v_mutations" not in df.columns:
        # Never: raise rather than silently skipping the guard -- without it every sequencing error
        # at a discriminating base becomes a confident observation of the wrong allele.
        raise ValueError(
            f"{airr}: `--unmutated` needs the v_mutations column (run `map` with `--shm framework`,"
            " the default), or pass `--no-unmutated` and accept that reads vote on their called"
            " allele's germline regardless of what the read itself says")
    if loci:
        df = df.filter(pl.col("locus").is_in(list(loci)))

    germ = segment_germlines(organism, "v")
    clonotypes = _clonotypes(df, TieResolver(germ), unmutated_only=unmutated_only, scope=scope)
    eps = error_rate(clonotypes)

    by_gene: dict[str, list[Clonotype]] = {}
    for c in clonotypes:
        by_gene.setdefault(c.gene, []).append(c)

    universe: dict[str, tuple[str, ...]] = {}
    for allele in germ:
        universe[gene_of(allele)] = universe.get(gene_of(allele), ()) + (allele,)

    assigned = sum(1 for c in clonotypes if c.allele is not None)
    rows, report = [], {
        "genes": 0, "called": 0, "het": 0, "hom": 0,
        "clonotypes": len(clonotypes), "reads": sum(c.reads for c in clonotypes),
        "assigned": assigned,
        # Never: report how DISCRIMINATING the library is, not just what was called. A library whose
        # reads are too short to separate a gene's alleles cannot genotype it at any depth, and
        # this is the number that says so before anyone reads a call.
        "assigned_fraction": assigned / len(clonotypes) if clonotypes else 0.0,
        "error_rate": eps,
    }
    for gene in sorted(by_gene):
        group = by_gene[gene]
        alleles = tuple(sorted(universe.get(gene, ())))
        if not alleles:
            continue
        report["genes"] += 1
        chosen, log10_bf, note = _call_gene(
            alleles, group, eps, min_log10_bf=min_log10_bf, min_clonotypes=min_clonotypes,
            ploidy=ploidy, fraction_to_explain=fraction_to_explain)
        counts = Counter(c.allele for c in group if c.allele is not None)
        reads = Counter()
        for c in group:
            if c.allele is not None:
                reads[c.allele] += c.reads
        bf = "" if log10_bf == float("inf") else f"{log10_bf:.2f}"
        if chosen:
            report["called"] += 1
            report["het" if len(chosen) > 1 else "hom"] += 1
            for allele in chosen:
                rows.append((group[0].locus, gene, allele, str(counts.get(allele, 0)),
                             str(reads.get(allele, 0)), str(len(group)), bf, note))
        else:
            # Never: a gene that could not be called is REPORTED with its reason, not dropped.
            rows.append((group[0].locus, gene, "", "", "", str(len(group)), bf, note))
    return rows, report
