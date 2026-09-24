"""Recombination scenarios: the reconstruction invariant, and the degeneracy that breaks EM.

`arda.scenarios` estimates the generative model `arda.dpost` currently borrows from OLGA. Two
things can go wrong silently and both are pinned here:

1. **A scenario that does not reproduce the junction** is not a scenario. Every tuple the
   enumerator emits must reassemble the input exactly, or the counts are of something else.
2. **An insertion that costs only its length makes EM degenerate.** The junction's 5' end reads
   either as templated V or as an insertion that happens to match V germline; without a
   per-base cost the second is free, and EM walks into the corner where everything is insertion.

No mmseqs and no network: the human/mouse references are committed.
"""

from __future__ import annotations

import pytest

from arda import scenarios as sc

# A real human TRB junction, C...F, in frame.
TRB = ("TGCAGTGCTAGAGATTTAGCGGGAGGGACTGAAGCTTTCTTT", "TRBV20-1*01", "TRBJ2-1*01")
# TRA -- a VJ locus, so no D and no insDJ.
TRA = ("TGTGCTGTGAGAGACAGCAACTATCAGTTAATCTGG", "TRAV1-1*01", "TRAJ33*01")


def _reconstruct(s: sc.Scenario, g: sc.Germlines) -> str:
    """Rebuild the junction the scenario claims, from germline alone."""
    v = g.v_nt[: len(g.v_nt) - s.del_v]
    j = g.j_nt[s.del_j:]
    d = ""
    if s.d_call:
        dg = dict(g.d_germlines)[s.d_call]
        d = dg[s.del_dl: len(dg) - s.del_dr]
    return v, s.ins_vd, d, s.ins_dj, j


@pytest.mark.parametrize("junction,v,j", [TRB, TRA], ids=["TRB", "TRA"])
def test_every_scenario_reconstructs_the_junction_exactly(junction, v, j):
    """The invariant the whole module rests on."""
    g = sc.germlines_for(v, j, "human")
    found = sc.enumerate_scenarios(junction, v, j, "human")
    assert found, "no scenario at all"
    for s in found:
        vpart, ins1, d, ins2, jpart = _reconstruct(s, g)
        assert len(vpart) + ins1 + len(d) + ins2 + len(jpart) == len(junction)
        assert junction.startswith(vpart), s
        assert junction.endswith(jpart), s
        if d:
            at = len(vpart) + ins1
            assert junction[at: at + len(d)] == d, s


def test_a_vj_locus_has_no_d_and_no_ins_dj():
    """TRA/TRG/IGK/IGL are not skipped -- they carry most of the trimming evidence."""
    found = sc.enumerate_scenarios(*TRA, "human")
    assert found
    assert all(s.d_call == "" and s.ins_dj == 0 for s in found)

    stats = sc.estimate([(*TRA, 1.0)], organism="human", iterations=1)
    kinds = {k for _, k, _ in stats.counts}
    assert "insVD" in kinds and "delV" in kinds and "delJ" in kinds
    assert "insDJ" not in kinds and "dlen" not in kinds


def test_a_fully_trimmed_d_is_enumerated():
    """`dlen == 0` is the modal outcome, not a failure.

    Dropping the records where no D survived would truncate `dlen` at 1 and inflate every
    insertion distribution by the length of the D that was really there.
    """
    found = sc.enumerate_scenarios(*TRB, "human")
    assert any(s.d_call and s.del_dl + s.del_dr == len(dict(sc.germlines_for(*TRB[1:], "human")
                                                            .d_germlines)[s.d_call])
               for s in found)


def test_each_record_contributes_one_in_total():
    """A junction with hundreds of admissible readings must not outvote one with three."""
    stats = sc.estimate([(*TRB, 1.0)], organism="human", iterations=1)
    # `d_marginal` is touched exactly once per scenario, so its total is the record's total mass.
    total = sum(v for (_, kind, _), v in stats.counts.items() if kind == "d_marginal")
    assert total == pytest.approx(1.0, rel=1e-9)

    # And the weight scales it linearly.
    heavy = sc.estimate([(*TRB, 7.0)], organism="human", iterations=1)
    assert sum(v for (_, kind, _), v in heavy.counts.items()
               if kind == "d_marginal") == pytest.approx(7.0, rel=1e-9)


def test_insertion_base_cost_is_load_bearing(monkeypatch):
    """Without `0.25^len`, EM inflates insertions -- the degeneracy the term exists to break.

    Measured on 503 real human TRB junctions before the term existed: three iterations moved
    `insVD` mass onto 10-11 nt with the log-likelihood rising the whole way. Here, one junction
    and one iteration are enough to show the direction.
    """
    def mean_ins(records):
        stats = sc.estimate(records, organism="human", iterations=2)
        rows = [(int(k), v) for _, kind, k, v in stats.rows() if kind == "insVD"]
        return sum(k * v for k, v in rows) / sum(v for _, v in rows)

    records = [(*TRB, 1.0)]
    with_cost = mean_ins(records)
    monkeypatch.setattr(sc, "_P_BASE", 1.0)          # an insertion costs only its length
    without_cost = mean_ins(records)
    assert without_cost > with_cost + 1.0, (
        f"free insertions must pull the mean up: {without_cost:.2f} vs {with_cost:.2f}")


def test_output_is_a_drop_in_for_the_shipped_prior(tmp_path, monkeypatch):
    """The point of the module: `dpost.load_d_prior` must read what this writes.

    Same long `locus/kind/key/value` shape, same key grammar (`<allele>:<n>` for `dlen`,
    `<d>|<j>` for `d_given_j`). If this drifts, the estimate is unusable by the one consumer it
    exists for.
    """
    from arda import dpost

    stats = sc.estimate([(*TRB, 1.0)], organism="human", iterations=2)
    rows = [r for r in stats.rows() if r[1] in ("insVD", "insDJ", "dlen",
                                                "d_marginal", "d_given_j")]
    assert rows
    out = tmp_path / "human"
    out.mkdir()
    with open(out / "d_prior.tsv", "w") as fh:
        fh.write("\t".join(sc.PRIOR_COLUMNS) + "\n")
        for locus, kind, key, value in rows:
            fh.write(f"{locus}\t{kind}\t{key}\t{value:.8g}\n")

    monkeypatch.setattr(dpost, "vdj_dir", lambda *a, **k: out)
    prior = dpost.load_d_prior("human")
    assert "TRB" in prior
    p = prior["TRB"]
    assert p.ins_vd and sum(p.ins_vd) == pytest.approx(1.0, rel=1e-3)
    assert p.dlen, "dlen did not parse back into per-allele curves"
    assert any(a.startswith("TRBD") for a in p.dlen)


def test_unresolvable_calls_are_skipped_not_guessed():
    stats = sc.estimate([("TGCAGTGCTAGAGATTTAGCGGG", "NOSUCHV*01", "NOSUCHJ*01", 1.0)],
                        organism="human", iterations=1)
    assert stats.records == 0 and stats.skipped == 1
    assert not stats.counts
