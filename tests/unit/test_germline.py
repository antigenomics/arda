"""Per-allele germlines recovered from the committed reference.

The claim under test is narrow and load-bearing: ``scaffold[:v_sequence_end]`` is the V allele's
own germline, in the SAME coordinate space as the ``v_germline_start``/``v_germline_end`` columns
arda emits. Everything that consumes a germline by allele name -- tie lists today, genotype
restriction next -- rests on it.
"""

from __future__ import annotations

import polars as pl
import pytest

from arda.germline import SEGMENTS, segment_germlines


def _reference_available() -> bool:
    from arda.paths import vdj_dir
    try:
        base = vdj_dir("human")
    except Exception:
        return False
    return (base / "markup.tsv").exists() and (base / "alleles.fasta").exists()


needs_reference = pytest.mark.skipif(not _reference_available(), reason="no committed reference")


# --- argument handling --------------------------------------------------------------------------

def test_an_unknown_segment_is_refused():
    with pytest.raises(ValueError, match="unknown segment"):
        segment_germlines("human", "d")


def test_a_missing_reference_raises_rather_than_returning_nothing(tmp_path, monkeypatch):
    """Never: an empty dict would leave every call untouched and report success."""
    import arda.germline as g
    monkeypatch.setattr(g, "vdj_dir", lambda org: tmp_path / org)
    g.segment_germlines.cache_clear()
    with pytest.raises(FileNotFoundError, match="no reference"):
        segment_germlines("nowhere", "v")
    g.segment_germlines.cache_clear()


# --- against the real reference -----------------------------------------------------------------

@needs_reference
def test_every_v_germline_is_the_prefix_of_a_scaffold_that_calls_it():
    """The derivation itself, checked against the artifacts rather than restated."""
    from arda.paths import vdj_dir
    from arda.refbuild.imgt import read_fasta

    base = vdj_dir("human")
    seqs = dict(read_fasta(base / "alleles.fasta"))
    markup = pl.read_csv(base / "markup.tsv", separator="\t", infer_schema_length=0)
    germ = segment_germlines("human", "v")

    checked = 0
    for row in markup.iter_rows(named=True):
        allele, end = (row["v_call"] or "").strip(), row["v_sequence_end"]
        scaffold = seqs.get(row["scaffold_id"])
        if not allele or scaffold is None or not str(end).isdigit() or int(end) <= 0:
            continue
        for member in (a.strip() for a in allele.split(",")):
            assert scaffold.startswith(germ[member][:200]), member
        checked += 1
    assert checked > 10_000, "the human reference should have ~15k V-bearing scaffolds"


@needs_reference
def test_the_human_v_and_j_sets_are_the_ones_the_reference_can_call():
    """⚠ Assert the COUNT, not the flag: a reference swap can silently be a no-op.

    801 V / 132 J is the built human reference, and it is deliberately NOT the 884 functional V
    alleles IMGT ships -- ``drop_unanchorable`` / ``drop_truncated`` and incomplete IgBLAST markup
    remove the rest. Offering one of those in a tie list would name an allele that could never
    have been the call.
    """
    v, j = segment_germlines("human", "v"), segment_germlines("human", "j")
    assert len(v) == 801 and len(j) == 132
    assert min(map(len, v.values())) >= 200 and max(map(len, v.values())) <= 400
    assert min(map(len, j.values())) >= 30 and max(map(len, j.values())) <= 100
    assert all(set(s) <= set("ACGTN") for s in v.values())


@needs_reference
def test_grouped_v_calls_are_split_into_their_member_alleles():
    """Never: ``markup.tsv``'s ``v_call`` is not always ONE allele.

    Scaffolds are deduplicated by assembled sequence, so alleles that produce byte-identical
    scaffolds collapse into one row whose ``v_call`` is a comma-joined group. Human has 23 such
    groups hiding 49 alleles, and J has 3 more groups. Keyed on the group STRING they are
    invisible to every consumer: ``TieResolver.expand`` takes ``call.split(",")[0]``, looks it
    up, misses, and returns the call untouched -- a tie list that silently never fires.
    """
    from arda.annotate.ties import TieResolver

    v = segment_germlines("human", "v")
    assert not [k for k in v if "," in k]
    assert not [k for k in segment_germlines("human", "j") if "," in k]
    assert "IGKV1-33*01" in v and "IGKV1D-33*01" in v
    assert TieResolver(v).expand("IGKV1-33*01", 1, 200) == "IGKV1-33*01,IGKV1D-33*01"


@needs_reference
def test_some_v_alleles_are_indistinguishable_at_any_span():
    """The limit of what a genotype can ever claim, pinned so nobody assumes otherwise.

    49 human V alleles in 23 groups have byte-identical germlines, so no read length separates
    them -- allele ambiguity is USUALLY a coverage problem, but not always. Two of the groups are
    TR (``TRBV6-2*01``/``TRBV6-3*01`` and ``TRBV24-1*02``/``TRBV24/OR9-2*01``), which is why a
    genotype caller has to carry such a group whole rather than pick a member.
    """
    from collections import defaultdict

    by_seq: dict[str, list[str]] = defaultdict(list)
    for allele, seq in segment_germlines("human", "v").items():
        by_seq[seq].append(allele)
    groups = sorted(sorted(g) for g in by_seq.values() if len(g) > 1)
    assert sum(len(g) for g in groups) == 49
    assert ["TRBV6-2*01", "TRBV6-3*01"] in groups


@needs_reference
@pytest.mark.parametrize("organism", ["human", "mouse", "rat", "rabbit", "rhesus_monkey"])
@pytest.mark.parametrize("segment", SEGMENTS)
def test_every_shipped_organism_yields_germlines(organism, segment):
    from arda.paths import vdj_dir
    if not (vdj_dir(organism) / "markup.tsv").exists():
        pytest.skip(f"{organism} not built")
    got = segment_germlines(organism, segment)
    assert got and all(got.values())


# --- the wobble rule ----------------------------------------------------------------------------

def _fake_reference(tmp_path, rows, seqs):
    d = tmp_path / "vdj"
    d.mkdir(parents=True, exist_ok=True)
    (d / "alleles.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    pl.DataFrame(rows).write_csv(d / "markup.tsv", separator="\t", quote_style="never")
    return d


@pytest.fixture
def fake(tmp_path, monkeypatch):
    import arda.germline as g

    def _install(rows, seqs):
        d = _fake_reference(tmp_path, rows, seqs)
        monkeypatch.setattr(g, "vdj_dir", lambda org: d)
        g.segment_germlines.cache_clear()
        return d

    yield _install
    g.segment_germlines.cache_clear()


def test_an_allele_whose_markup_wobbles_takes_the_mode(fake):
    """Never: ``v_sequence_end`` is IgBLAST's markup of ONE assembled scaffold, not a property of
    the allele, and it wobbles -- mouse disagrees on 4 of 897 V alleles, always by 1-2 nt and
    always with a landslide majority (``TRAV16*02``: 288 nt on 58 scaffolds, 290 nt on 1).

    First-wins or last-wins would hand the whole reference's tie resolution to whichever scaffold
    happened to sort first or last.
    """
    seq = "ACGT" * 30                                   # 120 nt
    fake(
        rows=[{"scaffold_id": f"L_{i}", "v_call": "V*01",
               "v_sequence_end": 100 if i < 5 else 98} for i in range(6)],
        seqs={f"L_{i}": seq for i in range(6)},
    )
    assert segment_germlines("x", "v")["V*01"] == seq[:100]


def test_a_dead_heat_breaks_by_longer_then_by_sequence(fake):
    """A TOTAL order -- two runs over one reference cannot disagree."""
    seq = "ACGT" * 30
    fake(
        rows=[{"scaffold_id": "L_0", "v_call": "V*01", "v_sequence_end": 90},
              {"scaffold_id": "L_1", "v_call": "V*01", "v_sequence_end": 100}],
        seqs={"L_0": seq, "L_1": seq},
    )
    assert segment_germlines("x", "v")["V*01"] == seq[:100]


def test_blank_cells_and_unusable_rows_are_skipped(fake):
    """``markup.tsv`` writes a blank as the two characters ``""``; J+C scaffolds carry no v_call."""
    seq = "ACGT" * 30
    fake(
        rows=[{"scaffold_id": "L_0", "v_call": "V*01", "v_sequence_end": 100},
              {"scaffold_id": "L_1", "v_call": "", "v_sequence_end": ""},
              {"scaffold_id": "L_2", "v_call": "V*02", "v_sequence_end": 0},
              {"scaffold_id": "absent", "v_call": "V*03", "v_sequence_end": 50}],
        seqs={"L_0": seq, "L_1": seq, "L_2": seq},
    )
    assert segment_germlines("x", "v") == {"V*01": seq[:100]}


def test_j_is_sliced_from_its_own_start(fake):
    seq = "ACGT" * 30
    fake(
        rows=[{"scaffold_id": "L_0", "v_call": "V*01", "v_sequence_end": 60,
               "j_call": "J*01", "j_sequence_start": 71, "vj_end": 110}],
        seqs={"L_0": seq},
    )
    assert segment_germlines("x", "j")["J*01"] == seq[70:110]
