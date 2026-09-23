"""The FIRST run against a fresh reference must not be the odd one out.

`segments.fasta` + `segments.markup.tsv` are GENERATED, not shipped, and they are generated on the
map path — after `load_reference` has already read `entries`. So the run that generated them used
to carry a reference that knew nothing about them: `segment_j_call` handed raw target names to the
combination lookup, `_cannot_reach_cys104` found no CDR3 markup, and `--v-only-on-segment` routed
nothing. Measured on `tests/data/rnaseq_real`, `arda amplicon` with both segment artifacts
removed: the generating run differed from the next run on **268 of 453** Stage-1 AIRR rows, read
`fast_fraction` **0.08** against **0.1736**, and omitted the `v_only_on_segment` counter (283) from
its report entirely. Exit 0 throughout — which is why it survived a release: the only symptom is
that run #1 and run #2 disagree, and nobody runs the same sample twice.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from arda.annotate import mapper
from arda.annotate.reference import Reference, _load_markup
from arda.paths import vdj_dir

from ..conftest import requires_human_db, requires_mmseqs


def _bare_reference(tmp_path: Path) -> Reference:
    """The real human reference, minus the generated segment artifacts."""
    src = vdj_dir("human")
    base = tmp_path / "vdj" / "human"
    base.mkdir(parents=True)
    for name in ("alleles.fasta", "markup.tsv"):
        (base / name).symlink_to(src / name)
    entries: dict = {}
    _load_markup(base / "markup.tsv", entries)
    return Reference("human", "nt", base / "alleles.fasta", entries, {}, {})


@requires_human_db
@requires_mmseqs
def test_the_run_that_generates_the_segment_reference_also_loads_its_markup(tmp_path, monkeypatch):
    monkeypatch.setattr(mapper, "data_dir", lambda: tmp_path / "data")
    ref = _bare_reference(tmp_path)
    assert not any("|" in k for k in ref.entries), "fixture is not a fresh reference"

    db = mapper._cached_segment_db(ref, "human")

    assert db is not None, "the segment reference must be generated, not tolerated"
    assert (ref.target_fasta.parent / "segments.markup.tsv").exists()
    v = [k for k in ref.entries if k.startswith("V|")]
    j = [k for k in ref.entries if k.startswith("J|")]
    assert v and j, "the generated segment markup was not loaded into the live reference"
    # The two things the mapper actually reads off a segment entry. `--v-only-on-segment` gates
    # on the V segment's own CDR3 start (`_cannot_reach_cys104` bails when it is <= 0, which is
    # exactly what a missing entry produced), and the J->C contest resolves through `c_call`.
    from arda.annotate.reference import REGIONS
    assert ref.get(v[0]).starts[REGIONS.index("cdr3")] > 0, \
        "no CDR3 markup on a V segment: --v-only-on-segment silently routes nothing"
    assert any(ref.get(k).c_call for k in ref.entries if k.startswith("C|")), \
        "a C| target must resolve to a constant allele, or every J->C read is left unpaired"


@requires_human_db
@requires_mmseqs
def test_a_surviving_fasta_with_no_markup_sibling_still_loads(tmp_path, monkeypatch):
    """The generated pair can be broken up — a `written` flag would miss exactly this case."""
    monkeypatch.setattr(mapper, "data_dir", lambda: tmp_path / "data")
    ref = _bare_reference(tmp_path)
    from arda.refbuild.segments import build_segment_reference
    build_segment_reference("human", out_dir=ref.target_fasta.parent)

    mapper._cached_segment_db(ref, "human")

    assert any(k.startswith("V|") for k in ref.entries)


@requires_human_db
def test_a_stale_pre_2_8_0_segment_key_does_not_outlive_the_file_it_came_from(tmp_path):
    """Regeneration renames `JC|<scaffold>` to `C|<allele>`; the old keys must go with it."""
    ref = _bare_reference(tmp_path)
    ref.entries["JC|TRB_stale"] = next(iter(ref.entries.values()))
    base = ref.target_fasta.parent
    from arda.refbuild.segments import build_segment_reference
    build_segment_reference("human", out_dir=base)

    ref.load_segment_markup()

    assert "JC|TRB_stale" not in ref.entries
    assert any(k.startswith("C|") for k in ref.entries)


@requires_human_db
def test_reloading_keeps_every_base_scaffold(tmp_path):
    """Dropping the `|` key space must not touch the V x J scaffolds it lives beside."""
    ref = _bare_reference(tmp_path)
    before = {k: v for k, v in ref.entries.items() if "|" not in k}
    from arda.refbuild.segments import build_segment_reference
    build_segment_reference("human", out_dir=ref.target_fasta.parent)

    ref.load_segment_markup()

    assert {k: v for k, v in ref.entries.items() if "|" not in k} == before
    assert ref._jc_combos is None, "the cached (j, c) -> scaffold map is derived and must be dropped"


def test_load_segment_markup_is_a_no_op_without_the_file(tmp_path):
    ref = Reference("human", "nt", tmp_path / "alleles.fasta", {}, {}, {})
    ref.load_segment_markup()
    assert ref.entries == {}


@pytest.mark.parametrize("seqtype", ["aa"])
def test_the_amino_acid_reference_has_no_segment_markup(tmp_path, seqtype):
    """`load_reference` only loads it for nt; the reload must agree, not invent a second rule."""
    (tmp_path / "segments.markup.tsv").write_text("scaffold_id\n")
    ref = Reference("human", seqtype, tmp_path / "alleles.aa.fasta", {}, {}, {})
    ref.load_segment_markup()
    assert ref.entries == {}
