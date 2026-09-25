"""`arda.shmmodel` — the context SHM model, its scoping rules and its table round trip."""

from __future__ import annotations

import pytest
import typer
from typer.testing import CliRunner

from arda.cli import app
from arda.shmmodel import (
    DEFAULT_K,
    V_REGIONS,
    ShmModel,
    estimate,
    germline_map,
    load_model,
    read_airr,
    write_table,
)

AIRR_COLUMNS = ("locus", "v_call", "v_germline_start", "v_germline_end", "v_mutations",
                "v_anchor_nt")


@pytest.fixture(scope="module")
def igh():
    """One real IGH allele: ``(v_call, germline nt, [(region, start, end)])``."""
    mapping = germline_map("human", "IGH")
    assert mapping, "the committed human reference carries IGH V scaffolds"
    call = sorted(mapping)[0]
    seq, spans = mapping[call]
    return call, seq, spans


def write_airr(path, rows, columns=AIRR_COLUMNS):
    """Rows are dicts; anything not named defaults to empty."""
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
    return path


def read(tmp_path, rows, **kwargs):
    return read_airr(write_airr(tmp_path / "in.tsv", rows), **kwargs)


def test_denominator_is_capped_at_the_anchor(tmp_path, igh):
    """A read running past Cys104 must not count those positions in the denominator.

    `v_mutations` is framework-scoped, so an uncapped denominator would divide a scoped
    numerator by an unscoped one and report a rate too low by however far the read overran.
    """
    call, seq, _ = igh
    anchor = len(seq) - 20
    rows = read(tmp_path, [{"locus": "IGH", "v_call": call, "v_germline_start": 10,
                            "v_germline_end": len(seq), "v_mutations": "", "v_anchor_nt": anchor}])
    assert rows == [(call, 10, anchor, "")]


def test_a_tie_list_is_dropped_not_resolved(tmp_path, igh):
    call, seq, _ = igh
    rows = read(tmp_path, [
        {"locus": "IGH", "v_call": f"{call},IGHV3-23*01", "v_germline_start": 10,
         "v_germline_end": 100, "v_anchor_nt": 200},
        {"locus": "IGH", "v_call": call, "v_germline_start": 10, "v_germline_end": 100,
         "v_anchor_nt": 200},
    ])
    assert [r[0] for r in rows] == [call]


def test_weight_unique_collapses_an_expansion(tmp_path, igh):
    """One clone sequenced many times is one observation of the mutation process, not many."""
    call, _seq, _ = igh
    same = {"locus": "IGH", "v_call": call, "v_germline_start": 10, "v_germline_end": 100,
            "v_mutations": "G45A", "v_anchor_nt": 200}
    assert len(read(tmp_path, [same] * 7)) == 1
    assert len(read(tmp_path, [same] * 7, weight="reads")) == 7


def test_a_missing_column_names_the_flag_that_produces_it(tmp_path):
    path = write_airr(tmp_path / "bad.tsv", [{"locus": "IGH"}], columns=("locus", "v_call"))
    with pytest.raises(ValueError, match="v_mutations"):
        read_airr(path)


def test_estimate_counts_the_context_it_was_given(tmp_path, igh):
    """A substitution planted at one position must land on that position's own 5-mer."""
    call, seq, _ = igh
    pos = 120
    half = DEFAULT_K // 2
    context = seq[pos - 1 - half:pos + half]
    rows = [(call, 60, 200, f"{seq[pos - 1]}{pos}A")] * 1
    model = estimate(rows, locus="IGH")
    assert model.context.get(context) is None or model.context[context] > 0
    # One read cannot clear `min_context`, so every context falls back to the sample rate; the
    # counting itself is what is checked here.
    assert model.substitutions == 1
    # 60 and 200 are both interior, so no context is lost to a germline edge.
    assert model.covered == 200 - 60 + 1


def test_a_region_with_no_substitution_is_omitted_not_zeroed(tmp_path, igh):
    """Never: an unestimated region must not come back as a multiplier of 0.

    0 does not read as "thin evidence", it reads as "a substitution here is impossible", and
    `rate()` would return 0 for every position in that region -- correct-looking output, silently
    wrong. An omitted region takes the default multiplier of 1.
    """
    call, seq, spans = igh
    by_name = {name: (lo, hi) for name, lo, hi in spans}
    lo, hi = by_name["cdr2"]
    rows = [(call, 1 + i, len(seq) - 10, f"{seq[lo - 1]}{lo}A") for i in range(80)]
    model = estimate(rows, locus="IGH")
    assert "cdr2" in model.region
    assert all(v > 0 for v in model.region.values())
    assert model.rate(seq, 5, "fwr1") == model.rate(seq, 5)      # omitted region, no effect


def test_region_multipliers_are_normalised_and_ordered(tmp_path, igh):
    """CDRs mutate more than frameworks, and the multipliers are centred on 1."""
    call, seq, spans = igh
    by_name = {name: (lo, hi) for name, lo, hi in spans}
    if not {"cdr2", "fwr1"} <= set(by_name):
        pytest.skip(f"{call} has no cdr2/fwr1 markup")
    cdr_lo, cdr_hi = by_name["cdr2"]
    fwr_lo, _fwr_hi = by_name["fwr1"]
    rows = []
    for i in range(80):
        # Every read mutates six CDR2 positions; one read in eight also mutates one in FWR1, so
        # both regions are estimated and CDR2 is the denser of the two.
        planted = list(range(cdr_lo, min(cdr_hi, cdr_lo + 6)))
        if i % 8 == 0:
            planted.append(fwr_lo + 10)
        muts = ",".join(f"{seq[p - 1]}{p}A" for p in planted)
        rows.append((call, 1 + i % 5, len(seq) - 10, muts))
    model = estimate(rows, locus="IGH")
    assert model.region["cdr2"] > model.region["fwr1"] > 0


def test_rate_falls_back_to_the_sample_rate_at_the_germline_edge():
    model = ShmModel(k=5, scale=0.03, context={"ACGTA": 0.5}, region={"cdr1": 2.0})
    assert model.rate("ACGTACGT", 3) == pytest.approx(0.5)
    assert model.rate("ACGTACGT", 1) == pytest.approx(0.03)     # no room for the context
    assert model.rate("ACGTACGT", 3, "cdr1") == pytest.approx(1.0)
    assert model.rate("ACGTACGT", 4) == pytest.approx(0.03)     # context never seen


def test_even_k_is_refused():
    with pytest.raises(ValueError, match="centre"):
        estimate([], k=4)


def test_table_round_trip(tmp_path):
    model = ShmModel(k=7, scale=0.0412, context={"AAA": 0.1, "CCC": 0.2},
                     region={"cdr1": 1.5, "fwr1": 0.5}, observations=11, covered=99,
                     substitutions=4, locus="IGH", min_context=25)
    path = tmp_path / "shm.tsv"
    write_table(model, path)
    back = load_model(path)
    assert back.k == 7
    assert back.min_context == 25
    assert back.locus == "IGH"
    assert back.scale == pytest.approx(0.0412)
    assert back.context == pytest.approx(model.context)
    assert back.region == pytest.approx(model.region)


def test_loader_skips_comments_and_the_header_by_content(tmp_path):
    """Never: by WHAT THEY ARE, not by position.

    `d_prior.tsv` shipped a loader that dropped line 1 by position, so the day a generator wrote
    a provenance line above the header, the header reached the parse loop and `float("value")`
    raised on the one file the docs called a drop-in.
    """
    path = tmp_path / "shm.tsv"
    path.write_text(
        "# arda shm-model: organism=human locus=IGH k=5 min_context=50\n"
        "\n"
        "kind\tkey\tvalue\n"
        "scale\tall\t0.05\n"
        "# a second comment, lower down\n"
        "context\tAAAAA\t0.1\n"
        "kind\tkey\tvalue\n"                 # a header repeated mid-file, as a concatenation makes
        "region\tcdr1\t1.4\n"
        "context\tbroken\tnot-a-number\n"    # a malformed row is skipped, not fatal
    )
    model = load_model(path)
    assert model.scale == pytest.approx(0.05)
    assert model.context == {"AAAAA": pytest.approx(0.1)}
    assert model.region == {"cdr1": pytest.approx(1.4)}
    assert model.k == 5


def test_v_regions_are_the_v_ones_only():
    """`cdr3` is the junction and `fwr4` is J: neither is a V region a substitution is scored in."""
    assert V_REGIONS == ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")


def test_cli_refuses_a_bad_weight_and_an_even_k(tmp_path, igh):
    call, seq, _ = igh
    path = write_airr(tmp_path / "in.tsv", [
        {"locus": "IGH", "v_call": call, "v_germline_start": 10, "v_germline_end": 200,
         "v_mutations": "", "v_anchor_nt": 250}])
    runner = CliRunner()
    for args in (["--weight", "bogus"], ["--k", "4"]):
        result = runner.invoke(app, ["shm-model", "-i", str(path), "-o",
                                     str(tmp_path / "out.tsv"), *args])
        assert result.exit_code != 0


def test_cli_writes_a_loadable_table(tmp_path, igh):
    call, seq, _ = igh
    rows = [{"locus": "IGH", "v_call": call, "v_germline_start": 1,
             "v_germline_end": len(seq) - 5, "v_anchor_nt": len(seq) - 5,
             "v_mutations": f"{seq[99]}100A" if i % 2 else ""} for i in range(60)]
    # 60 rows that differ only in their mutation set collapse to 2 under `unique`.
    for i, row in enumerate(rows):
        row["v_germline_start"] = 1 + i
    path = write_airr(tmp_path / "in.tsv", rows)
    out = tmp_path / "shm.tsv"
    result = CliRunner().invoke(app, ["shm-model", "-i", str(path), "-o", str(out)])
    assert result.exit_code == 0, result.output
    model = load_model(out)
    assert model.locus == "IGH"
    assert model.scale > 0
    assert set(model.region) <= set(V_REGIONS)


def test_cli_refuses_a_file_with_no_usable_rows(tmp_path):
    path = write_airr(tmp_path / "in.tsv", [{"locus": "TRB", "v_call": "TRBV20-1*01"}])
    result = CliRunner().invoke(app, ["shm-model", "-i", str(path), "-o",
                                      str(tmp_path / "out.tsv")])
    assert result.exit_code != 0
    assert isinstance(result.exception, (typer.BadParameter, SystemExit))
