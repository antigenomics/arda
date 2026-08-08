"""`arda export-ref` — the reference is only useful exported if its coordinates still cut correctly.

The reference is arda's most valuable offline artifact, and hand-joining its TSVs against its
FASTAs is exactly the kind of thing that goes wrong quietly: coordinates are 1-based closed, a
`J + C` scaffold has no V at all, and the aa reference has three frames per D allele. A join that
is off by one produces plausible nonsense — a browser would happily draw the shifted feature.

So these tests do not check that the exporter *ran*; they check that what it wrote **round-trips**:
every region's reported `*_seq` must equal the slice its own `*_start`/`*_end` imply, and the
declared junction must equal the slice the CDR3/FR4 coordinates imply.
"""

from __future__ import annotations

import csv
import io

import pytest

from arda.annotate.reference import REGIONS
from arda.refexport import FORMATS, KINDS, export_reference
from tests.conftest import requires_human_db

pytestmark = [requires_human_db]


def _tsv(tmp_path, **kw):
    out = tmp_path / "ref.tsv"
    n = export_reference("human", out=out, **kw)
    with open(out) as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    return n, rows


def test_every_exported_region_slice_matches_its_own_coordinates(tmp_path):
    """The round-trip that makes the export trustworthy: `*_seq` == sequence[start-1:end]."""
    n, rows = _tsv(tmp_path, kind="scaffolds", fmt="tsv", loci={"TRB"})
    assert n == len(rows) > 1000
    checked = 0
    for r in rows:
        seq = r["sequence"]
        assert len(seq) == int(r["sequence_length"])
        for reg in REGIONS:
            s, e = r[f"{reg}_start"], r[f"{reg}_end"]
            if not s or not e:
                assert r[f"{reg}_seq"] == "", f"{r['sequence_id']}: {reg} has no coords but a seq"
                continue
            assert r[f"{reg}_seq"] == seq[int(s) - 1:int(e)], f"{r['sequence_id']}: {reg} mismatch"
            checked += 1
    assert checked > 5000, f"only {checked} regions checked — the export looks empty"


def test_the_declared_junction_equals_what_the_cdr3_and_fwr4_coordinates_imply(tmp_path):
    """`junction` is Cys104..[FW]118, i.e. cdr3_start-3 .. fwr4_start+2 on the scaffold. If the
    export's coordinates and its junction column disagree, one of them is wrong."""
    _, rows = _tsv(tmp_path, kind="scaffolds", fmt="tsv", loci={"TRB"})
    checked = mismatched = 0
    for r in rows:
        if not (r["junction"] and r["cdr3_start"] and r["fwr4_start"]):
            continue
        js, je = int(r["cdr3_start"]) - 3, int(r["fwr4_start"]) + 2
        if js < 1 or je > len(r["sequence"]):
            continue
        checked += 1
        mismatched += r["sequence"][js - 1:je] != r["junction"]
    assert checked > 1000
    assert mismatched == 0, f"{mismatched} of {checked} junctions disagree with their coordinates"


def test_regions_are_ordered_and_do_not_overlap(tmp_path):
    _, rows = _tsv(tmp_path, kind="scaffolds", fmt="tsv", loci={"TRB"})
    bad = []
    for r in rows:
        prev_end = 0
        for reg in REGIONS:
            s, e = r[f"{reg}_start"], r[f"{reg}_end"]
            if not s:
                continue
            s, e = int(s), int(e)
            if s <= prev_end or e < s:
                bad.append((r["sequence_id"], reg, s, e, prev_end))
            prev_end = e
    assert not bad, f"{len(bad)} out-of-order/overlapping regions, e.g. {bad[:3]}"


@pytest.mark.parametrize("fmt", FORMATS)
def test_every_format_writes_something_well_formed(tmp_path, fmt):
    out = tmp_path / f"ref.{fmt}"
    n = export_reference("human", kind="scaffolds", fmt=fmt, loci={"TRB"}, out=out)
    assert n > 1000
    text = out.read_text()
    assert text
    if fmt == "fasta":
        assert text.startswith(">") and text.count(">") == n
    elif fmt == "gff3":
        assert text.startswith("##gff-version 3")
        # GFF3 is 1-based closed like arda, so a feature must never start at 0.
        starts = [int(l.split("\t")[3]) for l in text.splitlines()
                  if l and not l.startswith("#") and len(l.split("\t")) > 4]
        assert starts and min(starts) >= 1
    elif fmt == "airr":
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        assert len(rows) == n
        assert all(r["sequence_id"] and r["locus"] == "TRB" for r in rows)
        # A scaffold IS its own alignment, which is what makes an exported row valid AIRR input.
        assert all(r["sequence_alignment"] == r["sequence"] for r in rows)


def test_locus_filter_actually_filters(tmp_path):
    n_trb, rows = _tsv(tmp_path, kind="scaffolds", fmt="tsv", loci={"TRB"})
    assert {r["locus"] for r in rows} == {"TRB"}
    n_all = export_reference("human", kind="scaffolds", fmt="fasta",
                             out=tmp_path / "all.fa")
    assert n_all > n_trb


def test_segments_export_carries_the_per_segment_markup(tmp_path):
    n, rows = _tsv(tmp_path, kind="segments", fmt="tsv")
    assert n == len(rows) > 500
    kinds = {r["sequence_id"].split("|", 1)[0] for r in rows}
    assert kinds <= {"V", "J", "C", "JC"}, kinds
    v = [r for r in rows if r["sequence_id"].startswith("V|")]
    assert v, "no V segments exported"
    # A V segment must carry FR1-FR3; that is the property the two-pass path relies on.
    assert all(r["fwr1_start"] and r["fwr3_end"] for r in v)


def test_anchors_export_is_tsv_only_and_refuses_the_others(tmp_path):
    n = export_reference("human", kind="anchors", fmt="tsv", loci={"TRB"},
                         out=tmp_path / "anchors.tsv")
    assert n > 100
    with pytest.raises(ValueError, match="tsv"):
        export_reference("human", kind="anchors", fmt="fasta", out=tmp_path / "x.fa")


def test_bad_format_and_kind_are_rejected_by_name(tmp_path):
    with pytest.raises(ValueError, match="format"):
        export_reference("human", fmt="bed", out=tmp_path / "x")
    with pytest.raises(ValueError, match="kind"):
        export_reference("human", kind="germlines", out=tmp_path / "x")
    assert set(FORMATS) == {"tsv", "fasta", "gff3", "airr"}
    assert set(KINDS) == {"scaffolds", "segments", "anchors"}
