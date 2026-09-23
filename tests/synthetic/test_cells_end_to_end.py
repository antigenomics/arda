"""`arda cells` on a synthetic droplet library built from the committed receptor fixture.

The premise the whole module rests on is that a cell's molecules TILE its transcript even though
each molecule covers one window of it, so this test builds exactly that shape -- real receptor
sequences cut into 90 nt windows at random starts, named as migec molecules of one cell -- and
asks for the full-length contig and its junction back. Two cells are given two receptors each,
which is what a doublet is.
"""

from __future__ import annotations

import csv
import json
import random
from pathlib import Path

import pytest

from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "genbank_receptors.fa"
WINDOW = 90
DOUBLETS = 2


def _read_fasta(path: Path) -> list[str]:
    seqs, buf = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if buf:
                seqs.append("".join(buf))
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if buf:
        seqs.append("".join(buf))
    return seqs


@pytest.fixture(scope="module")
def droplet_library(tmp_path_factory):
    """One cell per receptor, plus `DOUBLETS` cells carrying two receptors of the SAME locus.

    Every molecule is a 90 nt window at a random start with a random read depth, which is the
    co-terminal pile a real UMI consensus is. Starts are uniform, so coverage is a coupon
    collector and the redundancy has to be generous before the tiling closes on the longest
    template -- at 3x the 1,620 nt receptor came back 1,499 nt. A real cell carries hundreds of
    molecules, so 6x is the conservative shape, not a favourable one.

    Never: **a doublet is two receptors of the same LOCUS.** A cell holding one TRA and one TRB
    is a normal paired T cell, and building the "doublets" from consecutive fixture records made
    exactly that -- the test then asked the pipeline to flag a healthy cell and it rightly did
    not.
    """
    from arda.annotate.contig import reannotate_contigs

    rng = random.Random(11)
    receptors = [s for s in _read_fasta(FIXTURE) if len(s) >= 300]
    assert len(receptors) >= 10, "fixture no longer holds enough long receptors"

    # What arda calls on the WHOLE template is the ceiling this test scores against. One of the
    # committed receptors gets no junction even from the intact sequence, so asserting "every
    # cell has a junction" would be a test of the fixture, not of the assembly.
    expected = dict(enumerate(reannotate_contigs(
        [(f"r{i}", s) for i, s in enumerate(receptors)], "human", threads=4)))

    by_locus: dict[str, list[int]] = {}
    for i, record in expected.items():
        if record.get("junction_aa"):
            by_locus.setdefault(record.get("locus") or "", []).append(i)
    same_locus = max(by_locus.values(), key=len)
    assert len(same_locus) >= 2 * DOUBLETS, "fixture no longer holds two receptors of one locus"

    cells: dict[str, list[int]] = {}
    for i in range(len(receptors)):
        cells["".join(rng.choice("ACGT") for _ in range(16))] = [i]
    for i in range(DOUBLETS):
        cells["".join(rng.choice("ACGT") for _ in range(16))] = [same_locus[2 * i],
                                                                 same_locus[2 * i + 1]]

    out = tmp_path_factory.mktemp("cells") / "consensus.fq"
    with open(out, "w") as fh:
        for cell, members in cells.items():
            for which in members:
                template = receptors[which]
                for _ in range(6 * (len(template) // 30)):
                    start = rng.randrange(0, len(template) - WINDOW + 1)
                    seq = template[start:start + WINDOW]
                    umi = "".join(rng.choice("ACGT") for _ in range(10))
                    depth = rng.choice([1, 1, 2, 5, 40])
                    fh.write(f"@PBMC.{cell}.{umi} cD:i:{depth}\n{seq}\n+\n{'I' * len(seq)}\n")
    return out, cells, receptors, expected


def test_a_cell_gets_its_receptor_back_full_length(tmp_path, droplet_library):
    from arda.singlecell import run, trim

    reads, cells, receptors, _expected = droplet_library
    summary = run(reads, tmp_path / "sc", cell_from="migec", threads=4)
    assert summary["cells"] == len(cells)

    contigs: dict[str, list[str]] = {}
    name, buf = None, []
    for line in (tmp_path / "sc.contigs.fasta").read_text().splitlines():
        if line.startswith(">"):
            if name:
                contigs.setdefault(name, []).append("".join(buf))
            name = line[1:].split("_contig_")[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        contigs.setdefault(name, []).append("".join(buf))

    # Never: the contig is a substring of the RAW template but its length is scored against the
    # TRIMMED one. Two behaviours, both correct, pull the two ends apart: molecules are windows
    # at random starts, so the first and last few bases are usually in no window (measured, 814
    # of 822); and one fixture record ends in a 120 nt polyA tail that `trim` removes on purpose,
    # while a window ending in eleven A's keeps them, so the contig is not a substring of the
    # trimmed template either. Scored against the raw template it reads 1,493 of 1,620 and looks
    # like an assembly failure.
    singletons = {c: m[0] for c, m in cells.items() if len(m) == 1}
    for cell, which in singletons.items():
        template = trim(receptors[which])
        best = max(contigs.get(cell, [""]), key=len)
        assert best in receptors[which], cell
        assert len(best) >= 0.95 * len(template), f"{cell}: {len(best)} of {len(template)}"


def test_the_junction_is_called_and_the_doublets_are_flagged(tmp_path, droplet_library):
    from arda.singlecell import run

    reads, cells, _receptors, expected = droplet_library
    run(reads, tmp_path / "sc", cell_from="migec", threads=4)

    with open(tmp_path / "sc.cells.tsv", newline="") as fh:
        rows = {r["cell_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == set(cells)
    # A cell holding one receptor gets back exactly the junction arda calls on that template,
    # and is not a doublet.
    scored = 0
    for cell, members in cells.items():
        if len(members) != 1:
            continue
        want = expected[members[0]].get("junction_aa")
        assert rows[cell]["cell_status"] != "doublet_candidate", cell
        if not want:
            continue
        scored += 1
        got = {rows[cell]["heavy_junction_aa"], rows[cell]["light_junction_aa"]} - {""}
        assert want in got, f"{cell}: expected {want}, got {got}"
    assert scored >= 25, f"only {scored} receptors carried a junction to score"

    doublets = [c for c, m in cells.items() if len(m) > 1]
    flagged = [c for c in doublets if rows[c]["cell_status"] == "doublet_candidate"]
    assert len(flagged) == len(doublets), f"{len(flagged)} of {len(doublets)} doublets flagged"


def test_the_report_and_every_table_are_written(tmp_path, droplet_library):
    from arda.scplot import PANELS
    from arda.singlecell import run

    reads, _cells, _receptors, _expected = droplet_library
    run(reads, tmp_path / "sc", cell_from="migec", threads=4, plot="svg", gnuplot="")

    report = json.loads((tmp_path / "sc.arda.json").read_text())
    assert report["contigs"] > 0 and report["contig_n50"] > 0
    for suffix in (".contigs.fasta", ".contigs.airr.tsv", ".chains.tsv", ".cells.tsv",
                   ".cell_rank.tsv", ".contig_lengths.tsv"):
        assert (tmp_path / f"sc{suffix}").exists(), suffix
    # No reference was given, so the two scored tables must NOT appear.
    assert not (tmp_path / "sc.sweep.tsv").exists()
    assert not (tmp_path / "sc.partition.tsv").exists()
    # `gnuplot=""` writes the scripts and renders nothing, which is the no-gnuplot path.
    for suffix, (stem, _title, _body) in PANELS.items():
        if (tmp_path / f"sc{suffix}").exists():
            assert (tmp_path / f"sc.{stem}.gp").exists(), stem
            assert not (tmp_path / f"sc.{stem}.svg").exists(), stem


def test_a_reference_adds_the_sweep_and_the_partition_scores(tmp_path, droplet_library):
    from arda.singlecell import run

    reads, cells, _receptors, _expected = droplet_library
    first = run(reads, tmp_path / "a", cell_from="migec", threads=4)

    # Score against our own primary calls: recall is then 1.0 by construction, which is what
    # makes this a test of the plumbing and not of the biology.
    with open(tmp_path / "a.chains.tsv", newline="") as fh:
        primary = [r for r in csv.DictReader(fh, delimiter="\t") if r["chain_rank"] == "1"]
    reference = tmp_path / "ref.tsv"
    with open(reference, "w") as fh:
        fh.write("cell_id\tlocus\tjunction_aa\n")
        for row in primary:
            fh.write(f"{row['cell_id']}\t{row['locus']}\t{row['junction_aa']}\n")

    second = run(reads, tmp_path / "b", cell_from="migec", threads=4, reference=reference)
    assert second["contigs"] == first["contigs"]
    assert second["reference_recall"] == pytest.approx(1.0)
    assert (tmp_path / "b.sweep.tsv").exists()
    with open(tmp_path / "b.partition.tsv", newline="") as fh:
        scores = {r["metric"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}
    # Homogeneity is 1 because no cluster mixes two reference clonotypes. Parsimony is NOT
    # asserted at 1: the doublet cells carry a superset label (two chains) of a class a singleton
    # cell also holds, so that class is genuinely split and the score says so. Asserting 1 here
    # would be asserting the fixture has no doublets.
    assert float(scores["homogeneity"]) == pytest.approx(1.0)
    assert 0.0 <= float(scores["parsimony"]) <= 1.0
    assert int(scores["classes_singleton"]) > 0


def test_restricting_to_called_cells_drops_the_rest(tmp_path, droplet_library):
    from arda.singlecell import run

    reads, cells, _receptors, _expected = droplet_library
    keep = sorted(cells)[:3]
    listing = tmp_path / "called.txt"
    listing.write_text("\n".join(keep) + "\n")
    summary = run(reads, tmp_path / "sc", cell_from="migec", threads=4, cells_file=listing)
    assert summary["cells"] == len(keep)


def test_a_panel_gnuplot_refuses_does_not_kill_the_run(tmp_path):
    """A library with no doublets has an empty doublet scatter, which gnuplot calls an invalid
    range. That must cost the one panel, not the run, and must not leave a 0-byte SVG behind."""
    import shutil

    from arda.scplot import draw

    if not shutil.which("gnuplot"):
        pytest.skip("gnuplot is not on PATH")
    prefix = tmp_path / "sc"
    # Every cell has one heavy chain, so every point in the scatter filters to 1/0.
    (tmp_path / "sc.cells.tsv").write_text(
        "cell\tcell_status\theavy_molecules_1\theavy_molecules_2\n"
        + "".join(f"C{i}\tpaired\t{10 + i}\t0\n" for i in range(20)))
    (tmp_path / "sc.contig_lengths.tsv").write_text(
        "length_bin\tcontigs\n300\t4\n325\t9\n350\t2\n")

    figures = draw(prefix, fmt="svg")
    assert (tmp_path / "sc.doublet_scatter.gp").exists()
    assert not (tmp_path / "sc.doublet_scatter.svg").exists()
    assert [f.name for f in figures] == ["sc.contig_lengths.svg"]
