"""`--cell-from` through a real `arda map` run: the S1 acceptance test.

The parser itself is covered in `tests/unit/test_cell.py`. What this adds is that the hook is in a
branch `map` actually reaches -- the obvious place for it, `mapper.py:1430`, is inside the
UNMAPPED-record branch, which `mapped_only=True` (what `map` always passes) skips, so a hook there
is dead code on the only path that matters and every test of the parser would still pass.
"""

from __future__ import annotations

import csv
import random
from pathlib import Path

import pytest

from tests.conftest import requires_human_db, requires_mmseqs

pytestmark = [requires_mmseqs, requires_human_db]

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "genbank_receptors.fa"


@pytest.fixture(scope="module")
def migec_named(tmp_path_factory):
    """The committed receptor fixture, renamed to migec's molecule-name grammar.

    Every third record carries a `.c<k>` contig suffix or a `.c<k>.<m>` contig+split pair, because
    both are emitted conditionally by `assemble` and a parser that assumes a fixed field count gets
    the barcode wrong on exactly those.
    """
    rng = random.Random(0)
    out = tmp_path_factory.mktemp("cellid") / "reads.fa"
    cells, n = {}, 0
    with open(FIXTURE) as fh, open(out, "w") as w:
        for line in fh:
            if not line.startswith(">"):
                w.write(line)
                continue
            n += 1
            cell = "".join(rng.choice("ACGT") for _ in range(16))
            umi = "".join(rng.choice("ACGT") for _ in range(10))
            suffix = "" if n % 3 else (".c2" if n % 2 else ".c1.2")
            name = f"PBMC.{cell}.{umi}{suffix}"
            cells[name] = cell
            w.write(f">{name}\n")
    assert n == 29
    return out, cells


def test_cell_id_is_emitted_and_is_the_barcode(tmp_path, migec_named):
    from arda.rnaseq.map import map_rnaseq

    reads, expected = migec_named
    airr = tmp_path / "out.airr.tsv"
    rep = map_rnaseq(reads, airr, organism="human", threads=4, cell_from="migec")
    assert rep.mapped_reads == 29

    with open(airr, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 29
    for r in rows:
        assert r["cell_id"] == expected[r["sequence_id"]]
        assert len(r["cell_id"]) == 16


def test_cell_id_goes_last_so_it_cannot_move_another_column(tmp_path, migec_named):
    """Every extra column is appended, in a fixed order, so a positional consumer is unaffected."""
    from arda.rnaseq.map import map_rnaseq

    reads, _ = migec_named
    plain, celled = tmp_path / "plain.tsv", tmp_path / "celled.tsv"
    map_rnaseq(reads, plain, organism="human", threads=4)
    map_rnaseq(reads, celled, organism="human", threads=4, cell_from="migec")

    a = plain.read_text().splitlines()[0].split("\t")
    b = celled.read_text().splitlines()[0].split("\t")
    assert b == a + ["cell_id"]


def test_without_the_flag_the_output_is_unchanged(tmp_path, migec_named):
    """A bulk run must not pay for, or notice, the single-cell path."""
    from arda.rnaseq.map import map_rnaseq

    reads, _ = migec_named
    one, two = tmp_path / "a.tsv", tmp_path / "b.tsv"
    map_rnaseq(reads, one, organism="human", threads=4)
    map_rnaseq(reads, two, organism="human", threads=4, cell_from="")
    assert one.read_text() == two.read_text()
    assert "cell_id" not in one.read_text().splitlines()[0]


def test_auto_finds_the_dialect(tmp_path, migec_named):
    from arda.rnaseq.map import map_rnaseq

    reads, expected = migec_named
    airr = tmp_path / "auto.airr.tsv"
    map_rnaseq(reads, airr, organism="human", threads=4, cell_from="auto")
    with open(airr, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows and all(r["cell_id"] == expected[r["sequence_id"]] for r in rows)
