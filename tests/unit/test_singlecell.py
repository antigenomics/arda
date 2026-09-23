"""The per-cell assembler and the chain-pairing rule, unit level.

The end-to-end acceptance test is `tests/synthetic/test_cells_end_to_end.py`, which tiles real
receptor sequences into molecule-shaped fragments and asks for the junction back.
"""

from __future__ import annotations

import random

import pytest

from arda.singlecell import (
    ADAPTER,
    Contig,
    assemble_cell,
    cell_rank,
    cell_summary,
    pair_chains,
    read_called_cells,
    read_molecules,
    read_records,
    threshold_sweep,
    trim,
)


def _random_dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# --- trimming -------------------------------------------------------------------------------


def test_the_adapter_and_everything_after_it_goes():
    assert trim("ACGTACGTAC" + ADAPTER + "GGGGTTTT") == "ACGTACGTAC"


def test_a_partial_adapter_stub_also_goes():
    # A consensus that stopped eight bases into the adapter still shares those eight with every
    # other such consensus, which is exactly enough to glue two transcripts together.
    assert trim("ACGTACGTAC" + ADAPTER[:8]) == "ACGTACGTAC"


def test_a_homopolymer_tail_goes():
    assert trim("ACGTACGTAC" + "A" * 20) == "ACGTACGTAC"


def test_a_sequence_with_neither_is_untouched():
    seq = _random_dna(200, 1)
    assert trim(seq) == seq


# --- assembly -------------------------------------------------------------------------------


def _tiles(template: str, width: int = 90, step: int = 30) -> list[str]:
    """Molecule-shaped windows that between them cover the whole template.

    The last window is anchored at the END, not at the last multiple of `step` -- a
    `range(0, len - width + 1, step)` alone stops short and the assembly is then correctly
    missing the tail, which reads as an assembler bug and is not one.
    """
    out = [template[i:i + width] for i in range(0, len(template) - width + 1, step)]
    if out[-1] != template[-width:]:
        out.append(template[-width:])
    return out


def test_tiles_of_one_template_assemble_back_into_it():
    template = _random_dna(500, 2)
    contigs = assemble_cell("CELL", _tiles(template))
    assert len(contigs) == 1
    assert contigs[0].sequence == template
    assert contigs[0].cell == "CELL"


def test_two_unrelated_templates_stay_two_contigs():
    a, b = _random_dna(400, 3), _random_dna(400, 4)
    contigs = assemble_cell("CELL", _tiles(a) + _tiles(b))
    assert len(contigs) == 2
    assert {c.sequence for c in contigs} == {a, b}


def test_an_untrimmed_adapter_wrecks_the_layout_and_trimming_restores_it():
    # An adapter suffix is a shared k-mer AND a 100%-identical overlap, so verification alone
    # cannot refuse it; it also sits at a mismatching offset against the REAL overlap of two
    # neighbouring molecules, which is what breaks the layout here. Either way the cure is the
    # trim, which is why `assemble_cells` trims before it assembles.
    a, b = _random_dna(400, 5), _random_dna(400, 6)
    dirty = [t + ADAPTER + "GGGGGGGGGGGGGGG" for t in _tiles(a) + _tiles(b)]
    assert len(assemble_cell("CELL", dirty)) < 2
    clean = assemble_cell("CELL", [trim(t) for t in dirty])
    assert {c.sequence for c in clean} == {a, b}


def test_the_deeper_molecule_wins_the_column():
    template = _random_dna(300, 7)
    tiles = _tiles(template)
    broken = tiles[0][:45] + ("A" if tiles[0][45] != "A" else "C") + tiles[0][46:]
    # One 1-read molecule disagrees with an 80-read one; depth weighting keeps the template.
    contigs = assemble_cell("CELL", [broken] + tiles, [1] + [80] * len(tiles))
    assert contigs[0].sequence == template
    assert contigs[0].reads == 1 + 80 * len(tiles)
    assert contigs[0].molecules == len(tiles) + 1


def test_a_lone_molecule_is_not_a_contig_by_default():
    assert assemble_cell("CELL", [_random_dna(200, 8)]) == []
    assert len(assemble_cell("CELL", [_random_dna(200, 8)], min_molecules=1)) == 1


def test_sequences_shorter_than_the_seed_are_skipped_not_crashed():
    template = _random_dna(300, 9)
    contigs = assemble_cell("CELL", ["ACGT"] + _tiles(template))
    assert contigs[0].sequence == template


def test_the_contig_id_is_the_cellranger_dialect():
    contig = Contig(cell="AAACCTGAGAAACCAT", index=2, sequence="ACGT", molecules=3, reads=9)
    assert contig.contig_id == "AAACCTGAGAAACCAT_contig_2"


# --- reading --------------------------------------------------------------------------------


def test_the_read_depth_comes_off_the_fastq_comment(tmp_path):
    path = tmp_path / "c.fq"
    path.write_text(
        "@PBMC.AAACCTGAGAAACCAT.CGTTTTTATC RX:Z:CGTTTTTATC\tcD:i:37\nACGT\n+\nIIII\n"
        "@PBMC.AAACCTGAGAAACCAT.TTTTTTTTTT\nACGA\n+\nIIII\n"
    )
    assert [d for _, _, d in read_records(path)] == [37, 1]


def test_a_quality_line_starting_with_an_at_sign_is_not_a_header(tmp_path):
    # Phred 31 is '@'. A reader that looks for '@' at a line start rather than counting four
    # lines per record loses the file from here on.
    path = tmp_path / "c.fq"
    path.write_text("@PBMC.AAACCTGAGAAACCAT.CGTTTTTATC\nACGT\n+\n@@@@\n"
                    "@PBMC.AAACCTGAGAAACCAT.TTTTTTTTTT\nACGA\n+\nIIII\n")
    assert len(list(read_records(path))) == 2


def test_grouping_by_cell_uses_the_named_dialect(tmp_path):
    path = tmp_path / "c.fa"
    path.write_text(">PBMC.AAACCTGAGAAACCAT.CGTTTTTATC\nACGT\n"
                    ">PBMC.AAACCTGAGAAACCAT.TTTTTTTTTT\nACGA\n"
                    ">PBMC.TTTCCTGAGAAACCAT.TTTTTTTTTT\nACGC\n")
    grouped = read_molecules(path, cell_from="migec")
    assert sorted(grouped) == ["AAACCTGAGAAACCAT", "TTTCCTGAGAAACCAT"]
    assert len(grouped["AAACCTGAGAAACCAT"]) == 2


def test_a_called_cell_list_reads_bare_or_from_a_migec_table(tmp_path):
    bare = tmp_path / "bare.txt"
    bare.write_text("AAACCTGAGAAACCAT\nTTTCCTGAGAAACCAT-1\n")
    assert read_called_cells(bare) == {"AAACCTGAGAAACCAT", "TTTCCTGAGAAACCAT"}

    table = tmp_path / "cells.tsv"
    table.write_text("cell\tmolecules\tcalled\n"
                     "AAACCTGAGAAACCAT\t120\ttrue\n"
                     "TTTCCTGAGAAACCAT\t2\tfalse\n")
    assert read_called_cells(table) == {"AAACCTGAGAAACCAT"}


# --- pairing --------------------------------------------------------------------------------


def _chain(cell, locus, junction, molecules, productive="T"):
    return {"cell_id": cell, "locus": locus, "junction_aa": junction,
            "molecules": molecules, "reads": molecules * 3, "productive": productive,
            "sequence": "A" * 500}


def test_one_heavy_and_one_light_chain_is_a_paired_cell():
    rows = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRA", "CAVRF", 20)]
    chains = pair_chains(rows)
    assert {r["status"] for r in chains} == {"primary"}
    assert cell_summary(chains)[0]["cell_status"] == "paired"


def test_a_well_supported_second_heavy_chain_is_a_doublet_candidate():
    rows = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRB", "CASSBF", 30),
            _chain("C1", "TRA", "CAVRF", 20)]
    chains = pair_chains(rows)
    second = next(r for r in chains if r["junction_aa"] == "CASSBF")
    assert second["status"] == "doublet_candidate"
    assert second["fraction_of_top"] == pytest.approx(0.75)
    assert cell_summary(chains)[0]["cell_status"] == "doublet_candidate"


def test_a_second_light_chain_is_allelic_inclusion_and_not_a_doublet():
    rows = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRA", "CAVRF", 20),
            _chain("C1", "TRA", "CAVSF", 18)]
    chains = pair_chains(rows)
    assert next(r for r in chains if r["junction_aa"] == "CAVSF")["status"] == "secondary"
    assert cell_summary(chains)[0]["cell_status"] == "paired"


def test_a_non_productive_extra_chain_is_never_believed():
    rows = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRB", "CASSBF", 30, "F")]
    chains = pair_chains(rows)
    assert next(r for r in chains if r["junction_aa"] == "CASSBF")["status"] == "extra"
    assert pair_chains(rows, require_productive=False)[1]["status"] == "doublet_candidate"


def test_a_thin_or_low_fraction_extra_chain_is_extra():
    thin = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRB", "CASSBF", 2)]
    assert pair_chains(thin)[1]["status"] == "extra"
    low = [_chain("C1", "TRB", "CASSAF", 400), _chain("C1", "TRB", "CASSBF", 10)]
    assert pair_chains(low)[1]["status"] == "extra"
    assert pair_chains(low, min_extra_fraction=0.0)[1]["status"] == "doublet_candidate"


def test_two_contigs_of_one_junction_are_one_chain_with_the_support_summed():
    rows = [_chain("C1", "TRB", "CASSAF", 12), _chain("C1", "TRB", "CASSAF", 8)]
    chains = pair_chains(rows)
    assert len(chains) == 1
    assert chains[0]["molecules"] == 20
    assert chains[0]["contigs"] == 2


def test_productive_anywhere_makes_the_chain_productive():
    rows = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRB", "CASSBF", 20, "F"),
            _chain("C1", "TRB", "CASSBF", 15, "T")]
    chains = pair_chains(rows)
    assert next(r for r in chains if r["junction_aa"] == "CASSBF")["productive"] == "T"


def test_a_row_missing_a_cell_a_locus_or_a_junction_is_dropped():
    rows = [_chain("", "TRB", "CASSAF", 40), _chain("C1", "", "CASSAF", 40),
            _chain("C1", "TRB", "", 40)]
    assert pair_chains(rows) == []


def test_the_knee_curve_is_molecules_and_is_ranked_descending():
    rows = cell_rank({"A": 3, "B": 100, "C": 20})
    assert [r["cell_id"] for r in rows] == ["B", "C", "A"]
    assert [r["rank"] for r in rows] == [1, 2, 3]
    assert [r["molecules"] for r in rows] == [100, 20, 3]


def test_the_sweep_leaves_the_primary_chain_out_and_reports_both_gates():
    rows = [_chain("C1", "TRB", "CASSAF", 40), _chain("C1", "TRB", "CASSBF", 6),
            _chain("C1", "TRB", "CASSCF", 4, "F")]
    chains = pair_chains(rows)
    sweep = threshold_sweep(chains, {("C1", "TRB", "CASSBF")}, thresholds=(1, 5))
    assert {r["require_productive"] for r in sweep} == {"F", "T"}
    loose = next(r for r in sweep if r["require_productive"] == "F" and r["min_molecules"] == 1)
    strict = next(r for r in sweep if r["require_productive"] == "T" and r["min_molecules"] == 1)
    assert loose["kept"] == 2 and loose["precision"] == pytest.approx(0.5)
    assert strict["kept"] == 1 and strict["precision"] == pytest.approx(1.0)
    # Recall shares the loose denominator, so turning the gate on can only lower it.
    assert strict["recall"] == pytest.approx(1.0)


def test_molecules_per_cell_is_the_input_count_not_the_contig_sum(tmp_path):
    # A molecule that joined no contig is missing from a contig-membership sum, and one shared
    # between two phased haplotypes is counted twice. Either way the knee plot's y axis stops
    # being unique molecules per cell, which is the one thing it claims to be.
    from arda.singlecell import assemble_cells

    template = _random_dna(400, 21)
    path = tmp_path / "c.fq"
    with open(path, "w") as fh:
        # The migec dialect needs the trailing field to be [ACGTN]+, so a UMI with digits in it
        # parses to no cell at all and the record is silently dropped.
        for n, seq in enumerate(_tiles(template)):
            umi = "".join("ACGT"[(n >> (2 * b)) & 3] for b in range(5)) + "ACGTT"
            fh.write(f"@PBMC.AAACCTGAGAAACCAT.{umi}\n{seq}\n+\n{'I' * len(seq)}\n")
        # One molecule that overlaps nothing: it is a real molecule of this cell and joins no contig.
        lone = _random_dna(200, 22)
        fh.write(f"@PBMC.AAACCTGAGAAACCAT.TTTTTTTTTT\n{lone}\n+\n{'I' * len(lone)}\n")

    contigs, report = assemble_cells(path, tmp_path / "out.fasta", cell_from="migec")
    counted = report.molecules_by_cell["AAACCTGAGAAACCAT"]
    assert counted == len(_tiles(template)) + 1
    assert counted > sum(c.molecules for c in contigs)


# --- the knee -------------------------------------------------------------------------------


def test_the_knee_is_found_on_a_curve_that_has_one():
    from arda.singlecell import find_knee

    # 100 cells at ~500 molecules, then 5,000 ambient barcodes at 1-3. The break is at 100.
    sizes = [500 - i for i in range(100)] + [3, 2, 1] * 1667
    rank, molecules = find_knee(sizes)
    assert 80 <= rank <= 130, rank
    assert molecules > 300


def test_an_ambient_only_curve_reports_no_knee():
    from arda.singlecell import find_knee

    # A global maximum always exists, so this returns a rank unless it is guarded. Its "knee"
    # would be a threshold that calls every barcode a cell.
    assert find_knee([3] * 1000 + [2] * 2000 + [1] * 5000) == (0, 0)


def test_a_flat_or_tiny_curve_reports_no_knee():
    from arda.singlecell import find_knee

    assert find_knee([]) == (0, 0)
    assert find_knee([7, 7]) == (0, 0)
    assert find_knee([5] * 500) == (0, 0)


def test_the_guard_can_be_relaxed_and_then_the_degenerate_answer_comes_back():
    from arda.singlecell import find_knee

    # Naming the floor is the point: at 0 the ambient-only curve reports its meaningless corner,
    # which is exactly what the default refuses.
    assert find_knee([3] * 1000 + [2] * 2000 + [1] * 5000, min_ratio_to_mean=0.0) != (0, 0)


def test_the_knee_row_is_marked_in_the_rank_table(tmp_path):
    from arda.singlecell import cell_rank, find_knee

    sizes = {f"C{i:05d}": v for i, v in
             enumerate([500 - i for i in range(100)] + [3, 2, 1] * 1667)}
    ranked = cell_rank(sizes)
    rank, _molecules = find_knee([r["molecules"] for r in ranked])
    assert rank > 0
    assert ranked[rank - 1]["molecules"] == find_knee([r["molecules"] for r in ranked])[1]


def test_two_samples_in_one_file_are_refused(tmp_path):
    # 10x barcodes come from a fixed whitelist, so two samples share them by design: grouping on
    # the barcode alone would silently merge two cells and the result reads as a doublet.
    path = tmp_path / "both.fa"
    path.write_text(">PBMC1.AAACCTGAGAAACCAT.CGTTTTTATC\nACGT\n"
                    ">PBMC2.AAACCTGAGAAACCAT.TTTTTTTTTT\nACGA\n")
    with pytest.raises(ValueError, match="2 samples"):
        read_molecules(path, cell_from="migec")


def test_one_sample_in_one_file_is_not_refused(tmp_path):
    path = tmp_path / "one.fa"
    path.write_text(">PBMC.AAACCTGAGAAACCAT.CGTTTTTATC\nACGT\n"
                    ">PBMC.TTTCCTGAGAAACCAT.TTTTTTTTTT\nACGA\n")
    assert len(read_molecules(path, cell_from="migec")) == 2


def test_an_unfiltered_droplet_library_is_refused_by_name(tmp_path):
    from arda.singlecell import assemble_cells

    # 40 cells at 30-70 molecules, then 4,000 ambient barcodes at 2. Without --cells every one of
    # those 4,000 is assembled and lands in the chain table looking like a cell.
    lines = []
    for i in range(40):
        for m in range(30 + i):
            lines.append(f">S.{_bc(i)}.{_bc(m)[:10]}\nACGTACGTACGT\n")
    for i in range(4000):
        for m in range(2):
            lines.append(f">S.{_bc(100000 + i)}.{_bc(m)[:10]}\nACGTACGTACGT\n")
    path = tmp_path / "all.fa"
    path.write_text("".join(lines))
    with pytest.raises(ValueError, match="--cells"):
        assemble_cells(path, tmp_path / "out.fa", cell_from="migec")
    # The same input with the called set passes.
    contigs, report = assemble_cells(path, tmp_path / "out.fa", cell_from="migec",
                                     cells={_bc(i) for i in range(40)})
    assert report.cells == 40


def _bc(n: int) -> str:
    """A distinct 16 nt barcode per integer."""
    out = []
    for _ in range(16):
        out.append("ACGT"[n & 3])
        n >>= 2
    return "".join(out)


def test_a_bulk_consensus_is_refused_rather_than_assembling_nothing(tmp_path):
    # `<sample>.<umi>` parses under the migec dialect with cell=None, so without this the whole
    # command succeeds having found zero cells.
    path = tmp_path / "bulk.fa"
    path.write_text(">PBMC.CGTTTTTATC\nACGT\n>PBMC.TTTTTTTTTT\nACGA\n")
    with pytest.raises(ValueError, match="no cell"):
        read_molecules(path, cell_from="migec")
