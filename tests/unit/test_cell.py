"""The single-cell identifier parser: every dialect, and the ways it must refuse."""

import pytest

from arda.cell import _MIN_CELL_LEN, DIALECTS, make_parser, parse, sniff

CELL = "AAACCTGCAGCCTGTT"
UMI = "CGTTTTTATC"


@pytest.mark.parametrize(
    "dialect,sequence_id,cell",
    [
        ("cellranger", "AAACCTGAGAAACCAT-1_contig_2", "AAACCTGAGAAACCAT"),
        ("cellranger", "AAACCTGAGAAACCAT_contig_1", "AAACCTGAGAAACCAT"),
        ("cellranger", "AAACCTGAGAAACCAT-1", None),          # a barcode is not a contig id
        ("migec", f"PBMC.{CELL}.{UMI}", CELL),
        ("migec", f"PBMC.{CELL}.{UMI}.c2", CELL),
        ("migec", f"PBMC.{CELL}.{UMI}.c2.3", CELL),
        ("migec", f"PBMC.{CELL}.{UMI}.3", CELL),
        ("migec", f"S1.{UMI}", None),                        # bulk: a UMI, no cell
        ("prefix", f"{CELL}_anything", CELL),
        ("prefix", "notdna_x", None),
    ],
)
def test_dialects(dialect, sequence_id, cell):
    key = parse(sequence_id, dialect)
    assert (key.cell if key else None) == cell


def test_migec_sample_id_may_contain_a_dot():
    """migec's `validate_sample_id` permits `.`, so a left-to-right split reads the wrong field."""
    key = parse(f"PBMC.rep2.{CELL}.{UMI}", "migec")
    assert key.cell == CELL
    assert key.umi == UMI
    assert key.sample == "PBMC.rep2"


def test_migec_suffixes_are_contig_then_molecule():
    key = parse(f"PBMC.{CELL}.{UMI}.c2.3", "migec")
    assert (key.contig, key.molecule) == (2, 3)
    # Emitted conditionally, so both defaults must survive a bare name.
    bare = parse(f"PBMC.{CELL}.{UMI}", "migec")
    assert (bare.contig, bare.molecule) == (1, 1)


def test_sniff_refuses_a_dna_like_bulk_sample_id():
    """The fabrication this guard exists for: a sample named over ACGT is not a cell.

    `TCGA.<umi>` parses consistently across every record, so neither the match rate nor a
    same-length check catches it -- only excluding `prefix` from the sniff set does.
    """
    assert sniff([f"TCGA.{UMI}", "TCGA.AAAAAAAAAA", "TCGA.CCCCCCCCCC"]) is None
    assert make_parser("prefix")(f"{CELL}_x") == CELL       # still reachable when NAMED


def test_sniff_refuses_a_barcode_shorter_than_the_floor():
    short = "A" * (_MIN_CELL_LEN - 1)
    assert sniff([f"S.{short}.{UMI}", f"S.{short}.AAAAAAAAAA"]) is None


def test_sniff_finds_the_real_dialects():
    assert sniff([f"PBMC.{CELL}.{UMI}", f"PBMC.{CELL}.AAAAAAAAAA"]) == "migec"
    assert sniff(["AAACCTGAGAAACCAT-1_contig_2", "TTTCCTGAGAAACCAT-1_contig_1"]) == "cellranger"
    assert sniff([]) is None


def test_auto_refuses_rather_than_guessing():
    with pytest.raises(ValueError, match="no dialect parses"):
        make_parser("auto", sample=[f"TCGA.{UMI}"])


def test_regex_needs_a_cell_group():
    with pytest.raises(ValueError, match="cell"):
        make_parser(regex=r"(.*)_x")
    assert make_parser(regex=r"^(?P<cell>[ACGT]+)_")(f"{CELL}_x") == CELL


def test_sampling_a_fastq_does_not_mistake_a_quality_line_for_a_header(tmp_path):
    """A FASTQ quality line can begin with `>` (Phred 29) or `@` (Phred 31).

    A per-line "starts with > or @" rule reads those as record headers, feeds `sniff` garbage, and
    surfaces as `--cell-from auto` refusing to detect a dialect it should have found.
    """
    from arda.rnaseq.map import _cell_sample

    fq = tmp_path / "r.fq"
    # Quality strings chosen so line 4 starts with '>' and line 8 starts with '@'.
    fq.write_text(
        f"@PBMC.{CELL}.{UMI}\nACGT\n+\n>III\n"
        f"@PBMC.{CELL}.AAAAAAAAAA\nACGT\n+\n@III\n"
    )
    ids = _cell_sample(fq, None)
    assert ids == [f"PBMC.{CELL}.{UMI}", f"PBMC.{CELL}.AAAAAAAAAA"]
    assert sniff(ids) == "migec"


def test_sampling_a_fasta_takes_every_record(tmp_path):
    from arda.rnaseq.map import _cell_sample

    fa = tmp_path / "r.fa"
    fa.write_text(f">PBMC.{CELL}.{UMI}\nACGT\nACGT\n>PBMC.{CELL}.AAAAAAAAAA\nACGT\n")
    assert _cell_sample(fa, None) == [f"PBMC.{CELL}.{UMI}", f"PBMC.{CELL}.AAAAAAAAAA"]


def test_unknown_dialect_is_named():
    with pytest.raises(ValueError, match="unknown cell-id dialect"):
        parse("x", "nonesuch")
    assert set(DIALECTS) == {"cellranger", "migec", "prefix"}
