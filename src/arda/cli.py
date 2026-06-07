"""arda command-line interface.

Subcommands:

* ``arda annotate``  — map input sequences and emit AIRR TSV (Phase 2).
* ``arda build-db``  — (re)build the curated reference DB from IMGT + IgBLAST (Phase 1).
* ``arda info``      — show resolved tool/data paths and versions.
"""

from __future__ import annotations

from pathlib import Path

import typer

from . import __version__

app = typer.Typer(add_completion=False, help="Antigen Receptor Domain Annotation")


@app.command()
def info() -> None:
    """Show resolved paths and external tool availability."""
    from .paths import project_root, bin_dir, data_dir, database_dir

    typer.echo(f"arda {__version__}")
    typer.echo(f"project_root : {project_root()}")
    typer.echo(f"bin_dir      : {bin_dir()}")
    typer.echo(f"data_dir     : {data_dir()}")
    typer.echo(f"database_dir : {database_dir()}")

    try:
        from .mmseqs import mmseqs_binary

        typer.echo(f"mmseqs       : {mmseqs_binary()}")
    except Exception as exc:  # noqa: BLE001
        typer.echo(f"mmseqs       : NOT FOUND ({exc})")

    try:
        from ._markup import __version__ as markup_version

        typer.echo(f"_markup ext  : {markup_version}")
    except Exception as exc:  # noqa: BLE001
        typer.echo(f"_markup ext  : NOT BUILT ({exc})")


@app.command("build-db")
def build_db(
    organism: str = typer.Option(
        "all", help="Organism to build (or 'all' for every supported organism)."
    ),
) -> None:
    """Build the curated reference database (Phase 1)."""
    from .refbuild.build import build

    build(organism)


@app.command()
def annotate(
    input: Path = typer.Option(..., "--input", "-i", help="Input FASTA/FASTQ."),
    output: Path = typer.Option(..., "--output", "-o", help="Output AIRR TSV."),
    organism: str = typer.Option("human", help="Reference organism."),
    seqtype: str = typer.Option("nt", help="Input sequence type: nt or aa."),
    threads: int = typer.Option(0, help="mmseqs threads (0 = all cores)."),
    strand: str = typer.Option("both", help="nt only: 'both' strands or 'forward'."),
    chunk_size: int = typer.Option(
        50000, help="Reads per streaming chunk (bounds memory for large FASTQ)."),
    map_d: bool = typer.Option(
        True, "--map-d/--no-map-d",
        help="Map D segments (d_call/d2_call/np*) for VDJ loci; nt input only."),
) -> None:
    """Annotate FR/CDR regions and write an AIRR TSV (streamed, memory-bounded)."""
    from .annotate.mapper import annotate_file

    annotate_file(input, output, organism=organism, seqtype=seqtype,
                  threads=threads, strand=strand, chunk_size=chunk_size, map_d=map_d)


if __name__ == "__main__":
    app()
