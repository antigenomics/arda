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


@app.command("build-index")
def build_index_cmd(
    organism: str = typer.Option("all", help="Organism (or 'all')."),
    force: bool = typer.Option(False, "--force", help="Rebuild even if up to date."),
) -> None:
    """(Re)build the precompiled mmseqs DBs shipped under database/.

    These let `arda annotate` run out of the box; they are regenerated here for the
    locally installed mmseqs version when it differs from the shipped one.
    """
    from .annotate.mapper import build_index

    build_index(organism, force=force)


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


@app.command()
def split(
    input: Path = typer.Argument(..., help="Input FASTA/FASTQ."),
    out_dir: Path = typer.Argument(..., help="Directory for shard FASTA files."),
    shards: int = typer.Option(..., "--shards", help="Number of shards."),
) -> None:
    """Round-robin split an input into N shard FASTA files (for cluster runs)."""
    from .cluster import split as _split

    paths = _split(input, out_dir, shards)
    typer.echo(f"wrote {len(paths)} shards to {out_dir}")


@app.command()
def merge(
    shard_dir: Path = typer.Argument(..., help="Directory of per-shard AIRR TSVs."),
    output: Path = typer.Argument(..., help="Combined AIRR TSV."),
) -> None:
    """Concatenate per-shard AIRR TSVs into one (single header)."""
    from .cluster import merge as _merge

    _merge(shard_dir, output)
    typer.echo(f"merged -> {output}")


@app.command()
def slurm(
    input: Path = typer.Option(..., "--input", "-i", help="Input FASTA/FASTQ."),
    output: Path = typer.Option(..., "--output", "-o", help="Combined AIRR TSV."),
    work_dir: Path = typer.Option(Path("arda_slurm"), help="Scratch dir for shards/outputs."),
    shards: int = typer.Option(..., "--shards", help="SLURM array size."),
    organism: str = typer.Option("human"),
    seqtype: str = typer.Option("nt"),
    threads: int = typer.Option(8, help="cpus-per-task per array task."),
    strand: str = typer.Option("both"),
    map_d: bool = typer.Option(True, "--map-d/--no-map-d"),
    partition: str = typer.Option(None, help="SLURM partition."),
    time: str = typer.Option("04:00:00"),
    mem: str = typer.Option("8G"),
    submit: bool = typer.Option(False, "--submit", help="Run the generated submit.sh now."),
) -> None:
    """Write (and optionally submit) a SLURM submit.sh: split → array-annotate → merge."""
    import os
    import subprocess
    from .cluster import render_submit_script

    script = render_submit_script(
        input, output, work_dir, shards=shards, organism=organism, seqtype=seqtype,
        threads=threads, strand=strand, map_d=map_d, partition=partition, time=time,
        mem=mem, arda_mmseqs=os.environ.get("ARDA_MMSEQS"))
    work_dir.mkdir(parents=True, exist_ok=True)
    submit_sh = work_dir / "submit.sh"
    submit_sh.write_text(script)
    submit_sh.chmod(0o755)
    typer.echo(f"wrote {submit_sh}")
    if submit:
        subprocess.run(["bash", str(submit_sh)], check=True)
    else:
        typer.echo(f"submit with: bash {submit_sh}")


if __name__ == "__main__":
    app()
