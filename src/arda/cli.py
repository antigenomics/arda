"""arda command-line interface.

Subcommands:

* ``arda annotate``  — map input sequences and emit AIRR TSV (Phase 2).
* ``arda build-db``  — (re)build the curated reference DB from IMGT + IgBLAST (Phase 1).
* ``arda info``      — show resolved tool/data paths and versions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

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


@app.command("igblast")
def igblast_cmd(
    input: Path = typer.Option(..., "--input", "-i", help="Query FASTA/FASTQ."),
    output: Path = typer.Option(..., "--output", "-o", help="Merged AIRR TSV (gold standard)."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="igblast threads (0 = all cores)."),
) -> None:
    """Gold-standard AIRR alignment with IgBLAST across all annotatable loci."""
    import os
    from .refbuild.gold import igblast_reads

    igblast_reads(input, output, organism=organism,
                  num_threads=threads or (os.cpu_count() or 1))
    typer.echo(f"[arda] igblast AIRR -> {output}")


rnaseq_app = typer.Typer(add_completion=False, help="RNA-seq filter/map/correct pipeline.")
app.add_typer(rnaseq_app, name="rnaseq")


@rnaseq_app.command("map")
def rnaseq_map(
    output: Path = typer.Option(..., "--output", "-o", help="AIRR TSV of mapped reads only."),
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="mmseqs threads (0 = all cores)."),
    sensitivity: Optional[float] = typer.Option(None, help="mmseqs -s (default: tuned)."),
    strand: str = typer.Option("both", help="'both' strands or 'forward'."),
    chunk_size: int = typer.Option(
        200000, help="Reads per streaming chunk (larger = faster, more RAM)."),
    map_d: bool = typer.Option(
        True, "--map-d/--no-map-d",
        help="Map D segments (VDJ loci). --no-map-d is ~19% faster and yields an identical read "
             "set, v_call and junction — only d_call is dropped. Use it for filtering or for the "
             "map|correct clonotype pipeline, which never keys on D."),
    reconstruct: bool = typer.Option(
        False, "--reconstruct", help="Merge overlapping paired mates into one fragment."),
    min_score: float = typer.Option(
        75.0, "--min-score",
        help="Min MMseqs2 bit score to keep a mapped read (0 = keep all). Once the reference carries "
             "J+C scaffolds the recall/score curve is flat over 0-75, so the exact value barely "
             "matters; precision still rises with it (0.933 -> 0.959 across 40 -> 75)."),
    max_seqs: int = typer.Option(
        300, "--max-seqs",
        help="MMseqs2 target hits per read. Does NOT change which reads are kept, only which "
             "V/J scaffold wins: 50 -> 300 lifts V-gene concordance 83%% -> 96%%. Use --kmer for "
             "memory, not this."),
    kmer: int = typer.Option(
        12, "--kmer", "-k",
        help="MMseqs2 -k. THE memory knob: the nucleotide prefilter allocates 4**k index entries, so "
             "peak RSS tracks 4**k*8 bytes almost exactly (k=15 -> 8.4 GB, 13 -> 697 MB, 12 -> 298 MB, "
             "11 -> 202 MB). Recall and precision are INVARIANT over k=11..14; 12 is also the fastest "
             "measured. 0 = MMseqs2 default."),
    drop_constant_only: bool = typer.Option(
        True, "--drop-constant-only/--keep-constant-only",
        help="Drop reads whose alignment lies wholly inside the constant region (tstart >= t_vjend). "
             "They are real receptor mRNA but carry no V(D)J, hence no clonotype. Dropping them takes "
             "precision 0.756 -> 0.968 at unchanged recall. Keep them if you are assembling contigs."),
    emit_reads: Optional[Path] = typer.Option(
        None, "--emit-reads", help="Also write the mapped reads as FASTA (for handoff)."),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Filter + map receptor reads from RNA-seq; write only mapped reads (AIRR TSV)."""
    from .rnaseq.map import map_rnaseq

    rep = map_rnaseq(r1, output, r2=r2, organism=organism, threads=threads,
                     sensitivity=sensitivity, strand=strand, chunk_size=chunk_size,
                     map_d=map_d, reconstruct=reconstruct, min_score=min_score,
                     max_seqs=max_seqs, kmer=(None if kmer == 0 else kmer),
                     drop_constant_only=drop_constant_only,
                     emit_reads=emit_reads, report_path=report)
    typer.echo(
        f"[arda] {rep.mapped_reads}/{rep.total_reads} reads mapped "
        f"({rep.mapped_fraction * 100:.2f}%) in {rep.wall_seconds:.1f}s "
        f"({rep.reads_per_second:.0f} reads/s); loci={rep.per_locus}")


@rnaseq_app.command("correct")
def rnaseq_correct(
    input: Path = typer.Option(..., "--input", "-i", help="Mapped-reads AIRR TSV (from `map`)."),
    output: Path = typer.Option(..., "--output", "-o", help="Corrected clonotype table TSV."),
    max_mismatches: int = typer.Option(2, help="Max CDR3 substitutions to collapse."),
    ratio: float = typer.Option(0.05, help="Parent:child count ratio (0.05 = 20x per mismatch)."),
    require_vj: bool = typer.Option(
        False, "--require-vj", help="Only collapse neighbours sharing V and J calls."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions",
        help="Keep only complete junctions (span C104..[FW]118, in frame, no stop). Reads that "
             "stop short of the anchor yield a prefix of a junction, not a clonotype."),
    read_map: Optional[Path] = typer.Option(
        None, "--read-map", help="Write sequence_id -> corrected junction map."),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Collapse CDR3 sequencing errors (parent:child ratio) into clonotypes."""
    from .rnaseq.correct import correct_airr

    rep = correct_airr(input, output, max_mismatches=max_mismatches, ratio=ratio,
                       require_vj=require_vj, complete_only=complete_only,
                       read_map=read_map, report_path=report)
    typer.echo(
        f"[arda] {rep.clonotypes_in} -> {rep.clonotypes_out} clonotypes "
        f"({rep.collapsed} collapsed) over {rep.reads} reads"
        + (f"; dropped {rep.reads_incomplete}/{rep.reads_with_junction} incomplete junctions"
           if rep.reads_incomplete else ""))


if __name__ == "__main__":
    app()
