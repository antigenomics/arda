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


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def _main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True,
        help="Show the arda version and exit."),
) -> None:
    """Antigen Receptor Domain Annotation."""


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
        help="Map D segments (d_call/d2_call/d_support/np*) for VDJ loci. Works on aa input "
             "too, against the D germlines' three translated frames -- informative for IGH, "
             "mostly silent for the TR loci, whose D is too short to survive in protein."),
) -> None:
    """Annotate FR/CDR regions and write an AIRR TSV (streamed, memory-bounded)."""
    from .annotate.mapper import annotate_file

    annotate_file(input, output, organism=organism, seqtype=seqtype,
                  threads=threads, strand=strand, chunk_size=chunk_size, map_d=map_d)


@app.command()
def markup(
    input: Path = typer.Option(..., "--input", "-i", help="Input TSV of CDR3/V/J records."),
    output: Path = typer.Option(..., "--output", "-o", help="Output TSV with markup + repair."),
    organism: str = typer.Option(
        "", help="Force one organism; default reads the species column."),
    cdr3_col: str = typer.Option("cdr3", help="Junction amino-acid column (C..[FW])."),
    v_col: str = typer.Option("v", help="V gene column."),
    j_col: str = typer.Option("j", help="J gene column."),
    species_col: str = typer.Option("species", help="Species column."),
    id_col: str = typer.Option("", help="Optional record-id column, used in the report."),
    vdjdb: bool = typer.Option(
        False, "--vdjdb", help="Read VDJdb column names (cdr3/v.segm/j.segm/species)."),
    max_replace: int = typer.Option(
        1, help="Repair edits at most this far from the conserved anchor; "
                "edits further in are reported but not applied."),
    d_posterior: bool = typer.Option(
        False, "--d-posterior",
        help="Also infer the D gene and its position from the junction length prior "
             "and the amino-acid match (human IGH/TRB/TRD, mouse TRB)."),
    report: Path = typer.Option(
        None, "--report", help="Write a human-readable fix log here ('-' for stdout)."),
    show_ok: bool = typer.Option(
        False, "--show-ok", help="List correct records in the report too, not just fixed/failed."),
) -> None:
    """Mark up and repair bare (junction_aa, V, J) records; emit vdjdb-style cdr3fix.

    The CDR3 column is the *junction*: Cys104 through Phe/Trp118, both included --
    the convention VDJdb's `cdr3` column uses.
    """
    import polars as pl

    from .cdr3fix import format_report, markup_records, to_frame

    if vdjdb:
        cdr3_col, v_col, j_col, species_col = "cdr3", "v.segm", "j.segm", "species"

    df = pl.read_csv(input, separator="\t", infer_schema_length=0)
    records = markup_records(df, cdr3=cdr3_col, v=v_col, j=j_col, species=species_col,
                             sequence_id=id_col or None, organism=organism or None,
                             max_replace=max_replace)
    out = to_frame(records)
    if d_posterior:
        from .dpost import posterior_d

        posts = [posterior_d(r.cdr3_repaired, r.v_call, r.j_call, r.species)
                 for r in records]
        out = out.with_columns([
            pl.Series("d_call", [p.d_call if p else "" for p in posts]),
            pl.Series("d_posterior", [round(p.posterior, 4) if p else None for p in posts]),
            pl.Series("d_entropy", [round(p.entropy, 3) if p else None for p in posts]),
            pl.Series("d_support_aa", [p.support_aa if p else None for p in posts]),
            pl.Series("d_start", [p.d_start if p else None for p in posts]),
            pl.Series("d_start_ci90",
                      [f"{p.d_start_ci90[0]}-{p.d_start_ci90[1]}" if p else "" for p in posts]),
        ])
    out.write_csv(output, separator="\t")

    if report is not None:
        text = format_report(records, show_ok=show_ok)
        if str(report) == "-":
            typer.echo(text)
        else:
            Path(report).write_text(text)

    n_fixed = sum(r.fix_needed for r in records)
    n_bad = sum(not r.good for r in records)
    typer.echo(f"{len(records)} records -> {output}  ({n_fixed} repaired, {n_bad} failed)")


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
    max_subs: int = typer.Option(2, help="Max substitutions between an error child and its parent."),
    max_indel: int = typer.Option(
        0, help="Max indel bases searched (default 0). 1-2 bp indel errors are frameshifts already "
                "dropped by --complete-only, so this only helps with --all-junctions; multi-bp SHM "
                "indels are kept regardless."),
    error_rate: float = typer.Option(
        0.001, help="Per-BASE substitution error rate (~Phred 30). Length-scaled: the per-sub "
                    "collapse prob is error_rate*junction_len, ~1/20 at a 45 nt (15 aa) junction."),
    indel_rate: float = typer.Option(
        0.001, help="Per-BASE indel error rate (length-scaled). Multi-bp (SHM) indels are kept."),
    require_vj: bool = typer.Option(
        True, "--require-vj/--no-require-vj",
        help="Only collapse neighbours sharing V and J (a true error keeps the germline V/J call)."),
    error_method: str = typer.Option(
        "simple", help="simple = spanning-read counts; binom|betabinom = per-position read-depth "
                       "pileup for very low coverage."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions",
        help="Keep only complete junctions (span C104..[FW]118, in frame, no stop). Reads that "
             "stop short of the anchor yield a prefix of a junction, not a clonotype."),
    read_map: Optional[Path] = typer.Option(
        None, "--read-map", help="Write sequence_id -> corrected junction map."),
    extra_airr: Optional[Path] = typer.Option(
        None, "--extra-airr",
        help="Stage-3 assembled-reads AIRR (from `assemble`) to fold into the clonotype table."),
    organism: str = typer.Option("human", help="Reference organism (used only to map D)."),
    map_d: bool = typer.Option(
        True, "--map-d/--no-map-d",
        help="Add d_call/d2_call/d_support to each clonotype, mapped into its error-corrected "
             "junction (once per clonotype, not per read)."),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Collapse CDR3 sequencing errors into clonotypes (per-substitution/indel error model)."""
    from .rnaseq.correct import correct_airr

    rep = correct_airr(input, output, organism=organism, map_d=map_d,
                       max_subs=max_subs, max_indel=max_indel, error_rate=error_rate,
                       indel_rate=indel_rate, require_vj=require_vj, error_method=error_method,
                       complete_only=complete_only, read_map=read_map, extra_airr=extra_airr,
                       report_path=report)
    typer.echo(
        f"[arda] {rep.clonotypes_in} -> {rep.clonotypes_out} clonotypes "
        f"({rep.collapsed} collapsed) over {rep.reads} reads"
        + (f"; dropped {rep.reads_incomplete}/{rep.reads_with_junction} incomplete junctions"
           if rep.reads_incomplete else ""))


@rnaseq_app.command("assemble")
def rnaseq_assemble(
    input: Path = typer.Option(..., "--input", "-i", help="Mapped-reads AIRR TSV (from `map`)."),
    output: Path = typer.Option(
        ..., "--output", "-o",
        help="Assembled-reads AIRR TSV: one row per rescued read carrying its contig's junction."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="mmseqs threads for re-annotation (0 = all cores)."),
    map_d: bool = typer.Option(
        True, "--map-d/--no-map-d",
        help="Call D (and tandem D-D) on each assembled contig and carry it onto its member "
             "reads. An ultralong CDR3 is where a D-D is most likely and least visible."),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Stage 3 — assemble long-CDR3 contigs the reads don't individually span.

    Reconstructs clonotypes (V(D)J ultralong CDR3s, ~20-40 aa) that ``map`` filters but that no
    single 100-150 bp read spans, by anchored greedy overlap-extension over the mapped reads.
    Feed the result to ``correct --extra-airr`` (``run`` does this automatically).
    """
    from .rnaseq.assemble import assemble_contigs

    rep = assemble_contigs(input, output, organism=organism, threads=threads, map_d=map_d,
                           report_path=report)
    typer.echo(
        f"[arda] assemble: {rep.contigs_complete}/{rep.contigs} complete contigs from "
        f"{rep.seeds} seeds; rescued {rep.reads_rescued} reads")


@rnaseq_app.command("run")
def rnaseq_run(
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
    out_prefix: str = typer.Option(
        ..., "--out-prefix", "-p",
        help="Output basename. Writes <prefix>.airr.tsv, <prefix>.clones.tsv, <prefix>.arda.json."),
    out_dir: Path = typer.Option(Path("."), "--out-dir", "-d", help="Directory for the outputs."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="mmseqs threads (0 = all cores)."),
    reconstruct: bool = typer.Option(
        False, "--reconstruct", help="Merge overlapping paired mates into one fragment."),
    min_score: float = typer.Option(
        75.0, "--min-score", help="Min MMseqs2 bit score to keep a mapped read (0 = keep all)."),
    kmer: int = typer.Option(
        12, "--kmer", "-k",
        help="MMseqs2 -k (the memory knob; k=12 ~= 298 MB peak RSS). 0 = MMseqs2 default (~8 GB)."),
    assemble: bool = typer.Option(
        True, "--assemble/--no-assemble",
        help="Stage 3: assemble long-CDR3 contigs the reads don't individually span (V(D)J "
             "ultralong, ~20-40 aa) and fold them into the clonotype table. Recovers ~95%% of the "
             "abundant long clones a filter-only pass misses; --no-assemble is faster (map+correct "
             "only)."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions",
        help="Keep only complete junctions when forming clonotypes."),
    map_d: bool = typer.Option(
        True, "--map-d/--no-map-d",
        help="Map D segments. Off skips D in all three stages; the clonotype table then carries "
             "no d_call/d2_call."),
) -> None:
    """One-shot RNA-seq -> clonotypes for pipeline integration: ``map`` -> ``assemble`` -> ``correct``.

    Runs the three RNA-seq stages with the shipped defaults, writing under ``--out-dir``:

    * ``<prefix>.airr.tsv``            -- Stage-1 mapped reads (AIRR Rearrangement)
    * ``<prefix>.assembled.airr.tsv``  -- Stage-3 assembled long-CDR3 reads (if ``--assemble``)
    * ``<prefix>.clones.tsv``          -- corrected clonotype table (folds in the assembled clones)
    * ``<prefix>.arda.json``           -- merged run report

    Use the ``map`` / ``assemble`` / ``correct`` commands separately to tune their knobs.
    """
    import json

    from .rnaseq.correct import correct_airr
    from .rnaseq.map import map_rnaseq

    out_dir.mkdir(parents=True, exist_ok=True)
    airr = out_dir / f"{out_prefix}.airr.tsv"
    clones = out_dir / f"{out_prefix}.clones.tsv"
    report = out_dir / f"{out_prefix}.arda.json"

    mrep = map_rnaseq(r1, airr, r2=r2, organism=organism, threads=threads,
                      reconstruct=reconstruct, min_score=min_score, map_d=map_d,
                      kmer=(None if kmer == 0 else kmer))
    typer.echo(
        f"[arda] map: {mrep.mapped_reads}/{mrep.total_reads} reads mapped "
        f"({mrep.mapped_fraction * 100:.2f}%); loci={mrep.per_locus}")

    arep = None
    extra = None
    if assemble:
        from .rnaseq.assemble import assemble_contigs
        extra = out_dir / f"{out_prefix}.assembled.airr.tsv"
        arep = assemble_contigs(airr, extra, organism=organism, threads=threads, map_d=map_d)
        typer.echo(
            f"[arda] assemble: {arep.contigs_complete}/{arep.contigs} complete contigs from "
            f"{arep.seeds} seeds; rescued {arep.reads_rescued} reads")

    crep = correct_airr(airr, clones, organism=organism, map_d=map_d,
                        complete_only=complete_only, extra_airr=extra)
    typer.echo(
        f"[arda] correct: {crep.clonotypes_in} -> {crep.clonotypes_out} clonotypes "
        f"({crep.collapsed} collapsed) over {crep.reads} reads")

    report.write_text(json.dumps(
        {"arda_version": __version__, "map": mrep.as_dict(),
         "assemble": arep.as_dict() if arep else None, "correct": crep.as_dict()},
        indent=2) + "\n")
    typer.echo(f"[arda] wrote {airr}, {clones}, {report}")


if __name__ == "__main__":
    app()
