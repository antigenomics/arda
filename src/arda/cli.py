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
    one_allele_per_gene: bool = typer.Option(
        False, "--one-allele-per-gene",
        help="Build scaffolds from one allele per gene (*01 where it exists, else the lowest). "
             "~4x smaller reference, no allele-level ambiguity.",
    ),
) -> None:
    """Build the curated reference database (Phase 1)."""
    from .refbuild.build import build

    build(organism, one_allele_per_gene=one_allele_per_gene)


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
    """Round-robin split an input into N shard FASTA files (amplicon / single-end).

    **Not for paired RNA-seq.** This writes FASTA, so quality is dropped, and it round-robins
    *records*, so a fragment's two mates land in different shards. Use ``arda rnaseq split``.
    """
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
    """Write (and optionally submit) a SLURM submit.sh: split → array-annotate → merge.

    **Not for paired RNA-seq.** This shards FASTA into ``arda annotate``: quality is dropped
    and mates are separated. It also has no Stage 2/3, so there is no clonotype table at the
    end. Use ``arda rnaseq slurm``.
    """
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
    limit: int = typer.Option(
        0, "--limit", "-n",
        help="Analyse only the first N reads (single-end) / read pairs (paired), then stop — a "
             "native head, so subsampling no longer needs an external `zcat | head | gzip`. "
             "0 = whole file."),
    two_pass: bool = typer.Option(
        False, "--two-pass/--one-pass",
        help="Shortlist ONE VxJ scaffold per read from a cheap segment search, then align only that one, instead of searching all 15,414 scaffolds. Reads it cannot resolve are realigned against the full reference, so none are dropped (measured: 0 lost at or above --min-score). **Reach for it when reads SPAN V INTO J, not by library type.** It needs a read to hit both a V and a J segment, and `fast_fraction` in the report is the predictor -- not whether the library is amplicon. Measured across regimes: 3.51x on a 48%-receptor TCR amplicon (fast path 85%), 2.96x on a 100%-receptor human TRB set (95.6%), 2.64x on mouse TRA (89%), but 1.03x SLOWER on a human IGH set of the SAME 100%-receptor data (16.3%, because those reads cover V and stop short of the short IGHJ target) and 0.762x on a 2.74%-receptor bulk library (5%), where the segment search is overhead on top of a rescue that is nearly the whole set. Needs `arda build-index`."),
    prefilter: bool = typer.Option(
        False, "--prefilter/--no-prefilter",
        help="Drop reads sharing no exact 16-mer with the reference BEFORE they reach MMseqs2. "
             "**Bulk only, and it is the bulk lever.** On a 0.024%-receptor library MMseqs2 "
             "spends 48.9s of a 82.6s run proving 4M reads are not receptor reads; the fitted "
             "cost model (wall ~ reads/46,353 + hits/350) says the dominant term is the READ "
             "COUNT, not the answer. Unlike MMseqs2's own prefilter this runs before `createdb`, "
             "so the FASTA write and DB build are skipped too. Costs ~0.5% of real reads "
             "(concentrated in J->C and hypermutated IGH), which is why it is OFF by default. "
             "Pointless on amplicon (46-49% receptor: almost everything passes)."),
    fast_segments: bool = typer.Option(
        False, "--fast-segments/--mmseqs-segments",
        help="With --two-pass, answer the segment pass structurally instead of with an MMseqs2 "
             "search. That pass exists only to learn each read's best V and best J with "
             "coordinates against a fixed 236kb germline reference -- a structural question, not a "
             "homology search. Measured on 100k amplicon reads: 74ms against MMseqs2's 2770ms "
             "(37x), agreeing with it on .9997 of V alleles and .9998 of J. It only NOMINATES: "
             "every candidate is still aligned against the full V+pad+J scaffold and scored by "
             "MMseqs2, so the AIRR output should not move. EXPERIMENTAL and off by default until "
             "that is proven end to end. Ignored without --two-pass."),
    indel_rescue: bool = typer.Option(
        False, "--indel-rescue/--no-indel-rescue",
        help="With --fast-segments, route reads that look like they carry an indel to the GAPPED "
             "rescue path instead of deciding them on the fast path. An ungapped extension follows "
             "ONE diagonal, so an indel-bearing read scores only up to the indel and its two "
             "halves land on two diagonals of the same target -- a signature visible in the seed "
             "votes before any extension runs. Measured on 341,294 real IGH mates: 3.18% of reads "
             "carry a V indel, and the rate tracks SHM load (0.74% at >=98% V identity, 8.00% "
             "below 90%), because AID makes indels and not only substitutions. These reads are "
             "REROUTED, never dropped, so a false positive costs a little speed and cannot cost a "
             "read. Ignored without --fast-segments."),
    segment_only_v: bool = typer.Option(
        False, "--v-only-on-segment/--no-v-only-on-segment",
        help="With --two-pass, align a read that hit a V but NO J against its own V segment "
             "instead of the full 15,414-scaffold reference. A `v_only` read carries no J -- that "
             "is the class, not a search failure -- so the scaffold search asks a question the "
             "read cannot answer, and it is 77%% of the amplicon rescue set at 338us/read against "
             "31us for a named-target alignment. MMseqs2 still does the alignment and still "
             "produces a real bit score, over exactly the nucleotides a whole-scaffold alignment "
             "of a J-less read would have covered, so --min-score keeps its meaning. Anything "
             "that fails falls through to the full-reference rescue: no read is lost. "
             "EXPERIMENTAL and off by default. Ignored without --two-pass."),
    adaptive: bool = typer.Option(
        False, "--adaptive/--no-adaptive",
        help="Cap alignments per read at --max-accept 40, then re-search UNCAPPED only the reads "
             "whose capped best score is under 90 bits. Attacks the align term, which is what is "
             "left once --prefilter has removed the scan term: on a 0.78%-receptor bulk library "
             "the search is 7.35s of a 12.25s map. Measured 2.17x on 1M bulk reads with zero "
             "reads lost. ⚠ Read survival is NOT the whole guarantee: on the real-read fixture it "
             "also moves junction_aa on 3 of 453 reads, two of them scoring 128 and 131 — far "
             "above the trigger — so a high score does not certify the best alignment was found. "
             "Opt in only where a junction-level difference is acceptable."),
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
                     limit=(limit or None), two_pass=two_pass, prefilter=prefilter,
                     fast_segments=fast_segments,
                     indel_rescue=indel_rescue,
                     segment_only_v=segment_only_v,
                     adaptive=adaptive,
                     emit_reads=emit_reads, report_path=report)
    typer.echo(
        f"[arda] {rep.mapped_reads}/{rep.total_reads} reads mapped "
        f"({rep.mapped_fraction * 100:.2f}%) in {rep.wall_seconds:.1f}s "
        f"({rep.reads_per_second:.0f} reads/s); loci={rep.per_locus}")


@rnaseq_app.command("correct")
def rnaseq_correct(
    input: Path = typer.Option(..., "--input", "-i", help="Mapped-reads AIRR TSV (from `map`)."),
    output: Path = typer.Option(..., "--output", "-o", help="Corrected clonotype table TSV."),
    max_subs: int = typer.Option(
        3, help="Max substitutions between an error child and its parent. A search RADIUS, not a "
                "threshold -- the length-scaled probability model still decides. Saturates at 3."),
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
    limit: int = typer.Option(
        0, "--limit", "-n",
        help="Analyse only the first N reads (single-end) / read pairs (paired), then stop. "
             "0 = whole file."),
    two_pass: bool = typer.Option(
        False, "--two-pass/--one-pass",
        help="Shortlist ONE VxJ scaffold per read from a cheap segment search, then align only that one, instead of searching all 15,414 scaffolds. Reads it cannot resolve are realigned against the full reference, so none are dropped (measured: 0 lost at or above --min-score). **Reach for it when reads SPAN V INTO J, not by library type.** It needs a read to hit both a V and a J segment, and `fast_fraction` in the report is the predictor -- not whether the library is amplicon. Measured across regimes: 3.51x on a 48%-receptor TCR amplicon (fast path 85%), 2.96x on a 100%-receptor human TRB set (95.6%), 2.64x on mouse TRA (89%), but 1.03x SLOWER on a human IGH set of the SAME 100%-receptor data (16.3%, because those reads cover V and stop short of the short IGHJ target) and 0.762x on a 2.74%-receptor bulk library (5%), where the segment search is overhead on top of a rescue that is nearly the whole set. Needs `arda build-index`."),
) -> None:
    """One-shot RNA-seq -> clonotypes for pipeline integration: ``map`` -> ``assemble`` -> ``correct``.

    Runs the three RNA-seq stages with the shipped defaults, writing under ``--out-dir``:

    * ``<prefix>.airr.tsv``            -- Stage-1 mapped reads (AIRR Rearrangement)
    * ``<prefix>.assembled.airr.tsv``  -- Stage-3 assembled long-CDR3 reads (if ``--assemble``)
    * ``<prefix>.clones.tsv``          -- corrected clonotype table (folds in the assembled clones)
    * ``<prefix>.arda.json``           -- merged run report

    Use the ``map`` / ``assemble`` / ``correct`` commands separately to tune their knobs.
    """
    from .rnaseq import pipeline

    # The body lives in `rnaseq.pipeline` so that `arda rnaseq reduce` -- the tail of a sharded
    # run -- calls the SAME Stage-2/3 function, not a copy of it. Two copies would drift, and
    # then "accuracy does not differ between run modes" would be a hope rather than a property.
    pipeline.run(r1, out_dir, out_prefix, r2=r2, organism=organism, threads=threads,
                 reconstruct=reconstruct, min_score=min_score,
                 kmer=(None if kmer == 0 else kmer), assemble=assemble,
                 complete_only=complete_only, map_d=map_d, limit=(limit or None),
                 two_pass=two_pass, echo=typer.echo)
    typer.echo(f"[arda] wrote {out_dir / f'{out_prefix}.airr.tsv'}, "
               f"{out_dir / f'{out_prefix}.clones.tsv'}, "
               f"{out_dir / f'{out_prefix}.arda.json'}")


@rnaseq_app.command("split")
def rnaseq_split(
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    out_dir: Path = typer.Option(..., "--out-dir", "-d", help="Directory for the shard FASTQs."),
    shards: int = typer.Option(..., "--shards", help="Number of contiguous blocks."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
) -> None:
    """Split paired FASTQ into contiguous blocks of read pairs, for a sharded Stage 1.

    Unlike ``arda split`` (FASTA/amplicon) this keeps the quality strings and never separates a
    fragment's two mates. Blocks are contiguous, so concatenating the per-shard AIRR in shard
    order reproduces the single-node row order exactly.
    """
    from .cluster import split_pairs

    written = split_pairs(r1, out_dir, shards=shards, r2=r2)
    typer.echo(f"[arda] wrote {len(written)} shard(s) to {out_dir}")


@rnaseq_app.command("reduce")
def rnaseq_reduce(
    shard_dir: Path = typer.Option(..., "--shard-dir", help="Directory of per-shard AIRR TSVs."),
    out_dir: Path = typer.Option(Path("."), "--out-dir", "-d", help="Directory for the outputs."),
    out_prefix: str = typer.Option(..., "--out-prefix", "-p", help="Output basename."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="mmseqs threads (0 = all cores)."),
    assemble: bool = typer.Option(True, "--assemble/--no-assemble", help="Stage 3."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions", help="Keep only complete junctions."),
    map_d: bool = typer.Option(True, "--map-d/--no-map-d", help="Map D segments."),
) -> None:
    """Merge a sharded Stage 1, then run Stages 2-3 ONCE over the whole thing.

    `assemble` and `correct` are global: run per shard, a clone split across shards is counted
    once per shard and contigs that tile across shards are never built. So the sharded path
    distributes only ``map``, and this is the step that finishes the job.
    """
    from .rnaseq import pipeline

    pipeline.reduce(shard_dir, out_dir, out_prefix, organism=organism, threads=threads,
                    assemble=assemble, complete_only=complete_only, map_d=map_d,
                    echo=typer.echo)
    typer.echo(f"[arda] wrote {out_dir / f'{out_prefix}.clones.tsv'}")


@rnaseq_app.command("slurm")
def rnaseq_slurm(
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    out_prefix: str = typer.Option(..., "--out-prefix", "-p", help="Output basename."),
    shards: int = typer.Option(..., "--shards", help="Array size = number of Stage-1 shards."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
    out_dir: Path = typer.Option(Path("."), "--out-dir", "-d", help="Directory for the outputs."),
    work_dir: Path = typer.Option(Path("arda_slurm"), help="Scratch for shards + submit.sh."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(8, help="--cpus-per-task, and mmseqs threads."),
    kmer: int = typer.Option(12, "--kmer", "-k", help="MMseqs2 -k (the memory knob)."),
    min_score: float = typer.Option(75.0, "--min-score", help="Min bit score to keep a read."),
    reconstruct: bool = typer.Option(False, "--reconstruct", help="Merge overlapping mates."),
    assemble: bool = typer.Option(True, "--assemble/--no-assemble", help="Stage 3."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions", help="Keep only complete junctions."),
    map_d: bool = typer.Option(True, "--map-d/--no-map-d", help="Map D segments."),
    partition: Optional[str] = typer.Option(None, help="SLURM partition."),
    time_limit: str = typer.Option("04:00:00", "--time", help="Walltime per array task."),
    mem: str = typer.Option("8G", help="Memory per array task (Stage 1 is flat, ~300-400 MB)."),
    reduce_time: str = typer.Option("08:00:00", help="Walltime for the reduce step."),
    reduce_mem: str = typer.Option(
        "16G", help="Memory for reduce. Stage 3 holds the clone set: budget ~4 GB, more for a "
                    "B-cell-rich sample (2.7 GB measured at 28k clonotypes)."),
    submit: bool = typer.Option(False, "--submit", help="Run the generated script."),
) -> None:
    """Write (and optionally submit) a SLURM script: split → array-``map`` → reduce.

    Only Stage 1 is distributed; Stages 2-3 run once over the merged AIRR, through the same
    code path ``arda rnaseq run`` uses. With contiguous pair shards the result is
    byte-identical to a single-node run.
    """
    import os
    import subprocess

    from .cluster import render_rnaseq_submit_script

    work_dir.mkdir(parents=True, exist_ok=True)
    script = work_dir / "submit.sh"
    script.write_text(render_rnaseq_submit_script(
        r1, out_prefix, work_dir, shards=shards, r2=r2, out_dir=out_dir, organism=organism,
        threads=threads, kmer=kmer, min_score=min_score, reconstruct=reconstruct,
        assemble=assemble, complete_only=complete_only, map_d=map_d, partition=partition,
        time=time_limit, mem=mem, reduce_time=reduce_time, reduce_mem=reduce_mem,
        arda_mmseqs=os.environ.get("ARDA_MMSEQS")))
    script.chmod(0o755)
    typer.echo(f"[arda] wrote {script}")
    if submit:
        subprocess.run(["bash", str(script)], check=True)


if __name__ == "__main__":
    app()
