"""arda command-line interface.

**Three MODES, named after the library.** Each owns the speed configuration that is right for its
regime, because the two configurations do not compose and picking the wrong one is not an error —
it is a silent 2-4x slowdown:

* ``arda rnaseq``     — bulk / whole-transcriptome RNA-seq  (``--prefilter``)
* ``arda amplicon``   — targeted RepSeq / 5'RACE  (``--two-pass --fast-segments --v-only-on-segment``)
* ``arda singlecell`` — reserved; not implemented yet

⛔ Until 2.16.0 the only entry point was ``arda rnaseq run``, which was used for amplicon too and
exposed the regime as four loose flags. ``--two-pass`` ALONE is a LOSS on both regimes (0.762x on
bulk, 0.87x on an IGH amplicon), so the one combination that was easy to reach was the dominated
one. Naming the mode is the fix; ``--exact`` opts out of every speedup.

The pipeline STAGES are separate commands, so any of them can be run, inspected or replaced on its
own: ``arda map`` -> ``arda assemble`` -> ``arda correct``, plus ``arda shm`` (SHM recount). Other
commands: ``arda annotate`` (FASTA/FASTQ -> AIRR), ``arda markup``, ``arda resolve-ties``,
``arda cluster`` (shard/submit), ``arda build-db`` / ``build-index`` / ``export-ref``,
``arda igblast``, ``arda info``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from . import __version__
from ._log import logger as log

app = typer.Typer(add_completion=False, help="Antigen Receptor Domain Annotation")

# Shared by every command that maps D. The gate is an E-value: the expected number of chance
# matches at least this good given the interior length and the D-database size (see
# `annotate.transfer._map_d`). 0.2 is the shipped operating point and stays the default -- this
# option only lets a caller ASK for the strict band, which is where D is worth trusting.
# ⛔ The default is None, NOT 0.2. The shipped operating point is alphabet-dependent -- 0.2 for
# nt, 0.05 for aa -- so a literal 0.2 here would silently LOOSEN `--seqtype aa` by 4x while
# looking like a no-op. None means "whatever this alphabet's calibrated value is".
_D_EVALUE_HELP = (
    "Max E-value for a D call (and for a tandem second D). Lower is stricter. Default: the "
    "calibrated operating point, 0.2 for nt and 0.05 for aa. Measured against IgBLAST at gene "
    "level on nt: 0.2 agrees .9765 on a TRB amplicon and .9417 on bulk IGH; 0.01 agrees .9985 "
    "and 1.0000, at roughly a third of the call rate.")

_EXACT_HELP = (
    "Turn OFF every speedup this mode would enable and run the shipped one-pass path. Use it when "
    "you want the exact output the default aligner produces -- notably to avoid --prefilter's "
    "~0.15 %% read cost on a bulk library, or to A/B a mode's preset against it.")

_SHM_HELP = (
    "SHM scoping. `framework` (default) keeps v_identity / v_mutations / j_mutations to positions "
    "OUTSIDE the junction, using the germline anchors arda emits per read. ⛔ Segment scoping "
    "alone is not junction exclusion -- the V germline's 3' tail and the J germline's 5' head are "
    "inside the junction, so chew-back and N/P bases used to enter both lists: measured on a TRA "
    "amplicon, where TCRs cannot hypermutate so every entry is spurious, 1.046 V and 1.658 J "
    "entries per read, 86.2 %% of the J ones at germline position <= 10. `both` also emits the old "
    "junction-inclusive values as v_identity_full / v_mutations_full / j_mutations_full. `off` "
    "emits no SHM fields. IGH/IGK/IGL are where SHM is real; on TR the surviving entries are "
    "allele mismatches in the templated framework, not hypermutation.")

_COMPLETE_JUNCTION_HELP = (
    "Finish a junction whose read reached Cys104 but stopped before [FW]118, taking up to N nt "
    "from the called J's germline. 0 (default) emits observed junctions only. The J's 5' chew-back "
    "and the N/P additions are all UPSTREAM of the read's last aligned J base, so what is missing "
    "is germline-TEMPLATED -- unlike the V side, where a short read is missing bases nothing "
    "templates. ⛔ The added bases are IMPUTED, not observed: every completed row carries the count "
    "in `junction_completed_nt`, so filter or weight on that column rather than trusting the "
    "junction. ⚠ On IG the imputed span can hide the SHM the read would have shown, biasing a "
    "completed junction's 3' end toward germline; TR does not hypermutate and has no such cost.")

_CHIMERA_HELP = (
    "Also emit `chimera_parents`: for each clonotype, two MORE ABUNDANT clonotypes of the same "
    "locus that explain it as prefix+suffix across one breakpoint -- the PCR template-switch "
    "signature. ⛔ FLAG ONLY, never a filter: measured 0.40 % of clonotypes / 0.18 % of reads on "
    "bulk RNA-seq (IG) against 0.01 % on a TRA amplicon, a 20x enrichment in the direction "
    "template-switch chemistry predicts but far too small to justify deleting clonotypes -- and the "
    "signature cannot separate a true chimera from two real clones sharing a prefix and a suffix. "
    "⛔ The breakpoint must sit in the NON-TEMPLATED core: a junction is V 3' tail + N/P/D + J 5' "
    "head, both tails germline, so the same test run on the raw junction calls 52 % of clonotypes "
    "chimeric. Requires the reference (no anchors -> no flags, never a germline-driven guess).")

_CALL_LEVEL_HELP = (
    "At what resolution a V/J call names a germline, BEFORE the clonotype key is formed. `allele` "
    "(default) = TRGJ1*01. `gene` = TRGJ1, which collapses allele-level CALL SPLITS -- Jurkat "
    "carries TRGJ1*01 at 64 reads against TRGJ1*02 at 140 on the same junction, and no error model "
    "can see that because two identical junctions have no discriminating base. It also collapses a "
    "tie list whose members differ only by allele. ⚠ IGH carries 4.33 alleles/gene, ~2x every "
    "other locus, so measure the cost on your own library before quoting one.")

_ISOTYPE_HELP = (
    "Resolve the IGH isotype: `c_call` (CH1 exon) and `c_class` (IGHG / IGHM / IGHA / ...) per "
    "read, voted once per FRAGMENT, then per clonotype. Report the CLASS, never the subclass -- "
    "IGHG1-4 are ~95 %% identical over CH1, so the top gene ties on 26.7 %% of real reads and the "
    "top class never does. --no-isotype drops the constant-region columns. IGH only in practice: "
    "TRA/TRD/IGK have one C allele and IGL's seven are one class.")



def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def _main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True,
        help="Show the arda version and exit."),
    verbose: int = typer.Option(
        0, "-v", "--verbose", count=True,
        help="Raise console verbosity. Default prints the stage and progress lines; -v adds "
             "DEBUG (per-chunk accounting, reference and index decisions, the fallbacks a run "
             "took silently) with the level and module on each line. Repeatable."),
    quiet: bool = typer.Option(
        False, "-q", "--quiet",
        help="Console warnings and errors only. Does NOT quieten --log-file, which stays at "
             "DEBUG -- so a cluster job can be silent and still leave a full record."),
    log_file: Optional[Path] = typer.Option(
        None, "--log-file",
        help="Also write a DEBUG log here, with a timestamp and the process peak RSS on every "
             "line. This is the artifact to keep from a long run: re-running a 10-hour bulk "
             "sample because the console was at the default level is not a diagnosis."),
) -> None:
    """Antigen Receptor Domain Annotation.

    Progress goes to **stderr** and results to **stdout**, so ``arda export-ref ... > out.tsv``
    and ``arda ... | head`` are safe. ``-v`` / ``-q`` / ``--log-file`` go BEFORE the subcommand:
    ``arda -v --log-file run.log amplicon --r1 ...``.
    """
    from ._log import setup

    setup(verbosity=verbose, quiet=quiet, log_file=log_file)


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
    allow_chimeras: bool = typer.Option(
        False, "--allow-chimeras",
        help="Also build TRDV x TRAJ scaffolds, which the default reference refuses as chimeric. "
             "TRDV1/2/3 are dedicated delta V genes, but they sit INSIDE the TRA locus between the "
             "TRAV genes and the TRAJ cluster. Measured on 48,030 TRA amplicon reads: IgBLAST "
             "calls TRDV1 + a TRAJ on 530 of them (1.10 % of the library, median v_score 93.8, all "
             "carrying a junction) and MiXCR independently agrees, while arda calls the same J as "
             "both and emits NO v_call -- 83 % of its whole remaining v_gene gap. Either the "
             "pairing is real and the default drops 1.1 % of a TRA repertoire, or it is a chimera "
             "and the other two report it because they call V and J independently. That is a "
             "domain judgement, so it is a flag. Off by default; assert the SCAFFOLD COUNT after "
             "building -- an earlier attempt added 7 scaffolds, not the ~483 the product implies.",
    ),
) -> None:
    """Build the curated reference database (Phase 1)."""
    from .refbuild.build import build

    build(organism, one_allele_per_gene=one_allele_per_gene, allow_chimeras=allow_chimeras)


@app.command("export-ref")
def export_ref(
    organism: str = typer.Option("human", help="Reference organism."),
    kind: str = typer.Option(
        "scaffolds", "--kind",
        help="scaffolds = the V x J (and J+C) reference the mapper aligns against; "
             "segments = the collapsed per-allele V/J/C reference the two-pass search uses; "
             "anchors = the per-allele CDR3 anchor table (germline_nt, templated_aa, status)."),
    fmt: str = typer.Option(
        "tsv", "--format",
        help="tsv (sequence + every region as its own column), fasta, gff3 (regions as features; "
             "GFF3 is 1-based closed like arda, so coordinates pass through unchanged), or airr "
             "(the same rows shaped as an AIRR Rearrangement, so a scaffold can be fed straight "
             "into anything that reads arda's own output)."),
    seqtype: str = typer.Option("nt", help="'nt' or 'aa' reference."),
    locus: Optional[str] = typer.Option(
        None, "--locus", help="Comma-separated loci to keep (e.g. 'TRB,IGH'). Default: all."),
    output: Optional[Path] = typer.Option(None, "-o", "--output", help="Write here (default: stdout)."),
) -> None:
    """Export reference sequences with their FR/CDR markup.

    The reference is arda's most valuable offline artifact -- every in-frame V.J germline scaffold
    with IgBLAST-quality FR1-4 / CDR1-3 coordinates -- and until now it was only reachable by
    hand-joining the build's TSVs against its FASTAs. Coordinates are **1-based closed** (AIRR).
    """
    from .refexport import export_reference

    loci = {x.strip() for x in locus.split(",") if x.strip()} if locus else None
    n = export_reference(organism, kind=kind, fmt=fmt, seqtype=seqtype, loci=loci, out=output)
    typer.echo(f"[arda] exported {n} {kind} record(s) as {fmt}"
               + (f" -> {output}" if output else ""), err=True)


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
    d_max_evalue: Optional[float] = typer.Option(
        None, "--d-max-evalue", help=_D_EVALUE_HELP),
) -> None:
    """Annotate FR/CDR regions and write an AIRR TSV (streamed, memory-bounded)."""
    from .annotate.mapper import annotate_file

    annotate_file(input, output, organism=organism, seqtype=seqtype,
                  threads=threads, strand=strand, chunk_size=chunk_size, map_d=map_d,
                  d_max_evalue=d_max_evalue)


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

    df = pl.read_csv(input, separator="\t", infer_schema_length=0, quote_char=None)
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


# ── CLUSTER ───────────────────────────────────────────────────────────────────────────────────
# One group for every sharded/SLURM helper. They used to be split across the top level (FASTA) and
# the `rnaseq` group (paired FASTQ) under near-identical names, which is exactly the shape that
# lets someone shard a paired library with the round-robin FASTA splitter and separate the mates.
cluster_app = typer.Typer(add_completion=False,
                          help="Shard a run across a cluster (split -> array map -> reduce).")
app.add_typer(cluster_app, name="cluster")


@cluster_app.command("split-fasta")
def split_fasta(
    input: Path = typer.Argument(..., help="Input FASTA/FASTQ."),
    out_dir: Path = typer.Argument(..., help="Directory for shard FASTA files."),
    shards: int = typer.Option(..., "--shards", help="Number of shards."),
) -> None:
    """Round-robin split an input into N shard FASTA files (amplicon / single-end).

    **Not for paired RNA-seq.** This writes FASTA, so quality is dropped, and it round-robins
    *records*, so a fragment's two mates land in different shards. Use ``arda cluster split``.
    """
    from .cluster import split as _split

    paths = _split(input, out_dir, shards)
    typer.echo(f"wrote {len(paths)} shards to {out_dir}")


@cluster_app.command("merge")
def merge(
    shard_dir: Path = typer.Argument(..., help="Directory of per-shard AIRR TSVs."),
    output: Path = typer.Argument(..., help="Combined AIRR TSV."),
) -> None:
    """Concatenate per-shard AIRR TSVs into one (single header)."""
    from .cluster import merge as _merge

    _merge(shard_dir, output)
    typer.echo(f"merged -> {output}")


@cluster_app.command("submit-fasta")
def slurm_fasta(
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
    end. Use ``arda cluster submit``.
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


# ── STAGES ────────────────────────────────────────────────────────────────────────────────────
# Flat, not under a `rnaseq` group: `map`/`assemble`/`correct` are the same three stages whatever
# the library is, and hiding them under one regime's name is what made `arda rnaseq run` the
# amplicon entry point in the first place.


@app.command("map")
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
    d_max_evalue: Optional[float] = typer.Option(
        None, "--d-max-evalue", help=_D_EVALUE_HELP),
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
    junction_quality: bool = typer.Option(
        False, "--junction-quality/--no-junction-quality",
        help="Also emit a `junction_quality` column: the read's Phred+33 string over exactly the "
             "bases of `junction`, same orientation. Stage 1 is the ONLY place the FASTQ quality "
             "is still in hand, and `correct --min-junction-q` needs it -- measured on a MIGEC "
             "spike-in library, the mismatching base of a real variant reads median Q 34-35 "
             "against median Q 16 for the sequencing-error cloud around it. OFF by default: it "
             "appends a non-schema column, so the default output is unchanged. Not usable with "
             "--reconstruct (a merged fragment has no single input quality string)."),
    mutation_quality: bool = typer.Option(
        False, "--mutation-quality/--no-mutation-quality",
        help="Also emit `v_mutation_quality` / `j_mutation_quality`: the read's Phred score at "
             "each entry of `v_mutations` / `j_mutations`, comma-joined, one-for-one and in the "
             "same order. A NOVEL ALLELE, somatic hypermutation and a base-miscall are the same "
             "string in the mutation list -- what separates them is how often the mutation recurs "
             "across an allele's reads and how good the base is, and only the second of those "
             "needs a column. `arda stats` reads it for `allele_candidate mean_quality`. OFF by "
             "default (non-schema columns, and it re-walks the alignment in Python); not usable "
             "with --reconstruct."),
    shm: str = typer.Option("framework", "--shm", help=_SHM_HELP),
    complete_junction_nt: int = typer.Option(
        0, "--complete-junctions", help=_COMPLETE_JUNCTION_HELP),
    emit_reads: Optional[Path] = typer.Option(
        None, "--emit-reads", help="Also write the mapped reads as FASTA (for handoff)."),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Stage 1 — filter + map receptor reads; write only the mapped ones (AIRR TSV)."""
    from .rnaseq.map import map_rnaseq

    map_rnaseq(r1, output, r2=r2, organism=organism, threads=threads,
                    sensitivity=sensitivity, strand=strand, chunk_size=chunk_size,
                    map_d=map_d, d_max_evalue=d_max_evalue,
                    reconstruct=reconstruct, min_score=min_score,
                    max_seqs=max_seqs, kmer=(None if kmer == 0 else kmer),
                    drop_constant_only=drop_constant_only,
                    limit=(limit or None), two_pass=two_pass, prefilter=prefilter,
                    fast_segments=fast_segments,
                    indel_rescue=indel_rescue,
                    segment_only_v=segment_only_v,
                    adaptive=adaptive, with_junction_quality=junction_quality,
                    with_mutation_quality=mutation_quality, shm=shm,
                    complete_junction_nt=complete_junction_nt,
                    emit_reads=emit_reads, report_path=report)
    # `map_rnaseq` already logged the summary line; only the output path belongs on stdout.
    typer.echo(str(output))


@app.command("correct")
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
    error_method: Optional[str] = typer.Option(
        None, help="simple = spanning-read counts; binom|betabinom = per-position read-depth "
                   "pileup for very low coverage. Default: whatever --ec-mode selects (simple). "
                   "⛔ binom/betabinom are ~270x slower AND more aggressive on a deep library "
                   "(MIGEC 302k reads: 0.73s/143 clonotypes vs 197s/79 and 254s/78) and "
                   "byte-identical on a monoclonal one -- neither is in a mode for that reason."),
    ec_mode: str = typer.Option(
        "fast", "--ec-mode",
        help="Denoising preset. `fast` (default) = shipped behaviour: the abundance error model on "
             "spanning read counts, no quality gate. `accurate` = the same plus --min-junction-q "
             "20, judging the one base that discriminates a clonotype from its parent on its Phred "
             "score. `amplicon` / `rnaseq` add the QUALITY-DIRECTED RESCUE: a clonotype 4+ subs "
             "from anything has no ladder of observed intermediates behind it (0 of 13 at k=4 on "
             "Jurkat vs 0.0019 predicted) and no discriminating base to gate on, but its reads are "
             "measurably bad -- those are routed to a much more abundant parent, NEVER deleted, and "
             "one with no parent keeps its reads. `amplicon` searches wide (12 subs, 50x ratio) "
             "because a real clonotype there is deep; `rnaseq` stays narrow (6 subs, 200x) because "
             "singletons are the norm in a 0.02-3 % receptor library. All need "
             "`map --junction-quality`. An explicit --error-method / --min-junction-q wins."),
    clonotype_key: str = typer.Option(
        "full", "--clonotype-key",
        help="`full` (default) = (locus, v_call, j_call, junction). `junction` = (locus, junction): "
             "V/J are canonicalised to the junction's majority first, so CALL SPLITS collapse -- a "
             "junction byte-identical to an abundant clone's under a different V or J call, which "
             "no error model can see because there is no discriminating base. On Jurkat that is "
             "the largest error class by reads (130 of 14,531, incl. an allele-level TRG split): "
             "TRB 35 -> 33 clonotypes at purity .99096 -> .99696, reads unchanged. Measured cost "
             "on a polyclonal TRA amplicon: 132 of 19,956 clonotypes merge (0.66 %), and the "
             "minority call there carries 1 read against 4-10 on a short junction."),
    call_level: str = typer.Option("allele", "--call-level", help=_CALL_LEVEL_HELP),
    flag_chimeras: bool = typer.Option(False, "--flag-chimeras", help=_CHIMERA_HELP),
    isotype: bool = typer.Option(True, "--isotype/--no-isotype", help=_ISOTYPE_HELP),
    min_junction_q: Optional[int] = typer.Option(
        None, "--min-junction-q",
        help="Reassign a read whose junction differs from its putative parent at ANY base below "
             "this Phred score ONTO that parent; matching bases are not evidence and are not "
             "looked at. The read is never discarded -- it came off a real rearrangement and is "
             "counted in the parent clonotype. 0 = off. "
             "Requires the `junction_quality` column (`map --junction-quality`) and RAISES "
             "without it rather than silently not gating. Measured on the MIGEC spike-ins, a Q20 "
             "gate takes the published-variant-to-error-cloud read ratio from 1.349 to 2.110, "
             "keeping 86 % of the real variant's reads and removing 89 % of the distinct error "
             "clonotypes; it plateaus over Q20-32 and starts eating real variants by Q35."),
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
        help="Add the D columns (d_call/d2_call/d_support/d2_support, the D and V/J spans, and "
             "np1-np3) to each clonotype, mapped into its error-corrected junction (once per "
             "clonotype, not per read). Coordinates are 1-based closed in JUNCTION space."),
    d_max_evalue: Optional[float] = typer.Option(
        None, "--d-max-evalue", help=_D_EVALUE_HELP),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Stage 2 — collapse CDR3 sequencing errors into clonotypes (per-substitution/indel model)."""
    from .rnaseq.correct import correct_airr

    rep = correct_airr(input, output, organism=organism, map_d=map_d,
                       d_max_evalue=d_max_evalue, max_subs=max_subs, max_indel=max_indel, error_rate=error_rate,
                       indel_rate=indel_rate, require_vj=require_vj, error_method=error_method,
                       ec_mode=ec_mode, min_junction_q=min_junction_q,
                       clonotype_key=clonotype_key, call_level=call_level, isotype=isotype,
                       flag_chimeras=flag_chimeras,
                       complete_only=complete_only, read_map=read_map, extra_airr=extra_airr,
                       report_path=report)
    log.info(
        f"correct: {rep.clonotypes_in} -> {rep.clonotypes_out} clonotypes "
        f"({rep.collapsed} collapsed) over {rep.reads} reads"
        + (f"; dropped {rep.reads_incomplete}/{rep.reads_with_junction} incomplete junctions"
           if rep.reads_incomplete else "")
        + (f"; moved {rep.reads_low_quality} reads onto their parent on --min-junction-q "
           f"({rep.clonotypes_low_quality} clonotypes emptied)"
           if rep.reads_low_quality else ""))
    typer.echo(str(output))


@app.command("assemble")
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
    d_max_evalue: Optional[float] = typer.Option(
        None, "--d-max-evalue", help=_D_EVALUE_HELP),
    report: Optional[Path] = typer.Option(None, "--report", help="Write a JSON run report."),
) -> None:
    """Stage 3 — assemble long-CDR3 contigs the reads don't individually span.

    Reconstructs clonotypes (V(D)J ultralong CDR3s, ~20-40 aa) that ``map`` filters but that no
    single 100-150 bp read spans, by anchored greedy overlap-extension over the mapped reads.
    Feed the result to ``correct --extra-airr`` (the modes do this automatically).
    """
    from .rnaseq.assemble import assemble_contigs

    rep = assemble_contigs(input, output, organism=organism, threads=threads, map_d=map_d,
                           d_max_evalue=d_max_evalue, report_path=report)
    log.info("assemble: %d/%d complete contigs from %d seeds; rescued %d reads",
             rep.contigs_complete, rep.contigs, rep.seeds, rep.reads_rescued)
    typer.echo(str(output))


@app.command("stats")
def stats_cmd(
    output: Path = typer.Option(..., "--output", "-o", help="QC TSV ('-' for stdout)."),
    airr: Optional[Path] = typer.Option(
        None, "--airr", "-i",
        help="Mapped-reads AIRR TSV (`map` or `annotate`): the per-read and per-chain rows, "
             "the per-gene read counts, and the candidate-allele shortlist."),
    clones: Optional[Path] = typer.Option(
        None, "--clones", "-c",
        help="Clonotype table (`correct`): the per-chain clonotype rows, chimera counts and "
             "per-gene clonotype counts."),
    report: Optional[Path] = typer.Option(
        None, "--report", "-r",
        help="`<prefix>.arda.json` (or a single-stage `--report` JSON). This is where total and "
             "mapped reads, threads, wall time and peak RSS come from -- the AIRR holds only the "
             "reads that mapped, so nothing in it can recover them."),
    r1: Optional[Path] = typer.Option(
        None, "--r1", help="Input FASTQ, read ONLY for its size on disk and to record that the "
                           "library is paired. Neither is recoverable from the AIRR."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ (marks the library paired)."),
    organism: str = typer.Option(
        "human", help="Reference organism, used for the V/J gene universe behind the coverage "
                      "fractions."),
    allele_min_frac: float = typer.Option(
        0.5, "--allele-min-frac",
        help="A V mutation is a CANDIDATE ALLELE at or above this frequency among the reads "
             "calling that allele. A germline the reference does not carry is in essentially "
             "every read of its allele; somatic hypermutation is per-clone and does not reach "
             "half. ⚠ A shortlist to look at, never a genotype call -- arda does not genotype."),
    allele_min_reads: int = typer.Option(
        10, "--allele-min-reads",
        help="...and in at least this many reads, so a 2-read allele cannot mint a candidate off "
             "one read."),
) -> None:
    """Run QC — reads, clonotypes, chains, genes and candidate alleles as one long-format TSV.

    Reads only what a run already wrote, adds no alignment, and never filters anything. Four
    columns, ``scope``/``key``/``metric``/``value``, so a metric can be grepped, joined across
    samples, or plotted without reshaping:

    \b
      run               map / correct / assemble   the run report, flattened verbatim
      sample            (blank)                    library-wide totals and coverage
      chain             TRB, IGH, ...              per locus, reads AND clonotypes
      v_gene / j_gene   TRBV19                     reads and clonotypes per germline gene
      allele_candidate  TRBV19*01:G45A             a recurrent, high-quality V mutation

    Every input is optional and contributes its own scopes, so this works on a bare ``annotate``
    output as well as on a full run directory::

        arda stats -i s.airr.tsv -c s.clones.tsv -r s.arda.json --r1 s_1.fq.gz --r2 s_2.fq.gz \\
                   -o s.stats.tsv

    ⚠ ``junction_quality_mean`` needs ``map --junction-quality``, and the
    ``allele_candidate mean_quality`` / ``shm_variant_mean_quality`` rows need
    ``map --mutation-quality``. Rows that have no input are omitted, never emitted as 0.
    """
    import sys

    from .stats import collect, write_stats

    if airr is None and clones is None and report is None:
        raise typer.BadParameter("give at least one of --airr / --clones / --report")
    rows = collect(airr=airr, clones=clones, report=report, r1=r1, r2=r2, organism=organism,
                   allele_min_frac=allele_min_frac, allele_min_reads=allele_min_reads)
    write_stats(rows, sys.stdout if str(output) == "-" else output)
    log.info("stats: %d rows over %d scopes", len(rows), len({r[0] for r in rows}))
    if str(output) != "-":
        typer.echo(str(output))


@app.command("shm")
def shm_cmd(
    input: Path = typer.Option(..., "--input", "-i", help="AIRR TSV from `map` or `annotate`."),
    output: Path = typer.Option(..., "--output", "-o", help="AIRR TSV with rescoped SHM fields."),
    mode: str = typer.Option("framework", "--mode", help=_SHM_HELP),
) -> None:
    """Recount somatic hypermutation OUTSIDE the junction — IGH / IGK / IGL.

    ``v_mutations`` / ``j_mutations`` / ``v_identity`` are scoped BY SEGMENT, and that is not the
    same as being outside the junction: a rearranged junction is *V 3' tail + N/P + J 5' head*, so
    the templated tails of both germlines lie inside it and every chew-back or non-templated base
    there reads as a substitution against a germline that does not template it. arda 2.14.0
    documented a guarantee that this could not happen; it was wrong, and this is the retraction.

    ⛔ **Needs no reference and no re-map.** ``v_anchor_nt`` / ``j_anchor_nt`` and the alignment
    strings are already in the file, so a table written by arda 2.14.0 or later can be recounted in
    place. A file older than that has no anchor columns and this RAISES rather than copying the
    input through with a success message.

    ⚠ IGH/IGK/IGL are where SHM is real. Scoping runs on every locus anyway, because the defect is
    not IG-specific — measured on a **TRA** amplicon, where TCRs cannot hypermutate so every entry
    is spurious by construction: 1.046 V and 1.658 J entries per read.
    """
    from .shm import recount_airr

    rep = recount_airr(input, output, mode=mode)
    log.info("shm (%s): %d rows, %d -> %d mutation entries (%d junction-internal)",
             rep["mode"], rep["rows"], rep["mutations_in"], rep["mutations_out"], rep["removed"])
    typer.echo(str(output))


# ── MODES ─────────────────────────────────────────────────────────────────────────────────────
# ⛔ The regime is the COMMAND NAME, not a flag combination. `arda rnaseq run` used to be the only
# entry point and was used for amplicon too, with the regime spelled out as four loose flags that
# do NOT compose: `--two-pass --fast-segments --v-only-on-segment` is amplicon, `--prefilter` is
# bulk, and `--two-pass` alone -- the single flag `run` exposed for four releases -- is a LOSS on
# both (0.762x bulk, 0.87x IGH amplicon). Naming the mode makes the dominated combination
# unreachable by accident.

#: What each mode turns on, and why. Read as: (two_pass, fast_segments, segment_only_v, prefilter).
#:
#: `amplicon` -- primer-anchored reads span V INTO J, so ~85 % hit both a V and a J segment and the
#: cheap segment pass answers them structurally. Measured on IGH RepSeq at 32 threads:
#: 316.44 s -> 76.25 s (4.15x) and 4,018 -> 1,479 MB.
#:
#: `rnaseq` -- 0.02-3 % of reads are receptor-derived, so MMseqs2 spends most of the run proving
#: reads are NOT receptor reads; the 16-mer prefilter removes them before `createdb`. 1.99x.
#: ⚠ It costs ~0.15 % of mapped reads (122 bulk datasets, up to 2.46 % on one library),
#: concentrated in J->C and hypermutated IGH. `--exact` turns it off.
_MODE_SPEED = {
    "amplicon": {"two_pass": True, "fast_segments": True, "segment_only_v": True,
                 "prefilter": False},
    "rnaseq": {"two_pass": False, "fast_segments": False, "segment_only_v": False,
               "prefilter": True},
}

def _mode_run(mode: str, *, exact: bool, indel_rescue: bool, **kw) -> None:
    """Body shared by `arda rnaseq` and `arda amplicon`: resolve the preset, then run the pipeline.

    ⛔ ONE body, not one per mode. The mode commands and `arda cluster reduce` already share
    `pipeline.finish` for the same reason: two copies drift in a parameter, and then "the modes
    only differ in their preset" is a hope rather than a property.
    """
    from .rnaseq import pipeline

    speed = {k: False for k in _MODE_SPEED[mode]} if exact else dict(_MODE_SPEED[mode])
    if indel_rescue:
        # It needs the fast segment pass, which only the amplicon preset turns on. arda would
        # otherwise ignore the flag silently -- the failure this project keeps hitting.
        if not speed["fast_segments"]:
            raise typer.BadParameter(
                "--indel-rescue requires the fast segment pass, which only `arda amplicon` "
                "enables (and --exact disables). It is ignored in every other configuration.")
        speed["indel_rescue"] = True
    pipeline.run(**speed, **kw)
    out_dir, out_prefix = Path(kw["out_dir"]), kw["out_prefix"]
    # One path per line on stdout, so `arda amplicon ... | tail -1` and `$(...)` work. Everything
    # else this run said went to stderr through the logger.
    for name in ("airr", "clones", "report", "stats"):
        path = out_dir / pipeline.OUTPUTS[name].format(prefix=out_prefix)
        if path.exists():
            typer.echo(str(path))


@app.command("rnaseq")
def rnaseq_mode(
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
    out_prefix: str = typer.Option(
        ..., "--out-prefix", "-p",
        help="Output basename. Writes <prefix>.airr.tsv, <prefix>.clones.tsv, <prefix>.arda.json and <prefix>.stats.tsv."),
    out_dir: Path = typer.Option(Path("."), "--out-dir", "-d", help="Directory for the outputs."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="mmseqs threads (0 = all cores)."),
    assemble: bool = typer.Option(
        True, "--assemble/--no-assemble",
        help="Stage 3: assemble long-CDR3 contigs no single 100-150 bp read spans (V(D)J "
             "ultralong, ~20-40 aa) and fold them into the clonotype table. Recovers ~95%% of the "
             "abundant long clones a filter-only pass misses. --no-assemble keeps the run on the "
             "flat mapping-only memory profile."),
    shm: str = typer.Option("framework", "--shm", help=_SHM_HELP),
    complete_junction_nt: int = typer.Option(
        0, "--complete-junctions", help=_COMPLETE_JUNCTION_HELP),
    isotype: bool = typer.Option(True, "--isotype/--no-isotype", help=_ISOTYPE_HELP),
    call_level: str = typer.Option("allele", "--call-level", help=_CALL_LEVEL_HELP),
    map_d: bool = typer.Option(
        True, "--map-d/--no-map-d",
        help="Map D segments in all three stages. --no-map-d is ~19%% faster and yields an "
             "identical read set, v_call and junction -- only d_call is dropped."),
    d_max_evalue: Optional[float] = typer.Option(None, "--d-max-evalue", help=_D_EVALUE_HELP),
    ec_mode: str = typer.Option(
        "rnaseq", "--ec-mode",
        help="Denoising preset; `rnaseq` is this mode's default. `fast` = the abundance model "
             "only (arda's historical behaviour). `accurate` adds --min-junction-q 20. `rnaseq` "
             "adds the quality-directed rescue kept NARROW (6 subs, 200x ratio) because bulk "
             "RNA-seq singletons are mostly real. ⛔ Nothing in any mode discards a read: an "
             "orphan with no qualifying parent keeps its reads. The Stage-1 quality column the "
             "gates need is turned on automatically here."),
    min_junction_q: Optional[int] = typer.Option(
        None, "--min-junction-q",
        help="Explicit Phred floor for the discriminating base; overrides --ec-mode's preset."),
    clonotype_key: str = typer.Option(
        "full", "--clonotype-key",
        help="`full` (default) = (locus, v_call, j_call, junction). `junction` = (locus, junction) "
             "with V/J canonicalised to the junction's majority first."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions",
        help="Keep only complete junctions (C104..[FW]118, in frame, no stop) when forming "
             "clonotypes."),
    reconstruct: bool = typer.Option(
        False, "--reconstruct", help="Merge overlapping paired mates into one fragment."),
    min_score: float = typer.Option(
        75.0, "--min-score", help="Min MMseqs2 bit score to keep a mapped read (0 = keep all)."),
    kmer: int = typer.Option(
        12, "--kmer", "-k",
        help="MMseqs2 -k (the memory knob; k=12 ~= 298 MB peak RSS). 0 = MMseqs2 default (~8 GB)."),
    limit: int = typer.Option(
        0, "--limit", "-n",
        help="Analyse only the first N reads / read pairs, then stop. 0 = whole file."),
    exact: bool = typer.Option(False, "--exact", help=_EXACT_HELP),
) -> None:
    """BULK RNA-seq -> clonotypes. map -> assemble -> correct, with the bulk speed preset.

    For whole-transcriptome libraries, where 0.02-3 % of reads are receptor-derived and a read
    lands anywhere in a transcript rather than spanning V into J. Enables ``--prefilter`` (1.99x);
    ``--exact`` turns it off.

    Writes under ``--out-dir``:

    * ``<prefix>.airr.tsv``            -- Stage-1 mapped reads (AIRR Rearrangement)
    * ``<prefix>.assembled.airr.tsv``  -- Stage-3 assembled long-CDR3 reads (if ``--assemble``)
    * ``<prefix>.clones.tsv``          -- corrected clonotype table (folds in the assembled clones)
    * ``<prefix>.arda.json``           -- merged run report
    * ``<prefix>.stats.tsv``           -- run QC, long format (see ``arda stats``)

    Run ``arda map`` / ``assemble`` / ``correct`` / ``shm`` separately to tune their own knobs.
    """
    _mode_run("rnaseq", exact=exact, indel_rescue=False,
              r1=r1, r2=r2, out_dir=out_dir, out_prefix=out_prefix, organism=organism,
              threads=threads, reconstruct=reconstruct, min_score=min_score,
              kmer=(None if kmer == 0 else kmer), assemble=assemble,
              complete_only=complete_only, map_d=map_d, d_max_evalue=d_max_evalue,
              limit=(limit or None), ec_mode=ec_mode, min_junction_q=min_junction_q,
              clonotype_key=clonotype_key, call_level=call_level,
              shm=shm, isotype=isotype,
              complete_junction_nt=complete_junction_nt)


@app.command("amplicon")
def amplicon_mode(
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
    out_prefix: str = typer.Option(
        ..., "--out-prefix", "-p",
        help="Output basename. Writes <prefix>.airr.tsv, <prefix>.clones.tsv, <prefix>.arda.json and <prefix>.stats.tsv."),
    out_dir: Path = typer.Option(Path("."), "--out-dir", "-d", help="Directory for the outputs."),
    organism: str = typer.Option("human", help="Reference organism."),
    threads: int = typer.Option(0, help="mmseqs threads (0 = all cores)."),
    assemble: bool = typer.Option(
        True, "--assemble/--no-assemble",
        help="Stage 3: assemble long-CDR3 contigs no single read spans and fold them into the "
             "clonotype table."),
    shm: str = typer.Option("framework", "--shm", help=_SHM_HELP),
    complete_junction_nt: int = typer.Option(
        0, "--complete-junctions", help=_COMPLETE_JUNCTION_HELP),
    isotype: bool = typer.Option(True, "--isotype/--no-isotype", help=_ISOTYPE_HELP),
    call_level: str = typer.Option("allele", "--call-level", help=_CALL_LEVEL_HELP),
    map_d: bool = typer.Option(True, "--map-d/--no-map-d", help="Map D segments in all stages."),
    d_max_evalue: Optional[float] = typer.Option(None, "--d-max-evalue", help=_D_EVALUE_HELP),
    ec_mode: str = typer.Option(
        "amplicon", "--ec-mode",
        help="Denoising preset; `amplicon` is this mode's default. It adds the quality-directed "
             "rescue searching WIDE (12 subs, 50x abundance ratio), because a real clonotype in a "
             "targeted library is deep, so a 1-read neighbour of an abundant clone is almost "
             "always error. `fast` = the abundance model only; `accurate` = + --min-junction-q 20. "
             "⛔ Nothing in any mode discards a read."),
    min_junction_q: Optional[int] = typer.Option(
        None, "--min-junction-q",
        help="Explicit Phred floor for the discriminating base; overrides --ec-mode's preset."),
    clonotype_key: str = typer.Option(
        "full", "--clonotype-key",
        help="`full` (default) = (locus, v_call, j_call, junction). `junction` = (locus, junction) "
             "with V/J canonicalised to the junction's majority first -- it collapses call splits."),
    complete_only: bool = typer.Option(
        True, "--complete-only/--all-junctions", help="Keep only complete junctions."),
    reconstruct: bool = typer.Option(
        False, "--reconstruct", help="Merge overlapping paired mates into one fragment."),
    min_score: float = typer.Option(
        75.0, "--min-score", help="Min MMseqs2 bit score to keep a mapped read (0 = keep all)."),
    kmer: int = typer.Option(12, "--kmer", "-k", help="MMseqs2 -k (the memory knob)."),
    limit: int = typer.Option(
        0, "--limit", "-n", help="Analyse only the first N reads / pairs. 0 = whole file."),
    indel_rescue: bool = typer.Option(
        False, "--indel-rescue",
        help="Route reads whose seed votes show TWO diagonals to the gapped rescue -- the "
             "signature of an indel, visible before any extension runs. Measured on 341,294 real "
             "IGH mates: 3.18 %% of reads carry a V indel and the rate tracks SHM load (0.74 %% at "
             ">=98 %% V identity, 8.00 %% below 90 %%), because AID makes indels and not only "
             "substitutions. Reads are REROUTED, never dropped. Its value is a per-library call "
             "(+181 reads on a hypermutated repertoire, -14 on a naive one), so it never rides "
             "the preset."),
    exact: bool = typer.Option(False, "--exact", help=_EXACT_HELP),
) -> None:
    """TARGETED RepSeq / 5'RACE amplicon -> clonotypes, with the amplicon speed preset.

    For primer-anchored libraries whose reads span V into J. Enables ``--two-pass
    --fast-segments --v-only-on-segment``: on an IGH RepSeq amplicon at 32 threads that is
    **316.44 s -> 76.25 s (4.15x)** at 4,018 -> 1,479 MB, and on a 100 k-read TRA amplicon it puts
    arda at **5.35 s wall against MiXCR's 5.90**, at 3.6x less CPU and 4.8x less RSS.
    ``--exact`` turns the preset off.

    Same four outputs as ``arda rnaseq``; see it for the list.
    """
    _mode_run("amplicon", exact=exact, indel_rescue=indel_rescue,
              r1=r1, r2=r2, out_dir=out_dir, out_prefix=out_prefix, organism=organism,
              threads=threads, reconstruct=reconstruct, min_score=min_score,
              kmer=(None if kmer == 0 else kmer), assemble=assemble,
              complete_only=complete_only, map_d=map_d, d_max_evalue=d_max_evalue,
              limit=(limit or None), ec_mode=ec_mode, min_junction_q=min_junction_q,
              clonotype_key=clonotype_key, call_level=call_level,
              shm=shm, isotype=isotype,
              complete_junction_nt=complete_junction_nt)


@app.command("singlecell")
def singlecell_mode() -> None:
    """RESERVED — single-cell (10x) mode is not implemented yet.

    The name is taken now so the three-mode surface is stable and a future release adds behaviour
    rather than a new command. arda has **no barcode or UMI concept at all**, which is the real
    single-cell gap -- not the assembler. See ROADMAP.md.
    """
    typer.echo(
        "arda singlecell is not implemented yet (scheduled -- see ROADMAP.md, single-cell).\n"
        "Today, for 10x contigs:\n"
        "  arda annotate -i all_contig.fasta -o contigs.airr.tsv\n"
        "  arda correct  -i contigs.airr.tsv -o contigs.clones.tsv\n"
        "Barcode demultiplexing and UMI consensus are not in arda.", err=True)
    raise typer.Exit(code=2)


@cluster_app.command("split")
def rnaseq_split(
    r1: Path = typer.Option(..., "--r1", help="FASTQ (single-end, or R1 of a pair)."),
    out_dir: Path = typer.Option(..., "--out-dir", "-d", help="Directory for the shard FASTQs."),
    shards: int = typer.Option(..., "--shards", help="Number of contiguous blocks."),
    r2: Optional[Path] = typer.Option(None, "--r2", help="R2 FASTQ for paired input."),
) -> None:
    """Split paired FASTQ into contiguous blocks of read pairs, for a sharded Stage 1.

    Unlike ``arda cluster split-fasta`` (FASTA/amplicon) this keeps the quality strings and never separates a
    fragment's two mates. Blocks are contiguous, so concatenating the per-shard AIRR in shard
    order reproduces the single-node row order exactly.
    """
    from .cluster import split_pairs

    written = split_pairs(r1, out_dir, shards=shards, r2=r2)
    typer.echo(f"[arda] wrote {len(written)} shard(s) to {out_dir}")


@cluster_app.command("reduce")
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
    d_max_evalue: Optional[float] = typer.Option(
        None, "--d-max-evalue", help=_D_EVALUE_HELP),
) -> None:
    """Merge a sharded Stage 1, then run Stages 2-3 ONCE over the whole thing.

    `assemble` and `correct` are global: run per shard, a clone split across shards is counted
    once per shard and contigs that tile across shards are never built. So the sharded path
    distributes only ``map``, and this is the step that finishes the job.
    """
    from .rnaseq import pipeline

    pipeline.reduce(shard_dir, out_dir, out_prefix, organism=organism, threads=threads,
                    assemble=assemble, complete_only=complete_only, map_d=map_d,
                    d_max_evalue=d_max_evalue)
    typer.echo(f"[arda] wrote {out_dir / f'{out_prefix}.clones.tsv'}")


@cluster_app.command("submit")
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
    mem: str = typer.Option("8G", help="Memory per array task (Stage 1 is flat, ~300-650 MB)."),
    reduce_time: str = typer.Option("08:00:00", help="Walltime for the reduce step."),
    reduce_mem: str = typer.Option(
        "16G", help="Memory for reduce. Stage 3 holds the clone set: budget ~4 GB, more for a "
                    "B-cell-rich sample (2,071.7 MB measured on 28,444 clonotypes)."),
    submit: bool = typer.Option(False, "--submit", help="Run the generated script."),
) -> None:
    """Write (and optionally submit) a SLURM script: split → array-``map`` → reduce.

    Only Stage 1 is distributed; Stages 2-3 run once over the merged AIRR, through the same
    code path the mode commands use. With contiguous pair shards the result is
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



@app.command("resolve-ties")
def resolve_ties_cmd(
    input: Path = typer.Option(..., "--input", "-i", help="AIRR TSV from `annotate` or `rnaseq map`."),
    output: Path = typer.Option(..., "--output", "-o", help="AIRR TSV with widened, ranked calls."),
    organism: str = typer.Option("human", "--organism"),
    segments: str = typer.Option("v,j", "--segments",
                                 help="Which calls to widen (comma-joined): v, j."),
    rank: bool = typer.Option(True, "--rank/--no-rank",
                              help="Second pass: put the allele the WHOLE LIBRARY supports first."),
) -> None:
    """Widen ``v_call``/``j_call`` to every germline the read's alignment cannot rule out.

    A read aligned over ``[germline_start, germline_end]`` is explained exactly as well by any
    germline carrying that same stretch, so naming one is a claim the data does not support.
    Measured against IgBLAST on a Ramos library: arda emitted **0 multi-gene tie lists in 504
    calls** where IgBLAST emitted **11.68 %**, and on the 104 reads where they disagreed, aligning
    each read to both germlines showed **59 of 60 fit identically** (identity 1.0000 over 63-70 nt).
    Neither was right; both were overconfident.

    ⛔ This ADDS no alignment. The tie is a string comparison against the reference over the span
    already aligned, so it costs neither the memory nor the time that keeping `top_hit` before
    `convertalis` was protecting (that collapse made the alignment TSV 2.88x smaller AND took
    allele agreement .9735 -> .9956; this does not undo it).

    ⛔ It is a SEPARATE COMMAND, not a flag on `map`, because the ranking needs every read before it
    can order any of them -- two passes over one file. And it changes `v_call`/`j_call` on every
    library: a consumer that splits on `,` and takes `[0]` sees the better answer, one that treats
    the field as a single gene sees a new shape.
    """
    from .annotate.ties import resolve_airr

    segs = tuple(x.strip() for x in segments.split(",") if x.strip())
    rep = resolve_airr(input, output, organism=organism, segments=segs, rank=rank)
    log.info("resolve-ties: %d rows", rep["rows"])
    typer.echo(str(output))


if __name__ == "__main__":
    app()
