// arda as a V(D)J assignment + clonotype step for nf-core/airrflow.
//
// ⛔ BREAKING vs the pre-2.28 module. That one invented its own `params.regime` and its own meta
// map. This one speaks airrflow: the library protocol comes from `params.library_generation_method`
// and the organism from `meta.species`, both of which an airrflow samplesheet already carries. A
// pipeline that set `--regime` must migrate -- see README.md ("Migrating from the regime module").
// The reason is not tidiness: a second, arda-only name for the protocol is a second place to get
// it wrong, and getting it wrong is a silent 2-4x slowdown rather than an error.
//
// Drop-in for CHANGEO_ASSIGNGENES + CHANGEO_MAKEDB: same `[meta, reads]` in, an AIRR Rearrangement
// TSV out. arda does the IgBLAST work ONCE, offline, when its reference is built, so there is no
// per-run germline database to stage and `--reference_igblast` is not consulted.

process ARDA_ASSIGN {
    tag "$meta.id"

    // arda is CPU-bound: the MMseqs2 search dominates. ~40-50k reads/s on 32 cores; a full-depth
    // ~100M-read sample takes ~45 min.
    //
    // MEMORY. Stage 1 (`map`) is FLAT -- 300-650 MB at any read depth, because it streams. What
    // scales is Stage 3 (`correct`), which holds the whole clone set in memory: it peaked
    // 2,071.7 MB on a B-cell-rich tumour with 28,444 clonotypes from 105 M reads, while a COLDER
    // 139 M-read sample -- more reads, almost no repertoire -- peaked 549 MB. So budget ~4 GB per
    // task and size it by expected repertoire RICHNESS, not by FASTQ size.
    label 'process_high'

    // Never: `conda` is declared HERE, as a directive inside the process body -- never as
    // `process { conda = ... }` at config scope. At process scope Nextflow applies it to EVERY
    // process in the pipeline, which silently builds a different arda *and* a different aligner
    // than this module was validated with. Override it under `withName: 'ARDA_ASSIGN' { ... }`.
    conda "${moduleDir}/environment.yml"
    container "arda-mapper:2.29.0"

    input:
    tuple val(meta), path(reads)

    output:
    // The AIRR Rearrangement TSV, which is what the rest of airrflow consumes. Not `*.airr.tsv`:
    // that also matches the assembler's `*.assembled.airr.tsv` and the two would arrive as a
    // 2-element list. Name it exactly.
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.airr.tsv"), emit: airr
    tuple val(meta), path("*.clones.tsv"),                           emit: clones
    tuple val(meta), path("*.assembled.airr.tsv"),                   emit: assembled_airr, optional: true
    tuple val(meta), path("*.arda.json"),                            emit: report
    tuple val(meta), path("*.stats.tsv"),                            emit: stats, optional: true
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // ── A SAMPLE MAY ARRIVE IN SEVERAL FILES ───────────────────────────────────────────────────
    // airrflow's samplesheet merges repeated `sample_id` rows, exactly as nf-core's does, so one
    // repertoire can reach us as several lanes. It is still ONE repertoire and must give ONE
    // clonotype table. arda maps each read group and concatenates before Stages 2-3, which is
    // byte-identical to the same reads in one file -- so pass them all rather than `cat`-ing.
    //
    // `reads` is the usual nf-core flat list, PAIR-ADJACENT: [r1, r2, r1, r2, ...] paired,
    // [r1, r1, ...] single. `collate` recovers the pairs; a plain 2-element list collates to
    // exactly one pair, so an ordinary one-lane sample is unchanged.
    //
    // Never: declared order, never sorted by filename. `A_L010` sorts before `A_L002`, and the
    // clonotype fold is not permutation-invariant -- `correct` collapses an error child onto the
    // parent it meets first. Let the channel's order stand.
    def groups = meta.single_end ? reads.collect { [it] } : reads.collate(2)
    def inputs = groups.collect { g ->
        meta.single_end ? "--r1 ${g[0]}" : "--r1 ${g[0]} --r2 ${g[1]}"
    }.join(' \\\n        ')
    def name_arg = groups.size() > 1 ? groups.collect { "--id ${prefix}" }.join(' ')
                                     : "--out-prefix ${prefix}"

    // ── PROTOCOL -> arda MODE ──────────────────────────────────────────────────────────────────
    // arda has two tuning paths and they do NOT compose; getting the choice backwards is a silent
    // 2-4x slowdown, not an error. So the choice is DERIVED from the protocol airrflow already
    // knows, and arda itself owns which flags that implies.
    //
    // Never: `params.getOrDefault(...)` throughout, never a bare `params.x`. Nextflow scans the
    // script STATICALLY for `params.<name>` tokens, so even a `containsKey`-guarded bare read
    // emits "Access to undefined parameter" -- a WARN normally, a HARD FAILURE under strict mode.
    // This module must stay correct when it is included WITHOUT its nextflow.config.
    def protocols = [
        // Targeted amplicon: the read already spans V into J.
        'specific_pcr'    : 'amplicon',
        'specific_pcr_umi': 'amplicon',
        'dt_5p_race'      : 'amplicon',
        'dt_5p_race_umi'  : 'amplicon',
        // Whole-transcriptome RNA-seq. airrflow spells this one after the tool it used to reach
        // for; here it selects arda's own bulk path, which is what the name describes.
        'trust4'          : 'rnaseq',
    ]
    def method = params.getOrDefault('library_generation_method', null)
    if (!method)
        throw new IllegalArgumentException(
            "ARDA_ASSIGN: --library_generation_method is required (sample '${meta.id}'). " +
            "This module covers ${protocols.keySet()}.")

    // ⛔ SINGLE CELL IS NOT THIS MODULE, and refusing is the honest answer.
    // `arda cells` exists, but its input is ONE PER-MOLECULE UMI CONSENSUS FASTQ with the cell
    // barcode in the record NAME -- what `migec assemble` or Cell Ranger writes -- not a raw 10x
    // read pair. arda does no barcode demultiplexing and no UMI collapse; both belong upstream.
    // Mapping `sc_10x_genomics` here would hand `arda cells` reads it cannot interpret, and the
    // failure would look like a bad repertoire rather than a wiring error. Run the single-cell
    // path separately -- `arda cells` after the UMI consensus step; see docs/singlecell.rst.
    if (method == 'sc_10x_genomics' || meta.single_cell)
        throw new IllegalArgumentException(
            "ARDA_ASSIGN: sample '${meta.id}' is single-cell (library_generation_method " +
            "'${method}', single_cell=${meta.single_cell}). This module covers BULK protocols " +
            "only: ${protocols.keySet()}. arda's single-cell entry point is `arda cells`, whose " +
            "input is a per-molecule UMI consensus FASTQ with the barcode in the record name -- " +
            "run it downstream of the UMI collapse rather than on raw reads.")

    if (!protocols.containsKey(method))
        throw new IllegalArgumentException(
            "ARDA_ASSIGN: library_generation_method '${method}' is not one of " +
            "${protocols.keySet()} (sample '${meta.id}').")
    def mode = protocols[method]

    // ── ORGANISM ───────────────────────────────────────────────────────────────────────────────
    // From the samplesheet's required `species` column. arda's reference is IMGT-derived, so only
    // the species matters -- there is no assembly to match and no iGenomes lookup to do.
    def species = (meta.species ?: params.getOrDefault('arda_organism', null) ?: 'human')
                  .toString().toLowerCase()

    // ── `params.productive_only` IS DELIBERATELY NOT CONSUMED HERE ─────────────────────────────
    // ⛔ arda HAS no productive filter and must not pretend to. It emits the AIRR `productive`
    // column and leaves the decision to the consumer -- the same stance as its QC surface, which
    // flags and never filters. airrflow's own `productive_only` step runs downstream of this
    // module and applies unchanged.
    //
    // Never: a parameter that is accepted and silently does nothing is the failure mode this
    // project keeps hitting, so it is named here rather than wired to a flag that does not exist.
    // (`arda cells --require-productive` is a DIFFERENT question -- whether a believed SECOND
    // chain must be productive -- and rides that mode's own default.)

    // Never: Pin the aligner to the one THIS task's environment provides. An mmseqs index is only
    // reusable by the release that built it, and a cluster's cache marker differs from the shipped
    // one -- unpinned, arda may reject the precompiled reference index and rebuild a private cache
    // per task, or auto-fetch a third build, with no error and results that are not comparable.
    def mmseqs_path = params.getOrDefault('arda_mmseqs', null)
    def mmseqs_pin  = mmseqs_path
        ? "export ARDA_MMSEQS='${mmseqs_path}'"
        : 'export ARDA_MMSEQS="${ARDA_MMSEQS:-$(command -v mmseqs)}"'
    """
    ${mmseqs_pin}

    arda ${mode} \\
        ${inputs} \\
        ${name_arg} \\
        --out-dir . \\
        --organism ${species} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arda: \$(arda --version)
        mmseqs2: \$(mmseqs version 2>/dev/null || echo unknown)
        library_generation_method: ${method}
        mode: ${mode}
        organism: ${species}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clones.tsv ${prefix}.airr.tsv ${prefix}.assembled.airr.tsv ${prefix}.arda.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arda: \$(arda --version)
    END_VERSIONS
    """
}
