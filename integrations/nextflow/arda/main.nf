// arda RNA-seq -> AIRR clonotypes, as a drop-in nf-core-style local module.
// One call to `arda rnaseq run` (map + assemble + correct) per sample; publishes to ${params.outdir}/arda/.
// See ./README.md for how to wire this into an nf-core/rnaseq (or similar) pipeline.

process ARDA {
    tag "$meta.id"
    // arda is CPU-bound: the MMseqs2 search dominates. ~40-50k reads/s on 32 cores; a full-depth
    // ~100M-read sample takes ~45 min.
    //
    // MEMORY. Stage 1 (`map`) is FLAT -- 300-650 MB at any read depth, because it streams. What
    // scales is Stage 3 (`correct`), which holds the whole clone set in memory: it peaked
    // 2,071.7 MB on a B-cell-rich tumour with 28,444 clonotypes from 105 M reads, while a COLDER
    // 139 M-read sample -- more reads, almost no repertoire -- peaked 549 MB. So budget ~4 GB per
    // task and size it by expected repertoire richness, not by FASTQ size. `--no-assemble` keeps
    // the run on the flat mapping-only profile if you must cap tightly.
    label 'process_high'

    // arda is pip-installable (PyPI: arda-mapper) and needs the mmseqs2 binary.
    //   -profile conda    -> works out of the box from environment.yml
    //   -profile docker   -> build the image from the Dockerfile beside this module and push it to
    //                        your registry, then point `container` at it (or override in a config).
    //
    // ⛔ `conda` is declared HERE, as a directive inside the process body -- never as
    // `process { conda = ... }` at config scope. At process scope Nextflow applies it to EVERY
    // process in the pipeline, which silently builds a different arda *and* a different aligner
    // than the one this module was validated with. If you must override it from a config, it goes
    // under `withName: 'ARDA' { conda = ... }` -- see nextflow.config beside this file.
    conda "${moduleDir}/environment.yml"
    container "arda-mapper:2.12.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.clones.tsv"),                            emit: clones
    // Not `*.airr.tsv`: that also matches the assembler's `*.assembled.airr.tsv`, and the two
    // would arrive on `airr` as a 2-element list. Name the mapped AIRR exactly.
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.airr.tsv"),  emit: airr
    tuple val(meta), path("*.assembled.airr.tsv"),                    emit: assembled_airr, optional: true
    tuple val(meta), path("*.arda.json"),                             emit: report
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2     = meta.single_end ? '' : "--r2 ${reads[1]}"

    // ── REGIME PRESET ──────────────────────────────────────────────────────────────────────────
    // arda has two tuning paths and they do NOT compose. Getting the choice backwards is not an
    // error -- it is a silent 2-4x slowdown -- so the combination is selected here, once, BY NAME:
    //
    //   'amplicon'  --two-pass --fast-segments --v-only-on-segment   targeted RepSeq / 5'RACE
    //   'bulk'      --prefilter                                      whole-transcriptome RNA-seq
    //   'default'   (nothing)                                        the shipped one-pass path
    //
    // ⛔ `--two-pass` ALONE is a LOSS -- 0.762x on bulk and 0.87x on an IGH amplicon -- which is
    // why no preset here ever emits it without `--fast-segments`.
    //
    // Set it globally with `--regime amplicon`, or per sample with a `regime` key in the meta map
    // (a sheet may legitimately mix the two); the meta map wins.
    def presets = [
        'bulk'    : '--prefilter',
        'amplicon': '--two-pass --fast-segments --v-only-on-segment',
        'default' : '',
    ]
    // ⛔ `params.getOrDefault(...)` throughout, never a bare `params.x`. This module must stay
    // correct when it is included WITHOUT its nextflow.config, and Nextflow scans the script
    // STATICALLY for `params.<name>` tokens -- so a `containsKey` guard around one still emits
    // "Access to undefined parameter" (a WARN normally, a hard failure under strict mode).
    // Measured: the guarded form warned, this one does not.
    def regime = meta.regime ?: params.getOrDefault('regime', 'bulk') ?: 'bulk'
    if (!presets.containsKey(regime))
        throw new IllegalArgumentException(
            "ARDA: regime must be one of ${presets.keySet()}, got '${regime}' (sample '${meta.id}')")
    def tuning = presets[regime]

    // `--indel-rescue` needs `--fast-segments`, i.e. the amplicon preset; arda would otherwise
    // ignore it silently. Its value tracks SHM load (+181 reads on a hypermutated repertoire,
    // -14 on a naive one), so it stays a deliberate per-library call and never rides a preset.
    if (params.getOrDefault('arda_indel_rescue', false)) {
        if (regime != 'amplicon')
            throw new IllegalArgumentException(
                "ARDA: --arda_indel_rescue needs regime 'amplicon' (it requires --fast-segments); " +
                "got regime '${regime}' for sample '${meta.id}'")
        tuning += ' --indel-rescue'
    }

    // ⛔ Pin the aligner to the one THIS task's environment provides. An mmseqs index is only
    // reusable by the release that built it, and a cluster's cache marker differs from the shipped
    // one -- unpinned, arda may reject the precompiled reference index and rebuild a private cache
    // per task, or auto-fetch a third build, with no error and results that are not comparable.
    // $ARDA_MMSEQS is arda's highest-precedence, never-second-guessed override. An outer value
    // still wins; an empty one falls through to arda's normal resolution.
    def mmseqs_path = params.getOrDefault('arda_mmseqs', null)
    def mmseqs_pin  = mmseqs_path
        ? "export ARDA_MMSEQS='${mmseqs_path}'"
        : 'export ARDA_MMSEQS="${ARDA_MMSEQS:-$(command -v mmseqs)}"'
    """
    ${mmseqs_pin}

    arda rnaseq run \\
        --r1 ${reads[0]} ${r2} \\
        --out-prefix ${prefix} \\
        --out-dir . \\
        --threads ${task.cpus} \\
        ${tuning} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arda: \$(arda --version)
        mmseqs2: \$(mmseqs version 2>/dev/null || echo unknown)
        regime: ${regime}
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
