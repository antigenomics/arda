// arda RNA-seq -> AIRR clonotypes, as a drop-in nf-core-style local module.
// One call to `arda rnaseq run` (map + assemble + correct) per sample; publishes to ${params.outdir}/arda/.
// See ./README.md for how to wire this into an nf-core/rnaseq (or similar) pipeline.

process ARDA {
    tag "$meta.id"
    // arda is CPU-bound (the MMseqs2 search dominates) but very low-memory (<400 MB, independent of
    // depth) -- give it cores, not RAM. ~40k reads/s on 32 cores; a full-depth ~100M-read sample ~45 min.
    label 'process_high'

    // arda is pip-installable (PyPI: arda-mapper) and needs the mmseqs2 binary.
    //   -profile conda    -> works out of the box from environment.yml
    //   -profile docker   -> build the image from the Dockerfile beside this module and push it to
    //                        your registry, then point `container` at it (or override in a config).
    conda "${moduleDir}/environment.yml"
    container "arda-mapper:2.5.0"

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2 = meta.single_end ? '' : "--r2 ${reads[1]}"
    """
    arda rnaseq run \\
        --r1 ${reads[0]} ${r2} \\
        --out-prefix ${prefix} \\
        --out-dir . \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arda: \$(arda --version)
        mmseqs2: \$(mmseqs version 2>/dev/null || echo unknown)
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
