// arda RNA-seq -> AIRR clonotypes, as a drop-in nf-core-style local module.
// One call to `arda <mode>` (map + assemble + correct) per sample; publishes to ${params.outdir}/arda/.
// See ./README.md for how to wire this into an nf-core/rnaseq (or similar) pipeline.
//
// ⛔ Requires arda >= 2.16.0. `arda rnaseq run` was REMOVED there: the regime is now the command
// name (`arda rnaseq` / `arda amplicon`), and each mode carries its own speed configuration. This
// module used to build that flag string itself, which is exactly the workaround the CLI absorbed.

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
    container "arda-mapper:2.20.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.clones.tsv"),                            emit: clones
    // Not `*.airr.tsv`: that also matches the assembler's `*.assembled.airr.tsv`, and the two
    // would arrive on `airr` as a 2-element list. Name the mapped AIRR exactly.
    tuple val(meta), path("${task.ext.prefix ?: meta.id}.airr.tsv"),  emit: airr
    tuple val(meta), path("*.assembled.airr.tsv"),                    emit: assembled_airr, optional: true
    tuple val(meta), path("*.arda.json"),                             emit: report
    // Run QC (arda >= 2.20.0), written by every mode run. `optional` so this module still runs
    // against an older arda in an environment someone pinned themselves.
    tuple val(meta), path("*.stats.tsv"),                             emit: stats, optional: true
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def r2     = meta.single_end ? '' : "--r2 ${reads[1]}"

    // ── REGIME = COMMAND NAME ──────────────────────────────────────────────────────────────────
    // arda has two tuning paths and they do NOT compose. Getting the choice backwards is not an
    // error -- it is a silent 2-4x slowdown -- so the regime is chosen here, once, BY NAME, and
    // arda itself owns which flags that implies:
    //
    //   'amplicon'  arda amplicon            targeted RepSeq / 5'RACE
    //   'bulk'      arda rnaseq              whole-transcriptome RNA-seq
    //   'default'   arda rnaseq --exact      the shipped one-pass path, no speedups
    //
    // ⛔ 'default' is NOT `arda rnaseq`: that mode enables `--prefilter`, which costs ~0.15 % of
    // mapped reads (122 bulk datasets, up to 2.46 % on one library). `--exact` is what reproduces
    // the pre-2.16.0 default output.
    //
    // Set it globally with `--regime amplicon`, or per sample with a `regime` key in the meta map
    // (a sheet may legitimately mix the two); the meta map wins.
    def modes = [
        'bulk'    : ['cmd': 'rnaseq',   'extra': ''],
        'amplicon': ['cmd': 'amplicon', 'extra': ''],
        'default' : ['cmd': 'rnaseq',   'extra': '--exact'],
    ]
    // ⛔ `params.getOrDefault(...)` throughout, never a bare `params.x`. This module must stay
    // correct when it is included WITHOUT its nextflow.config, and Nextflow scans the script
    // STATICALLY for `params.<name>` tokens -- so a `containsKey` guard around one still emits
    // "Access to undefined parameter" (a WARN normally, a hard failure under strict mode).
    // Measured: the guarded form warned, this one does not.
    def regime = meta.regime ?: params.getOrDefault('regime', 'bulk') ?: 'bulk'
    if (!modes.containsKey(regime))
        throw new IllegalArgumentException(
            "ARDA: regime must be one of ${modes.keySet()}, got '${regime}' (sample '${meta.id}')")
    def mode   = modes[regime].cmd
    def tuning = modes[regime].extra

    // `--indel-rescue` needs the fast segment pass, which only `arda amplicon` enables; arda now
    // REFUSES the flag elsewhere rather than ignoring it, but failing here names the sample. Its
    // value tracks SHM load (+181 reads on a hypermutated repertoire, -14 on a naive one), so it
    // stays a deliberate per-library call and never rides a preset.
    if (params.getOrDefault('arda_indel_rescue', false)) {
        if (regime != 'amplicon')
            throw new IllegalArgumentException(
                "ARDA: --arda_indel_rescue needs regime 'amplicon' (it requires the fast segment " +
                "pass); got regime '${regime}' for sample '${meta.id}'")
        tuning += ' --indel-rescue'
    }

    // SHM SCOPE. `framework` (arda's default) keeps v_identity / v_mutations / j_mutations OUT of
    // the junction; `both` also emits the pre-2.16.0 junction-inclusive values as *_full columns,
    // which is what a pipeline needs if it must reproduce previously published SHM numbers; `off`
    // emits no SHM fields.
    def shm = meta.shm ?: params.getOrDefault('arda_shm', 'framework') ?: 'framework'
    if (!(shm in ['framework', 'both', 'off']))
        throw new IllegalArgumentException(
            "ARDA: arda_shm must be framework|both|off, got '${shm}' (sample '${meta.id}')")
    tuning += " --shm ${shm}"

    // CALL LEVEL. `gene` drops the allele suffix before the clonotype key is formed, collapsing
    // allele-level call splits and allele-only tie lists. ⚠ IGH carries 4.33 alleles/gene, ~2x
    // every other locus, so its effect there is far larger than on TR -- a decision, not a knob.
    def call_level = params.getOrDefault('arda_call_level', 'allele') ?: 'allele'
    if (!(call_level in ['allele', 'gene']))
        throw new IllegalArgumentException(
            "ARDA: arda_call_level must be 'allele' or 'gene', got '${call_level}'")
    if (call_level != 'allele') tuning += " --call-level ${call_level}"

    // DENOISING. `--ec-mode` picks how Stage 2 decides what is an error. ⛔ Since 2.16.0 EACH MODE
    // HAS ITS OWN DEFAULT -- `arda rnaseq` defaults to `rnaseq`, `arda amplicon` to `amplicon` --
    // so leaving this unset no longer means `fast`. Set it explicitly to `fast` for arda's
    // historical behaviour. `accurate`/`amplicon`/`rnaseq` all need Stage 1 to carry per-read
    // junction quality, and the mode turns that on ITSELF -- do not pass `--junction-quality` here
    // as well, and never pass a mode without meaning it: a mode that is accepted and silently does
    // nothing is the failure this project keeps hitting.
    def ec_mode = meta.ec_mode ?: params.getOrDefault('arda_ec_mode', null)
    def ec_modes = ['fast', 'accurate', 'amplicon', 'rnaseq'] as Set
    if (ec_mode && !ec_modes.contains(ec_mode))
        throw new IllegalArgumentException(
            "ARDA: arda_ec_mode must be one of ${ec_modes}, got '${ec_mode}' (sample '${meta.id}')")
    if (ec_mode) tuning += " --ec-mode ${ec_mode}"

    // ⚠ `amplicon` and `rnaseq` differ because their clonotype-SIZE distributions differ, not by
    // taste: an amplicon clonotype is deep, so a 1-read neighbour of an abundant clone is almost
    // always error and the quality rescue can search wide; bulk RNA-seq is 0.02-3 % receptor and
    // its singletons are mostly real, so it stays narrow. Choosing the mode that does not match
    // the library is a real cost, which is why this warns rather than silently accepting it.
    // (Only reachable when it is set EXPLICITLY -- each mode's own default already matches.)
    if ((ec_mode == 'amplicon' && regime == 'bulk') || (ec_mode == 'rnaseq' && regime == 'amplicon'))
        log.warn "ARDA: ec_mode '${ec_mode}' with regime '${regime}' for sample '${meta.id}' -- " +
                 "these presets are tuned to opposite clonotype-size distributions."

    // The clonotype KEY. `junction` canonicalises V/J to the junction's majority so CALL SPLITS
    // collapse -- a junction byte-identical to an abundant clone's under a different V or J call,
    // which no error model can see because there is no discriminating base. Measured cost on a
    // polyclonal TRA amplicon: 0.66 % of clonotypes merge. Off by default: it changes the
    // clonotype identity, so it is a decision, not a tuning knob.
    def clonotype_key = params.getOrDefault('arda_clonotype_key', 'full') ?: 'full'
    if (!(clonotype_key in ['full', 'junction']))
        throw new IllegalArgumentException(
            "ARDA: arda_clonotype_key must be 'full' or 'junction', got '${clonotype_key}'")
    if (clonotype_key != 'full') tuning += " --clonotype-key ${clonotype_key}"

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

    arda ${mode} \\
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
        mode: ${mode}
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
