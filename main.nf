#!/usr/bin/env nextflow
/*
 * Digenome-seq (DSB) and nDigenome-seq (SSB) cleavage calling from Illumina
 * WGS, or from aligned long reads with --long_reads.
 *
 *   SAMPLESHEET -> PREPARE_INDEX -> FASTP -> ALIGN -> CALL_CHUNK (n per sample)
 *   -> FINALIZE -> MULTIQC
 *
 * Long-read samples are aligned BAMs, so they skip PREPARE_INDEX, FASTP, and
 * ALIGN. Calling rules live in bin/cleavage/ and docs/cleavage_algorithm.md.
 */

include { validateParameters } from 'plugin/nf-schema'

workflow {
    validateParameters()
    def genome = params.long_reads ? null : selectedGenome()
    def settings = callerSettings()
    def settings_json = groovy.json.JsonOutput.toJson(settings)
    writeRunInfo(genome, settings)

    def code = files("${projectDir}/bin/cleavage/*.py")
    def samples = SAMPLESHEET(file(params.input), code).splitJson()
    def aligned = params.long_reads ? LONG_READS(samples) : SHORT_READS(samples, genome)

    // Controls aren't called; each treated sample is called against its named
    // control, or with empty placeholders when it has none.
    def vcfs = samples.map { record ->
        tuple(record.sample, record.variant_vcf ? [file(record.variant_vcf), file(record.variant_index)] : [[], []])
    }
    def bams = aligned.bam.branch { meta, _bam, _bai ->
        control: meta.is_control
        treated: true
    }
    def controls = bams.control.map { meta, bam, bai -> tuple(meta.sample, bam, bai) }
    def treated = bams.treated
        .map { meta, bam, bai -> tuple(meta.sample, meta, bam, bai) }
        .join(vcfs)
        .map { _name, meta, bam, bai, vcf -> tuple(meta.control, meta, bam, bai, vcf) }
    def controlled = treated
        .filter { control, _meta, _bam, _bai, _vcf -> control }
        .combine(controls, by: 0)
        .map { _control, meta, bam, bai, vcf, control_bam, control_bai -> tuple(meta, bam, bai, control_bam, control_bai, vcf[0], vcf[1]) }
    def uncontrolled = treated
        .filter { control, _meta, _bam, _bai, _vcf -> !control }
        .map { _control, meta, bam, bai, vcf -> tuple(meta, bam, bai, [], [], vcf[0], vcf[1]) }
    def chunk_jobs = controlled.mix(uncontrolled).combine(channel.of(0..<params.cleavage_chunks))

    def blacklist = params.genome_blacklist ? file(params.genome_blacklist) : []
    CALL_CHUNK(chunk_jobs, blacklist, code, settings_json)
    FINALIZE(CALL_CHUNK.out.groupTuple(size: params.cleavage_chunks), code, settings_json)

    MULTIQC(aligned.qc.mix(FINALIZE.out.multiqc).collect())
}

// Long-read samples are aligned BAMs, called as given.
workflow LONG_READS {
    take:
    samples

    emit:
    bam = samples.map { record -> tuple(sampleMeta(record), file(record.bam), file(record.bam_index)) }
    qc = channel.empty()
}

// Short reads are trimmed and aligned to --genome; their QC goes to MultiQC.
workflow SHORT_READS {
    take:
    samples
    genome

    main:
    def index_prefix = PREPARE_INDEX(genome.name, genome.fasta, params.ref_cache, file("${projectDir}/bin/prepare_bwamem2_index.sh"))
        .map { output -> output.trim() }
    FASTP(samples.map { record ->
        def fastq_1 = record.fastq_1.collect { path -> file(path) }
        def fastq_2 = record.fastq_2.collect { path -> file(path) }
        tuple(sampleMeta(record) + [single_end: record.single_end], fastq_1, fastq_2)
    })
    ALIGN(FASTP.out.reads, index_prefix)

    emit:
    bam = ALIGN.out.bam
    qc = FASTP.out.qc.mix(ALIGN.out.qc.flatten())
}

// What every process needs to know about a sample.
def sampleMeta(record) {
    return [sample: record.sample, control: record.control, is_control: record.is_control]
}

// The configured genome that --genome names, by its name or an alias in any
// case, as [name, fasta]. Short reads are aligned to it; the index cache and
// the run record use its name, so every alias shares one index.
def selectedGenome() {
    if (!params.genome) {
        error("--genome is required unless --long_reads is set")
    }
    def matches = params.genomes.findAll { name, genome ->
        ([name] + (genome.aliases ?: [])).any { alias -> alias.toString().equalsIgnoreCase(params.genome.toString()) }
    }
    if (matches.size() > 1) {
        error("--genome '${params.genome}' matches more than one configured genome: ${matches.keySet().join(', ')}")
    }
    if (!matches) {
        def configured = params.genomes.collect { name, genome -> genome.aliases ? "${name} (${genome.aliases.join(', ')})" : name }
        error("Unknown --genome '${params.genome}'. Configured genomes: ${configured.join('; ')}")
    }
    def match = matches.entrySet().first()
    return [name: match.key, fasta: file(match.value.fasta).toString()]
}

// Settings for bin/cleavage, named like the parameters. --keep_multimappers
// lets MAPQ-0 primary alignments count, together with `bwa-mem2 -a` (ALIGN)
// and no fastp low-complexity filter (FASTP). --long_reads loosens the clip
// and indel limits, because nanopore reads routinely carry both.
def callerSettings() {
    def multimappers = params.keep_multimappers
    def long_reads = params.long_reads
    return [
        analysis: params.analysis,
        keep_multimappers: multimappers,
        long_reads: long_reads,
        digenome_overhang: params.digenome_overhang,
        digenome_pair_window: params.digenome_pair_window,
        digenome_min_mapq: multimappers ? 0 : params.digenome_min_mapq,
        digenome_forward_cutoff: params.digenome_forward_cutoff,
        digenome_reverse_cutoff: params.digenome_reverse_cutoff,
        digenome_depth_cutoff: params.digenome_depth_cutoff,
        digenome_fraction_cutoff: params.digenome_fraction_cutoff,
        digenome_pair_score_cutoff: params.digenome_pair_score_cutoff,
        ndigenome_min_count: params.ndigenome_min_count,
        ndigenome_min_fraction: params.ndigenome_min_fraction,
        ndigenome_min_mapq: multimappers ? 0 : params.ndigenome_min_mapq,
        ndigenome_opposite_window: params.ndigenome_opposite_window,
        ndigenome_ambiguous_min_count: params.ndigenome_ambiguous_min_count,
        ndigenome_ambiguous_min_fraction: params.ndigenome_ambiguous_min_fraction,
        cleavage_artifact_window: params.cleavage_artifact_window,
        cleavage_max_softclip_fraction: long_reads ? 1.0 : params.cleavage_max_softclip_fraction,
        cleavage_max_indel_fraction: long_reads ? 1.0 : params.cleavage_max_indel_fraction,
        cleavage_min_support_mean_mapq: multimappers ? 0 : params.cleavage_min_support_mean_mapq,
        cleavage_control_min_depth: params.cleavage_control_min_depth,
        cleavage_control_max_fraction: params.cleavage_control_max_fraction,
        cleavage_control_min_fold: params.cleavage_control_min_fold,
        cleavage_control_max_q: params.cleavage_control_max_q,
    ]
}

// Record what this run used before any task starts. Long reads arrive
// aligned, so they use no genome or index cache.
def writeRunInfo(genome, settings) {
    def info_dir = file("${params.outdir}/pipeline_info")
    info_dir.mkdirs()
    def info = [
        analysis: params.analysis,
        genome: genome?.name,
        fasta: genome?.fasta,
        ref_cache: genome ? file(params.ref_cache).toString() : null,
        genome_blacklist: params.genome_blacklist ? file(params.genome_blacklist).toString() : null,
        cleavage_chunks: params.cleavage_chunks,
        caller_settings: settings,
        nextflow_version: nextflow.version.toString(),
        session_id: workflow.sessionId.toString(),
    ]
    file("${info_dir}/analysis_parameters.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info)) + '\n'
    file("${projectDir}/containers/checksums.sha256").copyTo("${info_dir}/container_checksums.sha256")
    file("${projectDir}/containers/sources.tsv").copyTo("${info_dir}/container_sources.tsv")
}

// Single-quote a user-supplied value for bash.
def quote(value) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process SAMPLESHEET {

    input:
    path samplesheet, stageAs: 'input.csv'
    path code, stageAs: 'cleavage/*'

    output:
    path 'samples.json'

    // Stub runs check the sheet too: without a stub block, the script runs.
    script:
    """
    python3 -m cleavage samplesheet input.csv samples.json --analysis ${params.analysis} ${params.long_reads ? '--long-reads' : ''}
    """
}

process PREPARE_INDEX {
    tag "${genome}"
    // Always re-check the shared cache; a cached task could point at a
    // removed index.
    cache false

    input:
    val genome
    val fasta
    val ref_cache
    path helper

    output:
    stdout

    script:
    """
    bash ${helper} --genome ${quote(genome)} --fasta ${quote(fasta)} --cache-dir ${quote(ref_cache)}
    """

    stub:
    """
    echo ${quote("${ref_cache}/stub/${genome}")}
    """
}

process FASTP {
    tag "${meta.sample}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.fastp.*'
    publishDir "${params.outdir}/trimmed_fastqs", mode: 'copy', pattern: '*.trimmed.*', enabled: params.publish_trimmed_fastqs

    input:
    tuple val(meta), path(fastq_1, stageAs: 'in/r1_*.fastq.gz', arity: '1..*'), path(fastq_2, stageAs: 'in/r2_*.fastq.gz', arity: '0..*')

    output:
    tuple val(meta), path("${meta.sample}.trimmed.R*.fastq.gz"), emit: reads
    path "${meta.sample}.fastp.*", emit: qc

    script:
    def r1 = fastq_1.size() == 1 ? fastq_1[0] : 'r1.fastq.gz'
    def r2 = fastq_2.size() == 1 ? fastq_2[0] : 'r2.fastq.gz'
    def paired = meta.single_end ? '' : "--in2 ${r2} --out2 ${meta.sample}.trimmed.R2.fastq.gz --detect_adapter_for_pe"
    def low_complexity = params.analysis == 'digenome' && !params.keep_multimappers && params.fastp_low_complexity_filter
    """
    ${fastq_1.size() > 1 ? "cat ${fastq_1.join(' ')} > r1.fastq.gz" : ''}
    ${fastq_2.size() > 1 ? "cat ${fastq_2.join(' ')} > r2.fastq.gz" : ''}
    fastp \\
        --in1 ${r1} \\
        --out1 ${meta.sample}.trimmed.R1.fastq.gz \\
        ${paired} \\
        --html ${meta.sample}.fastp.html \\
        --json ${meta.sample}.fastp.json \\
        --thread ${task.cpus} \\
        --trim_poly_g \\
        --trim_poly_x \\
        --qualified_quality_phred ${params.fastp_qualified_quality_phred} \\
        --length_required ${params.fastp_length_required} \\
        ${low_complexity ? "--low_complexity_filter --complexity_threshold ${params.fastp_complexity_threshold}" : ''} \\
        ${params.fastp_extra_args}
    """

    stub:
    """
    touch ${meta.sample}.trimmed.R1.fastq.gz ${meta.sample}.fastp.html ${meta.sample}.fastp.json
    ${meta.single_end ? '' : "touch ${meta.sample}.trimmed.R2.fastq.gz"}
    """
}

process ALIGN {
    tag "${meta.sample}"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.txt'

    input:
    tuple val(meta), path(reads)
    val index_prefix

    output:
    tuple val(meta), path("${meta.sample}.sorted.markdup.bam"), path("${meta.sample}.sorted.markdup.bam.bai"), emit: bam
    tuple path("${meta.sample}.flagstat.txt"), path("${meta.sample}.stats.txt"), path("${meta.sample}.markdup.metrics.txt"), emit: qc

    script:
    if (task.cpus < 2) {
        error("ALIGN requires at least 2 CPUs")
    }
    // bwa-mem2 and the piped sort share the CPUs; later steps use all but one.
    def sort_threads = task.cpus >= 4 ? 2 : 0
    def bwa_threads = Math.max(1, task.cpus - sort_threads - 1)
    def threads = task.cpus - 1
    def all_alignments = params.keep_multimappers ? '-a' : ''
    // bwa-mem2 turns each \\t into a tab.
    def read_group = "@RG\\tID:${meta.sample}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}"
    // Paired reads are name-sorted so fixmate can add the mate tags markdup needs.
    def sort_for_markdup = meta.single_end
        ? "samtools sort -@ ${sort_threads} -o positionsort.bam -"
        : """samtools sort -@ ${sort_threads} -n -o namesort.bam -
    samtools fixmate -@ ${threads} -m namesort.bam fixmate.bam
    samtools sort -@ ${threads} -o positionsort.bam fixmate.bam"""
    """
    bwa-mem2 mem ${all_alignments} -t ${bwa_threads} -R '${read_group}' ${quote(index_prefix)} ${reads} \\
        | ${sort_for_markdup}
    samtools markdup -@ ${threads} -s positionsort.bam ${meta.sample}.sorted.markdup.bam \\
        2> ${meta.sample}.markdup.metrics.txt
    samtools index -@ ${threads} ${meta.sample}.sorted.markdup.bam
    samtools flagstat -@ ${threads} ${meta.sample}.sorted.markdup.bam > ${meta.sample}.flagstat.txt
    samtools stats -@ ${threads} ${meta.sample}.sorted.markdup.bam > ${meta.sample}.stats.txt
    """

    stub:
    """
    touch ${meta.sample}.sorted.markdup.bam ${meta.sample}.sorted.markdup.bam.bai
    touch ${meta.sample}.flagstat.txt ${meta.sample}.stats.txt ${meta.sample}.markdup.metrics.txt
    """
}

process CALL_CHUNK {
    tag "${meta.sample}:${chunk}"

    // Long-read BAMs keep the user's file names, so the treated and control
    // BAMs are staged apart, each beside its index.
    input:
    tuple val(meta),
        path(bam, stageAs: 'treated/*'), path(bai, stageAs: 'treated/*'),
        path(control_bam, stageAs: 'control/*'), path(control_bai, stageAs: 'control/*'),
        path(vcf), path(vcf_index), val(chunk)
    path blacklist
    path code, stageAs: 'cleavage/*'
    val settings

    output:
    tuple val(meta), path("chunk_*.json"), path("chunk_*.jsonl.gz")

    script:
    def prefix = "chunk_${chunk.toString().padLeft(3, '0')}"
    def control = control_bam ? "--control-bam ${quote(control_bam)} --control-sample ${meta.control}" : ''
    """
    cat > settings.json <<'JSON'
    ${settings}
    JSON
    python3 -m cleavage call \\
        --settings settings.json \\
        --bam ${quote(bam)} \\
        --sample ${meta.sample} \\
        --chunk ${chunk} \\
        --chunks ${params.cleavage_chunks} \\
        --out-prefix ${prefix} \\
        ${control} \\
        ${vcf ? "--vcf ${quote(vcf)}" : ''} \\
        ${blacklist ? "--blacklist ${quote(blacklist)}" : ''}
    """

    stub:
    def prefix = "chunk_${chunk.toString().padLeft(3, '0')}"
    """
    touch ${prefix}.json ${prefix}.jsonl.gz
    """
}

process FINALIZE {
    tag "${meta.sample}"
    publishDir "${params.outdir}/${params.analysis}", mode: 'copy', pattern: "*.${params.analysis}*"
    publishDir "${params.outdir}/pipeline_info/cleavage_chunks", mode: 'copy', pattern: '*.cleavage_chunks.tsv'

    input:
    tuple val(meta), path(summaries), path(records)
    path code, stageAs: 'cleavage/*'
    val settings

    output:
    path "${meta.sample}.${params.analysis}.*"
    path "${meta.sample}.cleavage_chunks.tsv"
    path "${meta.sample}.${params.analysis}_mqc.tsv", emit: multiqc

    script:
    """
    cat > settings.json <<'JSON'
    ${settings}
    JSON
    python3 -m cleavage finalize \\
        --settings settings.json \\
        --sample ${meta.sample} \\
        ${meta.control ? "--control-sample ${meta.control}" : ''} \\
        --out-prefix ${meta.sample} \\
        ${summaries}
    """

    stub:
    def stem = "${meta.sample}.${params.analysis}"
    """
    touch ${stem}.all.tsv ${stem}.high_confidence.tsv ${stem}.manual_review.tsv ${stem}.artifact.tsv
    touch ${stem}.bed ${stem}.qc.json ${stem}_mqc.tsv ${meta.sample}.cleavage_chunks.tsv
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    """
    multiqc .
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    """
}
