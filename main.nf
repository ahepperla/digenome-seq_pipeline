#!/usr/bin/env nextflow

import groovy.json.JsonOutput

nextflow.enable.dsl = 2

def shellQuote = { value ->
    "'${value.toString().replace("'", "'\"'\"'")}'"
}

def genomeLabel = { String genomeName, genomeConfig ->
    def aliases = genomeConfig.aliases ?: []
    aliases ? "${genomeName} (aliases: ${aliases.join(', ')})" : genomeName
}

def resolveGenome = { String requestedGenome ->
    def requestedLower = requestedGenome.toLowerCase()
    def matches = []
    params.genomes.each { genomeName, genomeConfig ->
        def names = [genomeName as String] +
            ((genomeConfig.aliases ?: []) as List).collect { it as String }
        if (names.collect { it.toLowerCase() }.contains(requestedLower)) {
            matches << (genomeName as String)
        }
    }
    matches = matches.unique()
    if (matches.size() > 1) {
        error "Genome name/alias '${requestedGenome}' is ambiguous and matches: " +
            "${matches.join(', ')}. Fix duplicate aliases in nextflow.config."
    }
    matches ? matches[0] : null
}

def requireContainer = { String name ->
    def configured = params.containers[name]
    if (!configured) {
        error "Missing required container configuration: params.containers.${name}"
    }
    def text = configured as String
    if (text.startsWith('/') && !file(text).exists()) {
        error "Configured ${name} container does not exist: ${text}"
    }
}

if (!params.input) {
    error "Missing required parameter: --input samplesheet.csv"
}
if (!params.genome) {
    error "Missing required parameter: --genome"
}
if (!params.ref_cache) {
    error "Missing required config value: params.ref_cache"
}
String genome_blacklist_path = params.genome_blacklist ?
    params.genome_blacklist.toString().trim() : ''
boolean has_genome_blacklist = genome_blacklist_path != ''
if (
    has_genome_blacklist
    && !file(genome_blacklist_path).exists()
) {
    error "Genome blacklist does not exist: ${genome_blacklist_path}"
}

String selected_analysis = (params.analysis ?: 'digenome').toString().toLowerCase()
if (!(selected_analysis in ['digenome', 'ndigenome'])) {
    error "Unknown --analysis '${params.analysis}'. Expected: digenome or ndigenome"
}

boolean keep_multimappers = params.keep_multimappers as Boolean
int cleavage_chunks = params.cleavage_chunks as int
if (cleavage_chunks < 1) {
    error "--cleavage_chunks must be at least 1"
}
// Cover a pair plus one-step neighboring endpoints that can compete for it.
int digenome_chunk_context =
    Math.abs(params.digenome_overhang as int) +
    3 * (params.digenome_pair_window as int)
int cleavage_chunk_padding = [
    params.cleavage_artifact_window as int,
    params.ndigenome_opposite_window as int,
    digenome_chunk_context
].max() as int
int effective_digenome_min_mapq = keep_multimappers ?
    0 : params.digenome_min_mapq as int
int effective_ndigenome_min_mapq = keep_multimappers ?
    0 : params.ndigenome_min_mapq as int
double effective_min_support_mean_mapq = keep_multimappers ?
    0.0 : params.cleavage_min_support_mean_mapq as double

if (keep_multimappers) {
    log.info(
        "--keep_multimappers enabled: using MAPQ 0 primary alignments and " +
        "disabling the support mean-MAPQ filter; secondary alignments remain diagnostic."
    )
}

['python', 'fastp', 'align', 'cleavage', 'multiqc'].each {
    requireContainer(it)
}

String requested_genome = params.genome as String
String selected_genome = resolveGenome(requested_genome)
if (!selected_genome) {
    def available = params.genomes.collect { genomeName, genomeConfig ->
        genomeLabel(genomeName as String, genomeConfig)
    }.join('; ')
    error "Unknown --genome '${requested_genome}'. Available genomes/aliases: ${available}"
}

String selected_fasta = params.genomes[selected_genome].fasta as String
if (!selected_fasta || selected_fasta == 'null') {
    error "Genome '${selected_genome}' does not have a configured FASTA"
}
if (!file(selected_fasta).exists()) {
    error "Configured FASTA for '${selected_genome}' is not visible: ${selected_fasta}"
}

def bundled_checksums = file("${baseDir}/containers/checksums.sha256").text
    .readLines()
    .findAll { it.trim() && !it.trim().startsWith('#') }
    .collectEntries { line ->
        def fields = line.trim().split(/\s+/, 2)
        [(fields[1]): fields[0]]
    }
def source_manifest_lines = file("${baseDir}/containers/sources.tsv").text
    .readLines()
    .findAll { it.trim() }
def source_manifest_header = source_manifest_lines.first().split('\t', -1)
def container_sources = source_manifest_lines
    .drop(1)
    .collectEntries { line ->
        def fields = line.split('\t', -1)
        def provenance = [:]
        source_manifest_header.drop(1).eachWithIndex { field_name, index ->
            provenance[field_name] = fields[index + 1]
        }
        [(fields[0]): provenance]
    }

def run_info = [
    schema_version: 8,
    analysis: selected_analysis,
    requested_genome: requested_genome,
    resolved_genome: selected_genome,
    fasta: selected_fasta,
    ref_cache: params.ref_cache as String,
    keep_multimappers: keep_multimappers,
    genome_blacklist: genome_blacklist_path,
    multimapper_counting: [
        primary_alignments: 'counted_once_at_the_bwa_selected_primary_placement',
        secondary_and_supplementary_alignments: 'diagnostic_only'
    ],
    automatic_overrides: keep_multimappers ? [
        'digenome_min_mapq=0',
        'ndigenome_min_mapq=0',
        'cleavage_min_support_mean_mapq=0'
    ] : [],
    cleavage: [
        chunks: cleavage_chunks,
        cpus_per_chunk: 1,
        interval_padding: cleavage_chunk_padding,
        artifact_window: params.cleavage_artifact_window,
        max_softclip_fraction: params.cleavage_max_softclip_fraction,
        max_indel_fraction: params.cleavage_max_indel_fraction,
        min_support_mean_mapq: effective_min_support_mean_mapq,
        control_min_depth: params.cleavage_control_min_depth,
        control_max_fraction: params.cleavage_control_max_fraction,
        control_min_fold: params.cleavage_control_min_fold,
        control_max_q: params.cleavage_control_max_q
    ],
    digenome: [
        overhang: params.digenome_overhang,
        pair_window: params.digenome_pair_window,
        min_mapq: effective_digenome_min_mapq,
        forward_cutoff: params.digenome_forward_cutoff,
        reverse_cutoff: params.digenome_reverse_cutoff,
        depth_cutoff: params.digenome_depth_cutoff,
        fraction_cutoff: params.digenome_fraction_cutoff,
        score_cutoff: params.digenome_score_cutoff
    ],
    ndigenome: [
        min_count: params.ndigenome_min_count,
        min_fraction: params.ndigenome_min_fraction,
        min_mapq: effective_ndigenome_min_mapq,
        opposite_window: params.ndigenome_opposite_window,
        ambiguous_min_count: params.ndigenome_ambiguous_min_count,
        ambiguous_min_fraction: params.ndigenome_ambiguous_min_fraction
    ],
    containers: params.containers,
    bundled_container_checksums: bundled_checksums,
    container_sources: container_sources
]
String run_info_json = JsonOutput.prettyPrint(JsonOutput.toJson(run_info))

process RUN_INFO {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    val run_json

    output:
    path "analysis_parameters.json"

    script:
    """
    cat > analysis_parameters.json <<'JSON'
${run_json}
JSON
    """
}

process VALIDATE_SAMPLESHEET {
    tag "samplesheet"

    input:
    path samplesheet
    path validator
    val analysis

    output:
    path "samplesheet.valid.csv", emit: csv

    script:
    """
    python3 ${shellQuote(validator)} ${shellQuote(samplesheet)} samplesheet.valid.csv ${shellQuote(analysis)}
    """
}

process PREPARE_BWAMEM2_INDEX {
    tag "${genome}"
    cache false

    input:
    tuple val(genome), val(fasta), val(ref_cache)
    path index_helper

    output:
    tuple val(genome), val(fasta), path("bwamem2_index.ready"), emit: index

    script:
    """
    bash ${shellQuote(index_helper)} \
        --genome ${shellQuote(genome)} \
        --fasta ${shellQuote(fasta)} \
        --cache-dir ${shellQuote(ref_cache)} \
        --ready bwamem2_index.ready \
        --lock-timeout-seconds ${params.index_lock_timeout_seconds} \
        --stale-lock-seconds ${params.index_stale_lock_seconds}
    """
}

process CONCAT_PE_FASTQS {
    tag "${meta.sample}"
    publishDir "${params.outdir}/concat_fastqs", mode: 'copy', enabled: params.publish_concat_fastqs

    input:
    tuple val(meta), path(r1s), path(r2s)

    output:
    tuple val(meta),
        path("${meta.sample}.R1.fastq.gz"),
        path("${meta.sample}.R2.fastq.gz"), emit: reads

    script:
    def r1_files = r1s.collect { shellQuote(it) }.join(' ')
    def r2_files = r2s.collect { shellQuote(it) }.join(' ')
    """
    cat ${r1_files} > ${shellQuote("${meta.sample}.R1.fastq.gz")}
    cat ${r2_files} > ${shellQuote("${meta.sample}.R2.fastq.gz")}
    """
}

process CONCAT_SE_FASTQS {
    tag "${meta.sample}"
    publishDir "${params.outdir}/concat_fastqs", mode: 'copy', enabled: params.publish_concat_fastqs

    input:
    tuple val(meta), path(r1s)

    output:
    tuple val(meta), path("${meta.sample}.fastq.gz"), emit: reads

    script:
    def r1_files = r1s.collect { shellQuote(it) }.join(' ')
    """
    cat ${r1_files} > ${shellQuote("${meta.sample}.fastq.gz")}
    """
}

process FASTP_PE {
    tag "${meta.sample}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.fastp.*", overwrite: true, failOnError: true
    publishDir "${params.outdir}/trimmed_fastqs", mode: 'copy', pattern: "*.trimmed.*.fastq.gz", enabled: params.publish_trimmed_fastqs

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta),
        path("${meta.sample}.trimmed.R1.fastq.gz"),
        path("${meta.sample}.trimmed.R2.fastq.gz"), emit: reads
    path "${meta.sample}.fastp.html", emit: html
    path "${meta.sample}.fastp.json", emit: json

    script:
    def low_complexity = (
        selected_analysis == 'digenome' &&
        !keep_multimappers &&
        params.fastp_low_complexity_filter
    ) ? "--low_complexity_filter --complexity_threshold ${params.fastp_complexity_threshold}" : ""
    """
    fastp \
        --in1 ${shellQuote(read1)} \
        --in2 ${shellQuote(read2)} \
        --out1 ${shellQuote("${meta.sample}.trimmed.R1.fastq.gz")} \
        --out2 ${shellQuote("${meta.sample}.trimmed.R2.fastq.gz")} \
        --html ${shellQuote("${meta.sample}.fastp.html")} \
        --json ${shellQuote("${meta.sample}.fastp.json")} \
        --thread ${task.cpus} \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --qualified_quality_phred ${params.fastp_qualified_quality_phred} \
        --length_required ${params.fastp_length_required} \
        ${low_complexity} \
        ${params.fastp_extra_args}
    """
}

process FASTP_SE {
    tag "${meta.sample}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.fastp.*", overwrite: true, failOnError: true
    publishDir "${params.outdir}/trimmed_fastqs", mode: 'copy', pattern: "*.trimmed.fastq.gz", enabled: params.publish_trimmed_fastqs

    input:
    tuple val(meta), path(read1)

    output:
    tuple val(meta), path("${meta.sample}.trimmed.fastq.gz"), emit: reads
    path "${meta.sample}.fastp.html", emit: html
    path "${meta.sample}.fastp.json", emit: json

    script:
    def low_complexity = (
        selected_analysis == 'digenome' &&
        !keep_multimappers &&
        params.fastp_low_complexity_filter
    ) ? "--low_complexity_filter --complexity_threshold ${params.fastp_complexity_threshold}" : ""
    """
    fastp \
        --in1 ${shellQuote(read1)} \
        --out1 ${shellQuote("${meta.sample}.trimmed.fastq.gz")} \
        --html ${shellQuote("${meta.sample}.fastp.html")} \
        --json ${shellQuote("${meta.sample}.fastp.json")} \
        --thread ${task.cpus} \
        --trim_poly_g \
        --trim_poly_x \
        --qualified_quality_phred ${params.fastp_qualified_quality_phred} \
        --length_required ${params.fastp_length_required} \
        ${low_complexity} \
        ${params.fastp_extra_args}
    """
}

process ALIGN_MARKDUP_PE {
    tag "${meta.sample}"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bam*", overwrite: true, failOnError: true
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.txt", overwrite: true, failOnError: true

    input:
    tuple val(meta), path(read1), path(read2), val(index_prefix)

    output:
    tuple val(meta),
        path("${meta.sample}.sorted.markdup.bam"),
        path("${meta.sample}.sorted.markdup.bam.bai"), emit: bam
    path "${meta.sample}.flagstat.txt", emit: flagstat
    path "${meta.sample}.stats.txt", emit: stats
    path "${meta.sample}.markdup.metrics.txt", emit: markdup_metrics

    script:
    def multimapper_opt = keep_multimappers ? "-a" : ""
    int total_cpus = task.cpus as int
    if (total_cpus < 2) {
        error "ALIGN_MARKDUP_PE requires at least 2 CPUs"
    }
    int piped_sort_threads = total_cpus >= 4 ? 2 : 0
    int alignment_threads = Math.max(
        1,
        total_cpus - piped_sort_threads - 1
    )
    int samtools_threads = Math.max(0, total_cpus - 1)
    """
    RG='@RG\\tID:${meta.sample}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}'

    bwa-mem2 mem \
        ${multimapper_opt} \
        -t ${alignment_threads} \
        -R "\${RG}" \
        ${shellQuote(index_prefix)} \
        ${shellQuote(read1)} \
        ${shellQuote(read2)} \
      | samtools sort -@ ${piped_sort_threads} -n \
            -o ${shellQuote("${meta.sample}.namesort.bam")} -

    samtools fixmate -@ ${samtools_threads} -m \
        ${shellQuote("${meta.sample}.namesort.bam")} \
        ${shellQuote("${meta.sample}.fixmate.bam")}
    samtools sort -@ ${samtools_threads} \
        -o ${shellQuote("${meta.sample}.positionsort.bam")} \
        ${shellQuote("${meta.sample}.fixmate.bam")}
    samtools markdup -@ ${samtools_threads} -s \
        ${shellQuote("${meta.sample}.positionsort.bam")} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")} \
        2> ${shellQuote("${meta.sample}.markdup.metrics.txt")}
    samtools index -@ ${samtools_threads} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")}
    samtools flagstat -@ ${samtools_threads} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")} \
        > ${shellQuote("${meta.sample}.flagstat.txt")}
    samtools stats -@ ${samtools_threads} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")} \
        > ${shellQuote("${meta.sample}.stats.txt")}
    """
}

process ALIGN_MARKDUP_SE {
    tag "${meta.sample}"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bam*", overwrite: true, failOnError: true
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.txt", overwrite: true, failOnError: true

    input:
    tuple val(meta), path(read1), val(index_prefix)

    output:
    tuple val(meta),
        path("${meta.sample}.sorted.markdup.bam"),
        path("${meta.sample}.sorted.markdup.bam.bai"), emit: bam
    path "${meta.sample}.flagstat.txt", emit: flagstat
    path "${meta.sample}.stats.txt", emit: stats
    path "${meta.sample}.markdup.metrics.txt", emit: markdup_metrics

    script:
    def multimapper_opt = keep_multimappers ? "-a" : ""
    int total_cpus = task.cpus as int
    if (total_cpus < 2) {
        error "ALIGN_MARKDUP_SE requires at least 2 CPUs"
    }
    int piped_sort_threads = total_cpus >= 4 ? 2 : 0
    int alignment_threads = Math.max(
        1,
        total_cpus - piped_sort_threads - 1
    )
    int samtools_threads = Math.max(0, total_cpus - 1)
    """
    RG='@RG\\tID:${meta.sample}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}'

    bwa-mem2 mem \
        ${multimapper_opt} \
        -t ${alignment_threads} \
        -R "\${RG}" \
        ${shellQuote(index_prefix)} \
        ${shellQuote(read1)} \
      | samtools sort -@ ${piped_sort_threads} \
            -o ${shellQuote("${meta.sample}.positionsort.bam")} -

    samtools markdup -@ ${samtools_threads} -s \
        ${shellQuote("${meta.sample}.positionsort.bam")} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")} \
        2> ${shellQuote("${meta.sample}.markdup.metrics.txt")}
    samtools index -@ ${samtools_threads} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")}
    samtools flagstat -@ ${samtools_threads} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")} \
        > ${shellQuote("${meta.sample}.flagstat.txt")}
    samtools stats -@ ${samtools_threads} \
        ${shellQuote("${meta.sample}.sorted.markdup.bam")} \
        > ${shellQuote("${meta.sample}.stats.txt")}
    """
}

process PLAN_CLEAVAGE_CHUNKS {
    tag "${meta.sample}"
    publishDir "${params.outdir}/pipeline_info/cleavage_chunks",
        mode: 'copy',
        overwrite: true,
        failOnError: true,
        pattern: "*.cleavage_chunks.tsv"

    input:
    tuple val(meta), path(bam), path(bai)
    path planner
    path blacklist_helper
    path genome_blacklist
    val blacklist_enabled
    val chunk_count
    val chunk_padding

    output:
    tuple val(meta), path("chunk_*.intervals.tsv"), emit: chunks
    path "${meta.sample}.cleavage_chunks.tsv", emit: plan

    script:
    def blacklist_args = blacklist_enabled ?
        "--genome-blacklist ${shellQuote(genome_blacklist)}" : ''
    """
    python3 ${shellQuote(planner)} \
        --bam ${shellQuote(bam)} \
        --chunks ${chunk_count} \
        --padding ${chunk_padding} \
        ${blacklist_args} \
        --output-dir . \
        --plan ${shellQuote("${meta.sample}.cleavage_chunks.tsv")}
    """
}

process CLEAVAGE_CALL_CHUNK {
    tag "${meta.sample}:${selected_analysis}:${chunk_id}"

    input:
    tuple val(meta),
        path(bam), path(bai), path(control_bam), path(control_bai),
        path(variant_vcf), path(variant_index),
        val(chunk_id), path(intervals_file)
    path caller
    path blacklist_helper
    path genome_blacklist
    val blacklist_enabled

    output:
    tuple val(meta),
        path("${meta.sample}.${chunk_id}.${selected_analysis}.raw.jsonl.gz"),
        path("${meta.sample}.${chunk_id}.${selected_analysis}.chunk.json"), emit: chunks

    script:
    def raw_output = "${meta.sample}.${chunk_id}.${selected_analysis}.raw.jsonl.gz"
    def summary_output = "${meta.sample}.${chunk_id}.${selected_analysis}.chunk.json"
    def control_args = meta.has_control ? [
        "--control-bam ${shellQuote(control_bam)}",
        "--control-sample ${shellQuote(meta.control_sample)}"
    ].join(' ') : ''
    def variant_args = meta.has_variant ?
        "--variant-vcf ${shellQuote(variant_vcf)}"
        : ''
    def multimapper_args = keep_multimappers ? '--keep-multimappers' : ''
    def blacklist_args = blacklist_enabled ?
        "--genome-blacklist ${shellQuote(genome_blacklist)}" : ''
    """
    python3 ${shellQuote(caller)} \
        --analysis ${shellQuote(selected_analysis)} \
        --bam ${shellQuote(bam)} \
        --sample ${shellQuote(meta.sample)} \
        --output-prefix ${shellQuote(meta.sample)} \
        ${control_args} \
        ${variant_args} \
        ${multimapper_args} \
        ${blacklist_args} \
        --intervals-file ${shellQuote(intervals_file)} \
        --chunk-id ${shellQuote(chunk_id)} \
        --raw-output ${shellQuote(raw_output)} \
        --summary-output ${shellQuote(summary_output)} \
        --ndigenome-min-count ${params.ndigenome_min_count} \
        --ndigenome-min-fraction ${params.ndigenome_min_fraction} \
        --ndigenome-min-mapq ${effective_ndigenome_min_mapq} \
        --ndigenome-opposite-window ${params.ndigenome_opposite_window} \
        --ndigenome-ambiguous-min-count ${params.ndigenome_ambiguous_min_count} \
        --ndigenome-ambiguous-min-fraction ${params.ndigenome_ambiguous_min_fraction} \
        --digenome-overhang ${params.digenome_overhang} \
        --digenome-pair-window ${params.digenome_pair_window} \
        --digenome-min-mapq ${effective_digenome_min_mapq} \
        --digenome-forward-cutoff ${params.digenome_forward_cutoff} \
        --digenome-reverse-cutoff ${params.digenome_reverse_cutoff} \
        --digenome-depth-cutoff ${params.digenome_depth_cutoff} \
        --digenome-fraction-cutoff ${params.digenome_fraction_cutoff} \
        --digenome-score-cutoff ${params.digenome_score_cutoff} \
        --artifact-window ${params.cleavage_artifact_window} \
        --max-softclip-fraction ${params.cleavage_max_softclip_fraction} \
        --max-indel-fraction ${params.cleavage_max_indel_fraction} \
        --min-support-mean-mapq ${effective_min_support_mean_mapq} \
        --control-min-depth ${params.cleavage_control_min_depth} \
        --control-max-fraction ${params.cleavage_control_max_fraction} \
        --control-min-fold ${params.cleavage_control_min_fold} \
        --control-max-q ${params.cleavage_control_max_q}
    """
}

process FINALIZE_CLEAVAGE_CALL {
    tag "${meta.sample}:${selected_analysis}"
    publishDir "${params.outdir}/${selected_analysis}",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(meta), path(raw_fragments), path(chunk_summaries)
    path finalizer
    path caller
    path blacklist_helper

    output:
    tuple val(meta),
        path("${meta.sample}.${selected_analysis}.all.tsv"), emit: all
    path "${meta.sample}.${selected_analysis}.high_confidence.tsv", emit: high_confidence
    path "${meta.sample}.${selected_analysis}.bed", emit: bed
    path "${meta.sample}.${selected_analysis}.qc.json", emit: qc
    path "${meta.sample}.${selected_analysis}_mqc.tsv", emit: multiqc

    script:
    def raw_args = raw_fragments.collect {
        "--raw-fragment ${shellQuote(it)}"
    }.join(' ')
    def summary_args = chunk_summaries.collect {
        "--chunk-summary ${shellQuote(it)}"
    }.join(' ')
    """
    python3 ${shellQuote(finalizer)} \
        --analysis ${shellQuote(selected_analysis)} \
        --sample ${shellQuote(meta.sample)} \
        --output-prefix ${shellQuote(meta.sample)} \
        ${raw_args} \
        ${summary_args}
    """
}

process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true, failOnError: true

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . --filename multiqc_report.html

    if [[ -d multiqc_report_data && ! -d multiqc_data ]]; then
        mv multiqc_report_data multiqc_data
    fi
    [[ -d multiqc_data ]] || {
        echo "ERROR: MultiQC data directory was not created" >&2
        exit 1
    }
    """
}

workflow {
    validator_ch = Channel.value(file("${baseDir}/bin/validate_samplesheet.py"))
    index_helper_ch = Channel.value(file("${baseDir}/bin/prepare_bwamem2_index.sh"))
    chunk_planner_ch = Channel.value(file("${baseDir}/bin/plan_cleavage_chunks.py"))
    blacklist_helper_ch = Channel.value(
        file("${baseDir}/bin/genome_blacklist.py")
    )
    cleavage_caller_ch = Channel.value(file("${baseDir}/bin/call_cleavage.py"))
    cleavage_finalizer_ch = Channel.value(
        file("${baseDir}/bin/finalize_cleavage_chunks.py")
    )
    samplesheet_ch = Channel.fromPath(params.input, checkIfExists: true)

    RUN_INFO(Channel.value(run_info_json))
    VALIDATE_SAMPLESHEET(
        samplesheet_ch,
        validator_ch,
        Channel.value(selected_analysis)
    )

    rows_ch = VALIDATE_SAMPLESHEET.out.csv.splitCsv(header: true)
    read_layout_ch = rows_ch.branch {
        pe: it.fastq_2?.trim()
        se: !it.fastq_2?.trim()
    }

    pe_grouped_ch = read_layout_ch.pe
        .map { row ->
            def meta = [
                sample: row.sample as String,
                is_control: row.is_control == 'true',
                control_sample: row.control as String,
                variant_vcf: row.variant_vcf as String,
                variant_index: row.variant_index as String
            ]
            tuple(
                meta,
                file(row.fastq_1 as String),
                file(row.fastq_2 as String)
            )
        }
        .groupTuple()
        .map { meta, r1s, r2s ->
            tuple(meta, r1s, r2s)
        }

    se_grouped_ch = read_layout_ch.se
        .map { row ->
            def meta = [
                sample: row.sample as String,
                is_control: row.is_control == 'true',
                control_sample: row.control as String,
                variant_vcf: row.variant_vcf as String,
                variant_index: row.variant_index as String
            ]
            tuple(meta, file(row.fastq_1 as String))
        }
        .groupTuple()
        .map { meta, r1s ->
            tuple(meta, r1s)
        }

    genome_ch = Channel.value(
        tuple(selected_genome, selected_fasta, params.ref_cache as String)
    )
    PREPARE_BWAMEM2_INDEX(genome_ch, index_helper_ch)

    resolved_index_ch = PREPARE_BWAMEM2_INDEX.out.index.map {
        genome, fasta, ready ->
            def values = ready.text.readLines()
                .findAll { it.contains('\t') }
                .collectEntries { line ->
                    def fields = line.split('\t', 2)
                    [(fields[0]): fields[1]]
                }
            values.index_prefix as String
    }

    CONCAT_PE_FASTQS(pe_grouped_ch)
    CONCAT_SE_FASTQS(se_grouped_ch)
    FASTP_PE(CONCAT_PE_FASTQS.out.reads)
    FASTP_SE(CONCAT_SE_FASTQS.out.reads)

    pe_align_ch = FASTP_PE.out.reads
        .combine(resolved_index_ch)
        .map {
            meta, read1, read2, index_prefix ->
                tuple(meta, read1, read2, index_prefix)
        }

    se_align_ch = FASTP_SE.out.reads
        .combine(resolved_index_ch)
        .map {
            meta, read1, index_prefix ->
                tuple(meta, read1, index_prefix)
        }

    ALIGN_MARKDUP_PE(pe_align_ch)
    ALIGN_MARKDUP_SE(se_align_ch)

    qc_ch = FASTP_PE.out.html
        .mix(FASTP_PE.out.json)
        .mix(FASTP_SE.out.html)
        .mix(FASTP_SE.out.json)
        .mix(ALIGN_MARKDUP_PE.out.flagstat)
        .mix(ALIGN_MARKDUP_PE.out.stats)
        .mix(ALIGN_MARKDUP_PE.out.markdup_metrics)
        .mix(ALIGN_MARKDUP_SE.out.flagstat)
        .mix(ALIGN_MARKDUP_SE.out.stats)
        .mix(ALIGN_MARKDUP_SE.out.markdup_metrics)

    no_control_bam = file("${baseDir}/assets/NO_CONTROL.bam")
    no_control_bai = file("${baseDir}/assets/NO_CONTROL.bam.bai")
    no_variant_vcf = file("${baseDir}/assets/NO_VARIANT.vcf.gz")
    no_variant_index = file("${baseDir}/assets/NO_VARIANT.vcf.gz.tbi")
    no_genome_blacklist = file("${baseDir}/assets/NO_BLACKLIST.bed")
    genome_blacklist_ch = Channel.value(
        has_genome_blacklist ?
            file(genome_blacklist_path) : no_genome_blacklist
    )

    aligned_bam_ch = ALIGN_MARKDUP_PE.out.bam.mix(ALIGN_MARKDUP_SE.out.bam)
    bam_type_ch = aligned_bam_ch.branch {
        control: it[0].is_control
        treated: !it[0].is_control
    }
    control_bam_ch = bam_type_ch.control
        .map { meta, bam, bai ->
            tuple(meta.sample, bam, bai)
        }

    treated_control_ch = bam_type_ch.treated.branch {
        matched: it[0].control_sample
        uncontrolled: !it[0].control_sample
    }

    matched_requests_ch = treated_control_ch.matched
        .map { meta, bam, bai ->
            tuple(meta.control_sample, meta, bam, bai)
        }
        .combine(control_bam_ch, by: 0)
        .map {
            control_sample, meta, bam, bai, control_bam, control_bai ->
                def request_meta = meta + [
                    has_control: true,
                    has_variant: meta.variant_vcf != ''
                ]
                tuple(
                    request_meta, bam, bai,
                    control_bam, control_bai,
                    meta.variant_vcf ?
                        file(meta.variant_vcf as String) : no_variant_vcf,
                    meta.variant_index ?
                        file(meta.variant_index as String) : no_variant_index
                )
        }

    uncontrolled_requests_ch = treated_control_ch.uncontrolled
        .map { meta, bam, bai ->
                def request_meta = meta + [
                    has_control: false,
                    has_variant: meta.variant_vcf != ''
                ]
                tuple(
                    request_meta, bam, bai,
                    no_control_bam, no_control_bai,
                    meta.variant_vcf ?
                        file(meta.variant_vcf as String) : no_variant_vcf,
                    meta.variant_index ?
                        file(meta.variant_index as String) : no_variant_index
                )
        }

    cleavage_requests_ch = matched_requests_ch.mix(uncontrolled_requests_ch)
    chunk_plan_requests_ch = cleavage_requests_ch.map {
        meta, bam, bai, control_bam, control_bai,
        variant_vcf, variant_index ->
            tuple(meta, bam, bai)
    }
    PLAN_CLEAVAGE_CHUNKS(
        chunk_plan_requests_ch,
        chunk_planner_ch,
        blacklist_helper_ch,
        genome_blacklist_ch,
        Channel.value(has_genome_blacklist),
        Channel.value(cleavage_chunks),
        Channel.value(cleavage_chunk_padding)
    )

    chunk_requests_ch = cleavage_requests_ch
        .join(PLAN_CLEAVAGE_CHUNKS.out.chunks, by: 0)
        .flatMap {
            meta, bam, bai, control_bam, control_bai,
            variant_vcf, variant_index,
            chunk_files ->
                def files = chunk_files instanceof List ?
                    chunk_files : [chunk_files]
                files.sort { it.name }.collect { intervals_file ->
                    def chunk_id = intervals_file.name.replace(
                        '.intervals.tsv',
                        ''
                    )
                    tuple(
                        meta, bam, bai,
                        control_bam, control_bai,
                        variant_vcf, variant_index,
                        chunk_id, intervals_file
                    )
                }
        }

    CLEAVAGE_CALL_CHUNK(
        chunk_requests_ch,
        cleavage_caller_ch,
        blacklist_helper_ch,
        genome_blacklist_ch,
        Channel.value(has_genome_blacklist)
    )
    finalize_requests_ch = CLEAVAGE_CALL_CHUNK.out.chunks.groupTuple()
    FINALIZE_CLEAVAGE_CALL(
        finalize_requests_ch,
        cleavage_finalizer_ch,
        cleavage_caller_ch,
        blacklist_helper_ch
    )
    qc_ch = qc_ch.mix(FINALIZE_CLEAVAGE_CALL.out.multiqc)

    MULTIQC(qc_ch.collect())
}
