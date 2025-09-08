// pipeline params
params.raw_data_dir = '' // Set default to an empty string
params.reads_suffix = '_R{1,2}.fastq.gz' // Default pattern for paired-end files

// PROCESSES

process PRE_FASTQC {

    publishDir "${params.publishDir}/pre_fastqc", mode: 'copy'

    input:
    tuple val(sampleID), file(reads)

    output:
    path "${sampleID}_R1_fastqc.{zip,html}", emit: R1_report
    path "${sampleID}_R2_fastqc.{zip,html}", emit: R2_report

    shell:
    '''
        module load bio/fastqc/0.12.1 # MAY NEED TO UPDATE IF HYDRA CHANGES
        fastqc !{reads[0]}
        fastqc !{reads[1]}
    '''
}

process FASTP {

    publishDir "${params.publishDir}/fastp", mode: 'copy'

    input:
    tuple val(sampleID), file(reads)

    output:
    path "${reads[0].simpleName}_trimmed.fastq.gz", emit: R1_trimmed
    path "${reads[1].simpleName}_trimmed.fastq.gz", emit: R2_trimmed

    shell:
    '''
        module load bio/fastp/0.23.4 # MAY NEED TO UPDATE IF HYDRA CHANGES
        # USER CAN MODIFY FASTP OPTIONS BELOW IF DESIRED
        fastp -i !{reads[0]} \
            -I !{reads[1]} \
            -o !{reads[0].simpleName}_trimmed.fastq.gz \
            -O !{reads[1].simpleName}_trimmed.fastq.gz
    '''
}

process POST_FASTQC {

    publishDir "${params.publishDir}/post_fastqc", mode: 'copy'

    input:
    path R1_trimmed
    path R2_trimmed

    output:
    path "${R1_trimmed.simpleName}_fastqc.{zip,html}", emit: R1_trimmed_report
    path "${R2_trimmed.simpleName}_fastqc.{zip,html}", emit: R2_trimmed_report

    shell:
    '''
        module load bio/fastqc/0.12.1 # MAY NEED TO UPDATE IF HYDRA CHANGES
        fastqc !{R1_trimmed}
        fastqc !{R2_trimmed}
    '''
}

process MULTIQC {

    publishDir "${params.publishDir}/multiqc", mode: 'copy'

    input:
    path pre_reports
    path post_reports

    output:
    path "pretrimming.html"
    path "posttrimming.html"

    shell:
    '''
        module load bio/multiqc/1.23 # MAY NEED TO UPDATE IF HYDRA CHANGES
        multiqc !{pre_reports} -n pretrimming.html
        multiqc !{post_reports} -n posttrimming.html
    '''
}

// WORKFLOW

workflow QC {
    // Validate that the user has provided the raw_data_dir
    if (!params.raw_data_dir) {
        exit(1, "Error: Please specify the input directory using --raw_data_dir /path/to/your/reads")
    }

    // check for directory existence
    def rawDir = file(params.raw_data_dir)
    if (!rawDir.exists() || !rawDir.isDirectory()) {
        exit(1, "Error: The provided raw data directory does not exist or is not a directory.\nPath: ${params.raw_data_dir}")
    }

    // New check to ensure the glob pattern finds files
    def inputFileList = file("${params.raw_data_dir}/*${params.reads_suffix}").toList()
    if (inputFileList.isEmpty()) {
        exit(1, "Error: No files found matching the pattern '${params.reads_suffix}' in the directory '${params.raw_data_dir}'. Please check the path and pattern.")
    }

    // channels
    reads = Channel.fromFilePairs("${params.raw_data_dir}/*${params.reads_suffix}")

    // run FASTQC and trimming
    PRE_FASTQC(reads)
    FASTP(reads)
    POST_FASTQC(FASTP.out.R1_trimmed, FASTP.out.R2_trimmed)

    // collect FASTQC reports
    pre_reports = PRE_FASTQC.out.R1_report.mix(PRE_FASTQC.out.R2_report).collect()
    post_reports = POST_FASTQC.out.R1_trimmed_report.mix(POST_FASTQC.out.R2_trimmed_report).collect()

    // run MULTIQC
    MULTIQC(pre_reports, post_reports)
}