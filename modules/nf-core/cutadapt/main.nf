process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.2--py311haab0aaa_0' :
        'biocontainers/cutadapt:5.2--py311haab0aaa_0'}"

    input:
    tuple val(meta), path(reads)
    path(adapter_fasta_r1)
    path(adapter_fasta_r2)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    tuple val("${task.process}"), val("cutadapt"), eval('cutadapt --version'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def adapter_list_r1 = adapter_fasta_r1 ? "-a file:${adapter_fasta_r1}" : ""
    def adapter_list_r2 = adapter_fasta_r2 ? "-A file:${adapter_fasta_r2}" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    def trimmed_interleaved = "-o ${prefix}.trim.fastq.gz"
    
    if ( task.ext.args?.contains('--interleaved') ) {
    """
    cutadapt \\
        --cores $task.cpus \\
        $adapter_list_r1 $adapter_list_r2 \\
        $args \\
        $trimmed_interleaved \\
        $reads \\
        > ${prefix}.cutadapt.log
    """
    } else {
    """
    cutadapt \\
        --cores $task.cpus \\
        $adapter_list_r1 $adapter_list_r2 \\
        $args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    """
    }
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.trim.fastq.gz ;"
    }
    else if ( task.ext.args?.contains('--interleaved') ) {
        output_command = "echo '' | gzip > ${prefix}.trim.fastq.gz ;"
    }
    else {
        output_command  = "echo '' | gzip > ${prefix}_1.trim.fastq.gz ;"
        output_command += "echo '' | gzip > ${prefix}_2.trim.fastq.gz ;"
    }
    """
    ${output_command}
    touch ${prefix}.cutadapt.log
    """
}
