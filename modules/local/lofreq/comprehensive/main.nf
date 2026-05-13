process LOFREQ_COMPREHENSIVE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4' :
        'biocontainers/lofreq:2.1.5--py38h588ecb2_4' }"

    input:
    tuple val(meta) , path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_filtered.vcf.gz")      , emit: vcf
    tuple val(meta), path("*_filtered.vcf.gz.tbi")  , emit: tbi
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def platform = "$meta.seq_chemistry"
    def args = task.ext.args ?: ''
    def args_indelqual = task.ext.args2 ?: ''
    def filter_args = task.ext.args3 ?: ''
    
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_intervals = intervals ? "-l ${intervals}" : ""

    def alignment_cram =  bam.Extension == "cram" ? true : false
    def alignment_bam = bam.Extension == "bam" ? true : false
    def alignment_out = alignment_cram ? bam.BaseName + ".bam" : "${bam}"

    def samtools_cram_convert = ''
    samtools_cram_convert += alignment_cram ? "    samtools view -T ${fasta} ${bam} -@ $task.cpus -o ${alignment_out}\n" : ''
    samtools_cram_convert += alignment_cram ? "    samtools index ${alignment_out}\n" : ''

    def samtools_cram_remove = ''
    samtools_cram_remove += alignment_cram ? "    rm ${alignment_out}\n" : ''
    samtools_cram_remove += alignment_cram ? "    rm ${alignment_out}.bai\n " : ''
    
    // don't want realigment with element chemistry
    if (meta.seq_chemistry == 'element') {
        """
        $samtools_cram_convert

        # add correct tags for indels and index
        lofreq \\
            indelqual \\
            $args_indelqual \\
            -f $fasta \\
            -o ${prefix}_lofreq.bam \\
            $alignment_out
     
            samtools index ${prefix}_lofreq.bam
    
        # variant calling
        lofreq \\
            call-parallel \\
            --pp-threads $task.cpus \\
            $args \\
            $options_intervals \\
            -f $fasta \\
            -o ${prefix}_unfiltered.vcf.gz \\
            ${prefix}_lofreq.bam

        # variant filtering via lofreq filter function

        lofreq \\
            filter \\
            $filter_args \\
            -o ${prefix}_filtered.vcf.gz \\
            -i ${prefix}_unfiltered.vcf.gz

        $samtools_cram_remove

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
        END_VERSIONS
        """

    } else {
        """
        $samtools_cram_convert

        # Realign reads and pipe directly to samtools sort
        lofreq \\
            viterbi \\
            -f $fasta \\
            $alignment_out \\
            | samtools sort -o ${prefix}_lofreq_realign.bam -

        # add correct tags for indels and index
        lofreq \\
            indelqual \\
            $args_indelqual \\
            -f $fasta \\
            -o ${prefix}_lofreq.bam \\
            ${prefix}_lofreq_realign.bam

        samtools index ${prefix}_lofreq.bam

        # variant calling
        lofreq \\
            call-parallel \\
            --pp-threads $task.cpus \\
            $args \\
            $options_intervals \\
            -f $fasta \\
            -o ${prefix}_unfiltered.vcf.gz \\
            ${prefix}_lofreq.bam

        # variant filtering via lofreq filter function

        lofreq \\
            filter \\
            $filter_args \\
            -o ${prefix}_filtered.vcf.gz \\
            -i ${prefix}_unfiltered.vcf.gz

        $samtools_cram_remove

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
