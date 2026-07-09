// Creates sample specific json file to be uploaded to SQVD
// Curretly imports fastqc, samtools stats & readgroupsbyumi data 

process EXPORT_TO_JSON_SQVD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yaml"
    container "seglh/bw_ngstools_java:1.0.1"

    publishDir "${params.outdir}/sqvd_json", mode: 'copy'

    input:
    tuple val(meta), path(files1), path(files2), path(files3)

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path  "versions.yml"            , emit: versions

    script:
    """
    aggregate_metrics.py --input . --output ${meta.id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        json_sqvd: python_env: \$(echo \$(/usr/bin/env python --version 2>&1) | sed -e "s/Python //g" )
    END_VERSIONS
    """
}