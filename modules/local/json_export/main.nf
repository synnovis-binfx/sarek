process EXPORT_TO_JSON_SQVD {
    tag "${meta.sample}"
    label 'process_single'

    conda "${moduleDir}/environment.yaml"
    container "seglh/bw_ngstools_java:1.0.1"

    publishDir "${params.outdir}/sqvd_json", mode: 'copy'

    input:
    tuple val(meta), path(files1), path(files2), path(files3)

    output:
    tuple val(meta), path("${meta.sample}.json"), emit: json

    script:
    """
    aggregate_metrics.py --input . --output ${meta.sample}.json
    """
}