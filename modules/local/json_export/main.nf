process EXPORT_TO_JSON_SQVD {
    tag "${meta.sample}"
    label 'process_single'

    conda "${moduleDir}/environment.yaml"
    container "seglh/bw_ngstools_java:1.0.1"

    publishDir "${params.outdir}/sqvd_json", mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.sample}.json"), emit: json

    script:
    """
    python ${projectDir}/modules/local/json_export/resources/usr/bin/aggregate_metrics.py --input . --output ${meta.sample}.json
    """
}