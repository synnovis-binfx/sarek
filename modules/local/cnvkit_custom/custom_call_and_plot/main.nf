process CNVKIT_CUSTOM_CALL_AND_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.10--pyhdfd78af_0':
        'biocontainers/cnvkit:0.9.10--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(cns), path(cnr), path(vcf)
    //path(cnvkit_plot_targets)

    output:
    tuple val(meta), path("*.cns"), emit: cns       , optional: true
    tuple val(meta), path("*.pdf"), emit: pdfs      , optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_cmd = vcf ? "-v $vcf" : ""
    """
    echo -e "figure.figsize: 18, 14\\naxes.titlepad: 50" >  matplotlibrc  && \\
        export MPLCONFIGDIR=${task.workDir}/tmp && \\
        cnvkit.py call \\
        $cns \\
        $args \\
        -o ${prefix}.cns \\
        && cnvkit.py scatter ${cnr} -s ${prefix}.cns -v ${vcf} -z 0 --min-variant-depth 50 -o ${prefix}.pdf 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e 's/cnvkit v//g')
    END_VERSIONS
    """
}
