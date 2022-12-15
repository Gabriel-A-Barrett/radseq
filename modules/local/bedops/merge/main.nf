process BEDOPS_MERGE_BED {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bedops=2.4.41" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    bedops --merge *.bed \\
    ${args} \\
    > ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops: \$(bedops --version | sed -n '/version:/p' | cut -d' ' -f 5)
    END_VERSIONS
    """
}