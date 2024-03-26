process CHECKQC {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkqc:4.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/fastqc:checkqc:4.0.1--pyhdfd78af_0' }"

    input:
    path(run_folder)

    output:
    path("*.json")          , emit: json
    path  'versions.yml'    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${params.run_name}"
    report = prefix + '.checkqc.json'
    """
    checkqc $args --json $report $run_folder
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkqc: \$( checkqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
