process READQC_JSON {

    label 'serial_short'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.json") , emit: json
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}" 

    """
    readqc_json.py --run_name ${meta.id} \\
    --output ${prefix}.readqc.json \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
    END_VERSIONS
    """
  
}
