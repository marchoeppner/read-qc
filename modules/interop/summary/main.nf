process INTEROP_SUMMARY {

    tag "${meta.sample_id}"
    
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/illumina-interop:1.5.0--h503566f_0' :
        'quay.io/biocontainers/illumina-interop:1.5.0--h503566f_0' }"

    input:
    tuple val(meta), path(folder)

    output:
    tuple val(meta), path('*interop.csv')  , emit: csv
    path('versions.yml')            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def summary = meta.sample_id + ".interop.csv"

    """
    interop_summary $args $folder > $summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        illumina-interop: \$(interop_summary -version 2>&1 | head -n  1 | sed -e "s/# Version: //g")
    END_VERSIONS
    """
   
}
