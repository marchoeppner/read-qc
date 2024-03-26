process BIOBLOOMTOOLS_CATEGORIZER {

    label 'short_parallel'

    tag "${meta.sample_id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biobloomtools:2.3.5--h4056dc3_2' :
        'quay.io/biocontainers/biobloomtools:2.3.5--h4056dc3_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    path('versions.yml'), emit: versions
    path("*summary.tsv"), emit: results

    script:

    """
    biobloomcategorizer -p $meta.sample_id -t ${task.cpus} -e -f "${params.references.bloomfilter}" $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BioBloomtools: \$(biobloomcategorizer -version 2>&1 | head -n1 | sed -e "s/.*) //g")
    END_VERSIONS

    """
}
