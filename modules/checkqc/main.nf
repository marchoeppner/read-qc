process CHECKQC {
    tag "$meta.sample_id"
    label 'short_parallel'

    container "community.wave.seqera.io/library/python_numpy_pip_checkqc_interop:b5301d9801b8e66b"

    input:
    tuple val(meta),path(run_folder)

    output:
    path("*.json")          , emit: json
    path  'versions.yml'    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CheckQC module does not support Conda yet. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    report = prefix + '.checkqc.json'

    """
    checkqc $args --json  $run_folder > $report || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkqc: \$( checkqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
