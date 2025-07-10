process SINGLEM_METAPACKAGE {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/singlem:0.18.3--pyhdfd78af_0' :
    'biocontainers/singlem:0.19.0--pyhdfd78af_0' }"

    output:
    path("*.smpkg.zb"), emit: singleM_metapackage

    script:
    def args = task.ext.args ?: ''
    """
    singlem data --output-directory .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singleM: \$(singlem pipe --version)
    END_VERSIONS
    """
}