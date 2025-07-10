process SINGLEM {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/singlem:0.18.3--pyhdfd78af_0' :
    'biocontainers/singlem:0.19.0--pyhdfd78af_0' }"

    
    input:
    tuple val(meta), path(genomes)
    path(metapackage)

    output:
    tuple val(meta), path("*.tsv"), emit: singleM_profile

    script:
    def args = task.ext.args ?: ''

    """
    singlem pipe \\
        --genome-fasta-files ${genomes} \\
        -p ${meta.id}_profile.tsv \\
        --threads ${task.cpus} \\
        --metapackage ${metapackage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singleM: \$(singlem pipe --version)
    END_VERSIONS
    """
}