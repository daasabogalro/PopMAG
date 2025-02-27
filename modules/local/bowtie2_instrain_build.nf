process BOWTIE2_INSTRAIN_BUILD {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1' :
        'biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1' }"

    input:
    tuple val(meta), path(concatenated_mags)
    
    output:
    tuple val(meta), path(concatenated_mags), path("${meta.id}_index*"), emit: index

    script:
    """
    ## Build index
    bowtie2-build \
    --threads ${task.cpus} \
    ${concatenated_mags} \
    "${meta.id}_index"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}