process MERGE_GENOME_REPORTS_INSTRAIN {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'docker.io/suyanan/python:pip-numpy-h5py-pandas-matplotlib-seaborn' }"

    input:
    tuple val(meta), path(genome_reports)
    
    output:
    tuple val(meta), path("*intradiv*"), emit: merged_instrain_genome_reports

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    
    python3 ${projectDir}/bin/combine_genome_info.py \\
        ${genome_reports} \\
        --output intradiv \\
        ${args}
    """
}