process TRANSFORM_CHECKM2_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'biocontainers/python:3.9--1' }"
    
    input:
    tuple val(meta), path(checkm2_report)
    val(extension)

    output:
    tuple val(meta), path("${meta.id}_transformed_checkm2_report.csv"), emit: transformed_report

    script:
    """
    transform_checkm2_report.py \
        ${checkm2_report} \
        ${meta.id}_transformed_checkm2_report.csv \
        ${extension}
    """
}
