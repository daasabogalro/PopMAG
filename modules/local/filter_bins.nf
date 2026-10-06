process FILTER_BINS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(checkm2_output_files), path(transformed_report)

    output:
    tuple val(meta), path("${meta.id}_filtered_bins"), emit: filtered_bins
    path "filtered_report_${meta.id}.csv", emit: filtered_report

    script:
    def completeness = params.min_completeness ?: 90
    def contamination = params.max_contamination ?: 5
    """
    filter_bins.py \
        "${checkm2_output_files}" \
        "${transformed_report}" \
        "${meta.id}" \
        ${completeness} \
        ${contamination}
    """
}