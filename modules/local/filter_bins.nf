process FILTER_BINS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(checkm2_output_files), path(transformed_report)
    val min_completeness
    val max_contamination

    output:
    tuple val(meta), path("${meta.id}_filtered_bins"), emit: filtered_bins
    path "filtered_report_${meta.id}.csv", emit: filtered_report

    script:
    """
    filter_bins.py \
        "${checkm2_output_files}" \
        "${transformed_report}" \
        "${meta.id}" \
        ${min_completeness} \
        ${max_contamination}
    """
}