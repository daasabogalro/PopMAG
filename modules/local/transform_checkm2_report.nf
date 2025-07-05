process TRANSFORM_CHECKM2_REPORT {
    label 'process_single'
    
    input:
    tuple val(meta), path(checkm2_report)
    val(extension)

    output:
    tuple val(meta), path("${meta.id}_transformed_checkm2_report.csv"), emit: transformed_report

    script:
    """
    ${projectDir}/bin/transform_checkm2_report.py \
        ${checkm2_report} \
        ${meta.id}_transformed_checkm2_report.csv \
        ${extension}
    """
}
