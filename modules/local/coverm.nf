process COVERM {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::coverm==0.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/coverm:0.7.0--hb4818e0_2' :
    'biocontainers/coverm:0.7.0--hb4818e0_2'}"

    input:
    tuple val(meta), path(contigs2bin), path(bams)

    output:
    tuple val(meta), path("*_coverage.txt"), emit: coverage

    script:
    """
    awk -F'\t' 'BEGIN {OFS = FS} {print \$2,\$1}' ${contigs2bin} > tr_c2b.txt

    coverm genome \
    --genome-definition tr_c2b.txt \
    --bam-files ${bams} \
    --methods trimmed_mean \
    --trim-min 0.1 \
    --trim-max 0.9 \
    --min-covered-fraction 0.5 \
    --threads $task.cpus \
    --output-file ${meta.id}_coverage.txt
    """
}
