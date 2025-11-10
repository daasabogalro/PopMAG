process METACERBERUS_ANNOTATION {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::metacerberus==1.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_0' :
    'biocontainers/metacerberus:1.4.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(prodigal_genes)
    tuple val(db_meta), path(DB)

    output:
    tuple val(meta), path("annotation/${meta.sample}/${meta.id}"), emit: annotation_dir
    tuple val(meta), path("annotation/${meta.sample}/${meta.id}/final/*/final_annotation_summary.tsv"), emit: final_annotation

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def database = params.annotation_db ?: ALL

    """
    gunzip -c ${prodigal_genes} > ${meta.id}.faa

    metacerberus.py --protein ${meta.id}.faa \\
    --hmm ${database} \\
    --illumina \\
    --meta \\
    --db-path ${DB} \\
    --dir-out annotation/${meta.sample}/${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacerberus: \$(echo \$(metacerberus.py --version))
    END_VERSIONS
    """
}
