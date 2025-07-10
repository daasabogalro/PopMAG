process MERGE_REPORTS_METACERBERUS_INSTRAIN {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'docker.io/suyanan/python:pip-numpy-h5py-pandas-matplotlib-seaborn' }"

    input:
    tuple val(meta), path(bin_metrics), path(annotation_file)
    
    output:
    tuple val(meta), path("*.tsv"), emit: merged_reports
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_name = meta.bin ?: meta.id
    """
    
    python3 ${projectDir}/bin/merge_reports_metacerberus_instrain.py \\
        --bin_metrics ${bin_metrics} \\
        --annotations ${annotation_file} \\
        --output_dir . \\
        --genome_name ${genome_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        merge_reports_metacerberus_instrain: "1.0.0"
    END_VERSIONS
    """
} 