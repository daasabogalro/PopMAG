process EXTRACT_BIN_METRICS {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'docker.io/suyanan/python:pip-numpy-h5py-pandas-matplotlib-seaborn' }"

    input:
    tuple val(meta), path(gene_metrics_files), path(contigs2bin_file)
    
    output:
    tuple val(meta), path("*_metrics.tsv"), emit: bin_metrics
    path("SNVs_${meta.id}_summary.tsv"), emit: combined_summary
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_id = meta.id
    """
    mkdir -p bin_metrics
    
    python3 ${projectDir}/bin/extract_bin_metrics.py \\
        ${contigs2bin_file} \\
        ${gene_metrics_files.join(' ')} \\
        -d . \\
        -s ${sample_id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        extract_bin_metrics: "1.0.0"
    END_VERSIONS
    """
} 