process SNVS_TO_VCF {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'docker.io/suyanan/python:pip-numpy-h5py-pandas-matplotlib-seaborn' }"

    input:
    tuple val(meta), path(snv_files), path(contigs2bin)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf_files
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def min_coverage = task.ext.min_coverage ?: 0
    def min_var_freq = task.ext.min_var_freq ?: 0
    def include_nonvariant = task.ext.include_nonvariant ?: false
    def vcf_prefix = task.ext.vcf_prefix ?: ''
    def vcf_suffix = task.ext.vcf_suffix ?: ''

    """
    python3 ${projectDir}/bin/snvs_to_vcf_by_reference.py \\
        --snv_files ${snv_files.join(' ')} \\
        --out_dir . \\
        --min_coverage ${min_coverage} \\
        --min_var_freq ${min_var_freq} \\
        --vcf_prefix "${vcf_prefix}" \\
        --vcf_suffix "${vcf_suffix}" \\
        ${include_nonvariant ? '--include_nonvariant' : ''} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        snv_converter: "1.0.0"
    END_VERSIONS
    """

    stub:
    """
    touch \${meta.id}_dummy.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        pandas: "1.5.3"
        snv_converter: "1.0.0"
    END_VERSIONS
    """
} 