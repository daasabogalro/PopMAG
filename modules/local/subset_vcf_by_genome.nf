process SUBSET_VCF_BY_GENOME {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::pandas"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/python:3.9--1' :
    'docker.io/suyanan/python:pip-numpy-h5py-pandas-matplotlib-seaborn' }"

    input:
    tuple val(meta), path(vcf_file), path(contigs2bin_file)

    output:
    tuple val(meta), path("genome_vcfs"), emit: genome_vcfs

    script:
    """
    python3 ${projectDir}/bin/subset_vcf_by_genome.py \
        --vcf_file ${vcf_file} \
        --contigs2bin_file ${contigs2bin_file} \
        --out_dir genome_vcfs

    ls -lh genome_vcfs
    """
}