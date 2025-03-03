process BOWTIE2_INSTRAIN_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' :
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' }"

    input:
    tuple val(meta), path(concatenated_mags), path(index), val(reads_meta), path(reads)

    output:
    tuple val(meta), path(concatenated_mags), path("${reads_meta.id}_to_${meta.id}_genomes.bam"), path("${reads_meta.id}_to_${meta.id}_genomes.bam.bai"), val(reads_meta), emit: mappings
    tuple val(meta), val(reads_meta), path("*.bowtie2.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def name = "${reads_meta.id}_to_${meta.id}_genomes"
    def input = "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
    bowtie2 \\
        -p "${task.cpus}" \\
        -x \$INDEX \\
        $args \\
        $input \\
        2> "${name}.bowtie2.log" | \\
        samtools view -@ "${task.cpus}" -bS | \\
        samtools sort -@ "${task.cpus}" -o "${name}.bam"
    samtools index "${name}.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}