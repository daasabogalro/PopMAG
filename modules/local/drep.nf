process DREP {
    tag "${meta.id}"
    label 'process_medium'    

    conda "bioconda::drep=3.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0':
    'biocontainers/drep:3.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(checkm2_report)

    output:
    tuple val(meta), path("dRep_${meta.id}"), emit: drep_output
    tuple val(meta), path("${meta.id}_concatenated.fa"), emit: concatenated_genome
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "dRep_${meta.id}"
    """
    if [ \$(ls -1 ${fasta} | wc -l) -gt 1 ]; then
        dRep dereplicate ${prefix} \
            -g ${fasta}/* \
            --genomeInfo ${checkm2_report} \
            -p ${task.cpus} \
            ${args}
    else
        mkdir -p ${prefix}/dereplicated_genomes
        cp ${fasta}/* ${prefix}/dereplicated_genomes/
    fi

    # Concatenate the dereplicated genomes
    cat ${prefix}/dereplicated_genomes/*.fa > ${meta.id}_concatenated.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drep: \$(dRep -h 2>&1 | head -2 | tail -1 | sed 's/^dRep v//')
    END_VERSIONS
    """
}
