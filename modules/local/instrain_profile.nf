process INSTRAIN_PROFILE {
    tag "${meta.id}_${reads_meta.id}"
    label 'process_high'

    conda "bioconda::instrain==1.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/instrain:1.9.0--pyhdfd78af_0' :
    'biocontainers/instrain:1.9.0--pyhdfd78af_0'}"

    input: 
    tuple val(meta), path(concatenated_mags), path(bam), path(bai), val(reads_meta), path(genes), path(scaffolds2bin)
    
    output:
    tuple val(meta), val(reads_meta), path("${reads_meta.id}_to_${meta.id}_reps.IS"), emit: instrain_output
    tuple val(meta), path("${meta.id}_concatenated.fna"), emit: genome_db
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    zcat ${genes} > ${meta.id}_concatenated.fna
    
    inStrain profile \
        $bam \
        $concatenated_mags \
        -o ${reads_meta.id}_to_${meta.id}_reps.IS \
        -p $task.cpus \
        -s $scaffolds2bin \
        #-g ${meta.id}_concatenated.fna \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(inStrain -v | head -2 | sed 's/inStrain //g')
    END_VERSIONS
    """
}