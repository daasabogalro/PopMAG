process INSTRAIN_COMPARE {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::instrain==1.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/instrain:1.9.0--pyhdfd78af_0' :
    'biocontainers/instrain:1.9.0--pyhdfd78af_0'}"

    input: 
    //TODO modify the input channel to include the bams used in the profile step. They have to ve in the same order
    //as the input files 
    tuple val(meta), path(profiles)
    
    output:
    tuple val(meta), path("*.IS_compare")                                               , emit: compare
    //tuple val(meta), path("*.IS_compare/output/*.IS_compare_comparisonsTable.tsv")      , emit: comparisons_table   , optional: true
    //tuple val(meta), path("*.IS_compare/output/*.IS_compare_pooled_SNV_data.tsv")       , emit: pooled_snv
    //tuple val(meta), path("*.IS_compare/output/*.IS_compare_pooled_SNV_data_keys.tsv")  , emit: snv_keys
    //tuple val(meta), path("*.IS_compare/output/*.IS_compare_pooled_SNV_info.tsv")       , emit: snv_info
    //path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    inStrain compare \
        -i $profiles \\
        -o ${meta.id}.IS_compare \\
        --processes $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(inStrain -h | head -2 | sed 's/inStrain //g')
    END_VERSIONS
    """
}