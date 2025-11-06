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
    tuple val(meta), path("*.IS")                               , emit: profile
    tuple val(meta), path("*.IS/output/*.IS_SNVs.tsv")          , emit: snvs
    tuple val(meta), path("*.IS/output/*.IS_gene_info.tsv")     , emit: gene_info       , optional: true
    tuple val(meta), path("*.IS/output/*.IS_genome_info.tsv")   , emit: genome_info
    tuple val(meta), path("*.IS/output/*.IS_linkage.tsv*")       , emit: linkage
    tuple val(meta), path("*.IS/output/*.IS_mapping_info.tsv")  , emit: mapping_info
    tuple val(meta), path("*.IS/output/*.IS_scaffold_info.tsv") , emit: scaffold_info
    tuple val(meta), path("${meta.id}_concatenated.fna")        , emit: genome_db
    path "versions.yml"                                         , emit: versions

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
        -g ${meta.id}_concatenated.fna \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(inStrain -h | head -2 | sed 's/inStrain //g')
    END_VERSIONS
    """
}
