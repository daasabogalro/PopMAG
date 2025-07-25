process POGENOM {
    tag "${meta.sample} - ${mag.baseName}"
    label 'process_low'

    container 'docker.io/daasabogalro/pogenom:latest'

    input:
    tuple val(meta),  path(mag), path(vcf_file)    
    
    output:
    tuple val(meta), path("*"), emit: pogenom_results
    tuple val(meta), path("*.fst.txt"), emit: pogenom_fst
    //tuple val(meta), path("*.intradiv*"), emit: pogenom_intradiv
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    genome_size=\$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length(\$0)} END {print seqlen}' ${mag} | awk '{sum+=\$1} END {print sum}' )

    pogenom.pl \
    --vcf_file ${vcf_file} \
    --out ${mag.baseName} \
    --genome_size \${genome_size} \
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pogenom: "latest"
    END_VERSIONS
    """
}