process CHECKM2_DATABASEDOWNLOAD {
    label 'process_single'

    conda "bioconda::checkm2=1.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.2--pyh7cba7a3_0':
        'biocontainers/checkm2:1.0.2--pyh7cba7a3_0' }"

    input:

    output:
    tuple val(meta), path("checkm2_database.dmnd"), emit: database
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args        = task.ext.args ?: ''
    meta            = [id: 'checkm2_db']    

    """
    checkm2 database --download --path .

    db_path=\$(find -name *.dmnd)
    mv \$db_path checkm2_database.dmnd
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version 2>&1 | sed 's/^.*CheckM2 v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p checkm2_database
    touch checkm2_database/checkm2_database.dmnd
    touch checkm2_database/checkm2_database.sdb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: 1.0.1
    END_VERSIONS
    """
}
