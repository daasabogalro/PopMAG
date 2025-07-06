process METACERBERUS_DATABASEDOWNLOAD {
    label 'process_high'

    conda "bioconda::metacerberus==1.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_0' :
    'biocontainers/metacerberus:1.4.0--pyhdfd78af_1' }"

    input:

    output:
    tuple val(meta), path("metacerberus_db"), emit: database
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    additional_dbs = task.ext.additional_dbs ?: ''
    meta           = [id: 'metacerberus_db']
    """
    metacerberus.py --download \
        COG \
        #KOFam \
        #AMRFinder \
        $additional_dbs

    echo "MetaCerberus databases downloaded" > metacerberus_db
    #cp -r /usr/local/lib/python3.9/site-packages/meta_cerberus/DB ./metacerberus_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacerberus: \$(echo \$(metacerberus.py --version))
    END_VERSIONS
    """
}