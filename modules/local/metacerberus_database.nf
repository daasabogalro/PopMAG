process METACERBERUS_DATABASEDOWNLOAD {
    label 'process_high'

    conda "bioconda::metacerberus==1.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/metacerberus:1.4.0--pyhdfd78af_1' :
    'biocontainers/metacerberus:1.4.0--pyhdfd78af_1' }"

    output:
    tuple val(meta), path("metacerberus_db"), emit: database
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta           = [id: 'metacerberus_db']
    // Only fetch the download group that METACERBERUS_ANNOTATION will search, instead of
    // all 13 groups. Note that databases.tsv keys rows by *group* (KOFam, CAZy, ...) while
    // --hmm takes the per-file name metacerberus.py derives from the filename
    // (KOFam_prokaryote, KOFam_eukaryote, KOFam_all, CAZy, ...). A group name is always a
    // prefix of the file names it contains, and no group name is a prefix of another, so
    // matching the group as a prefix of the selection resolves both spellings.
    def selected   = params.annotation_db ?: 'ALL'
    def db_filter  = selected == 'ALL' ? '1' : "index(\"${selected}\", \$1) == 1"
    """
    # --download writes databases.tsv but its own fetching is unreliable against the
    # current mirrors, so the files are retrieved explicitly below. Resolve the package
    # directory rather than hardcoding a Python version.
    metacerberus.py --download || true

    DB_DIR=\$(python3 -c "import meta_cerberus, os; print(os.path.join(os.path.dirname(meta_cerberus.__file__), 'DB'))")
    mkdir -p metacerberus_db
    cp -r "\$DB_DIR"/. metacerberus_db/

    tail -n +2 metacerberus_db/databases.tsv \\
        | awk -F'\\t' '${db_filter} { print \$2 "\\t" \$3 }' > db_files.tsv

    if [ ! -s db_files.tsv ]; then
        echo "ERROR: no database files matched '${selected}'. Valid values are the names" >&2
        echo "       in metacerberus_db/databases.tsv, or ALL." >&2
        exit 1
    fi

    # Read from a file rather than a pipe so the loop runs in this shell and the script's
    # `set -e` still applies to a failed download. The container ships BusyBox wget, which
    # has no retry option, so retries are explicit.
    while IFS=\$'\\t' read -r filename url; do
        echo "fetching \$filename"
        ok=0
        for attempt in 1 2 3; do
            if wget -q -T 60 -O "metacerberus_db/\$filename" "\$url"; then ok=1; break; fi
            echo "  attempt \$attempt failed, retrying"
            sleep 5
        done
        if [ "\$ok" != 1 ]; then
            echo "ERROR: could not download \$filename from \$url" >&2
            exit 1
        fi
    done < db_files.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metacerberus: \$(echo \$(metacerberus.py --version))
    END_VERSIONS
    """
}
