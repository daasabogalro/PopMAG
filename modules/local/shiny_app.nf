process LAUNCH_SHINY_APP {
    tag "${meta}"
    label 'process_medium'
    errorStrategy 'ignore' 

    container 'docker.io/daasabogalro/popmag-dashboard:latest'
    containerOptions '-p 3838:3838'

    input:
    path(SNVs_summary)
    path(metadata)
    path(profile)
    path(merged_reports)
    path(fst_files)
    
    when:
    task.ext.when == null || task.ext.when

    script:
    
    def args = task.ext.args ?: ''
    def timeout = params.shiny_timeout ?: 1500
    println "To access the shiny app please open 0.0.0.0:3838 or localhost:3838 in your browser"
    """
    mkdir -p /srv/shiny-server/app/
    cp *.tsv /srv/shiny-server/app/
    cp *.txt /srv/shiny-server/app/
    cp intradiv* /srv/shiny-server/app/
    if [ -f "${metadata}" ]; then
        cp ${metadata} /srv/shiny-server/app/metadata.csv
    fi
    
    timeout ${timeout} Rscript -e "shiny::runApp('/srv/shiny-server/app', host='0.0.0.0', port=3838)"
    """
}
