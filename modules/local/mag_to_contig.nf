process EXTRACT_CONTIG_NAMES {
    label 'process_single'
    
    input:
    tuple val(meta), path(mags)
    
    output:
    tuple val(meta), path("*combined_contigs.txt"), emit: contigs_to_bin_combined
    tuple val(meta), path("*individual_contigs.txt"), emit: contigs_to_bin_individual
    script:
    """
    ${projectDir}/bin/process_mags.py ${meta.id} ${mags}
    """
}