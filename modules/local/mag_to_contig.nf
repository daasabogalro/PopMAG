process EXTRACT_CONTIG_NAMES {
    label 'process_single'
    
    input:
    tuple val(meta), path(mags)
    
    output:
    tuple val(meta), path("${meta.id}_combined_contigs.txt"), emit: scaffold_to_bin_combined
    
    script:
    """
    process_mags.py ${meta.id} ${mags}
    """
}