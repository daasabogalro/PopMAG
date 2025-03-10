process EXTRACT_CONTIG_NAMES {
    label 'process_single'
    
    input:
    tuple val(sample_id), val(mag_ids), path(mags)
    
    output:
    tuple val(sample_id), path("${sample_id}_combined_contigs.txt"), emit: scaffold_to_bin_combined
    
    script:
    """
    process_mags.py ${sample_id} ${mags}
    """
}