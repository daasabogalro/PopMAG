process EXTRACT_CONTIG_NAMES {
    label 'process_single'
    
    input:
    tuple val(sample_id), val(mag_id), path(mag)
    
    output:
    tuple val(sample_id), path("${sample_id}_${mag.baseName}_contigs.txt"), emit: scaffold_to_bin
    
    script:
    """
    process_mags.py ${mag} ${sample_id}_${mag.baseName}_contigs.txt
    """
}