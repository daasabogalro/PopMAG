process EXTRACT_CONTIG_NAMES {
    tag "${meta.id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/pip_biopython:326a6be8fb21b301'
    
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