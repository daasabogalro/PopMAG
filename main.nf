#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CHECKM2_DATABASEDOWNLOAD } from './modules/local/checkm2_database'
include { CHECKM2_PREDICT } from './modules/nf-core/checkm2/predict/main'
include { TRANSFORM_CHECKM2_REPORT } from './modules/local/transform_checkm2_report'
include { DREP } from './modules/local/drep'
include { FILTER_BINS } from './modules/local/filter_bins'
include { BOWTIE2_INSTRAIN_BUILD } from './modules/local/bowtie2_instrain_build'
include { BOWTIE2_INSTRAIN_ALIGN } from './modules/local/bowtie2_instrain_align'
include { EXTRACT_CONTIG_NAMES } from './modules/local/mag_to_contig'
include { PROKKA } from './modules/nf-core/prokka/main'
include { PRODIGAL } from './modules/nf-core/prodigal/main'
include { INSTRAIN_PROFILE } from './modules/local/instrain_profile'

// Define the main workflow
workflow {
    // Create a channel for input data     
    mag_ch = Channel
        .fromPath(params.mag_paths)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            tuple([id: row.sample_id], file(row.mag_path))
    }
    
    //Download CheckM2 database
    if (!params.skip_checkm2) {
        if (params.checkm2_db) {
            ch_checkm2_db = [[:], file(params.checkm2_db, checkIfExists: true)]
        }
        else {
            CHECKM2_DATABASEDOWNLOAD()
            ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
            }
    }
    
    // Calculate CheckM2 metrics for all MAGs
    CHECKM2_PREDICT(mag_ch.groupTuple(), ch_checkm2_db)
    ch_checkm2_report = CHECKM2_PREDICT.out.checkm2_tsv

    //Channel to get the extensions of the MAGs
    extension_ch = Channel
    .fromPath(params.mag_paths)
    .splitCsv(header:true, sep:'\t')
    .map { row -> file(row.mag_path).extension }
    .first()

    TRANSFORM_CHECKM2_REPORT(ch_checkm2_report, extension_ch)
    ch_transformed_report = TRANSFORM_CHECKM2_REPORT.out.transformed_report

    //Prepare input for filtering
    filter_input = mag_ch.groupTuple()
        .join(TRANSFORM_CHECKM2_REPORT.out.transformed_report)

    FILTER_BINS(filter_input)
    ch_filtered_bins = FILTER_BINS.out.filtered_bins

    //Prepare input for dRep
    drep_input = ch_filtered_bins
        .join(ch_transformed_report)

    DREP(drep_input)
    ch_concatenated_mags = DREP.out.concatenated_mags

    //InStrain setup

    BOWTIE2_INSTRAIN_BUILD(ch_concatenated_mags)
    ch_bowtie2_instrain_index = BOWTIE2_INSTRAIN_BUILD.out.index

    reads_ch = Channel
	.fromPath(params.reads_paths)
	.splitCsv(header:true, sep:'\t')
	.map { row -> tuple([id: row.sample_id], [file(row.forward), file(row.reverse)]) } 
    .unique()

// Cartesian product of MAGs+index with reads
ch_bowtie2_align_input = ch_bowtie2_instrain_index.combine(reads_ch)
    .map { mag_meta, mag, index, reads_meta, reads -> 
        tuple(mag_meta, mag, index, reads_meta, reads)
    }

    BOWTIE2_INSTRAIN_ALIGN(ch_bowtie2_align_input)
    ch_bowtie2_mapping = BOWTIE2_INSTRAIN_ALIGN.out.mappings

ch_complete_mags = Channel
    .fromPath(params.mag_paths)
    .splitCsv(header:true, sep:'\t')
    .map { row -> 
        def mag_id = file(row.mag_path).name.replaceFirst(/\.fa$/, '')
        tuple(row.sample_id, mag_id, file(row.mag_path)) 
    }

    EXTRACT_CONTIG_NAMES(ch_complete_mags.groupTuple())
    ch_contig_names = EXTRACT_CONTIG_NAMES.out.scaffold_to_bin_combined
    .map { meta, contigs ->
        tuple([id: meta], file(contigs))
    }

    // TODO: Find a way to simplify the usage of channels that manipulate params.mag_paths 
    ch_prokka_mags = Channel
    .fromPath(params.mag_paths)
    .splitCsv(header:true, sep:'\t')
    .map { row -> 
        def mag_id = file(row.mag_path).name.replaceFirst(/\.fa$/, '')
        tuple([id: mag_id, sample: row.sample_id], file(row.mag_path)) 
    }

//    PROKKA(ch_prokka_mags, [], [])
//    ch_prokka_fna = PROKKA.out.fna.map { meta, fna_files -> 
//        [[id: meta.sample], fna_files]
//    }
//    .groupTuple(by: 0)

    PRODIGAL(ch_prokka_mags,'gff')
    ch_prokka_fna = PRODIGAL.out.nucleotide_fasta

    ch_instrain_input = ch_bowtie2_mapping.combine(ch_prokka_fna, by: 0).combine(ch_contig_names, by: 0)

    //TODO: Check if every input file for InStrain is being produced correctly and, if so, why it's not running. 
    //INSTRAIN_PROFILE(ch_instrain_input)
}