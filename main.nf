#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CHECKM2_DATABASEDOWNLOAD } from './modules/local/checkm2_database'
include { CHECKM2_PREDICT } from './modules/nf-core/checkm2/predict/main'
include { TRANSFORM_CHECKM2_REPORT } from './modules/local/transform_checkm2_report'
include { DREP } from './modules/local/drep'
include { FILTER_BINS } from './modules/local/filter_bins'
include { BOWTIE2_INSTRAIN_BUILD } from './modules/local/bowtie2_instrain_build'

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

    ch_filtered_bins.view()

    //Prepare input for dRep
    drep_input = ch_filtered_bins
        .join(ch_transformed_report)

    DREP(drep_input)
    ch_concatenated_mags = DREP.out.concatenated_mags

    //InStrain setup

    ch_instrain_genes = CHECKM2_PREDICT.out.checkm2_output
    ch_instrain_genes.view()

    BOWTIE2_INSTRAIN_BUILD(ch_concatenated_mags)

//    INSTRAIN()
}

