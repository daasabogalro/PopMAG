#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CHECKM2_DATABASEDOWNLOAD } from './modules/local/checkm2_database'
include { CHECKM2_PREDICT } from './modules/nf-core/checkm2/predict/main'
include { TRANSFORM_CHECKM2_REPORT } from './modules/local/transform_checkm2_report'

// Define the main workflow
workflow {
    // Create a channel for input data     
    mag_ch = Channel
        .fromPath(params.mag_paths)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            tuple([id: row.sample_id], file(row.mag_path))
    }
    
    //CheckM2
    if (!params.skip_checkm2) {
        if (params.checkm2_db) {
            ch_checkm2_db = [[:], file(params.checkm2_db, checkIfExists: true)]
        }
        else {
            CHECKM2_DATABASEDOWNLOAD()
            ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
            }
    }

    CHECKM2_PREDICT(mag_ch.groupTuple(), ch_checkm2_db)
    ch_checkm2_report = CHECKM2_PREDICT.out.checkm2_tsv

    //Channel to get the extensions of each MAG
    extension_ch = Channel
    .fromPath(params.mag_paths)
    .splitCsv(header:true, sep:'\t')
    .map { row -> file(row.mag_path).extension }
    .first()

    TRANSFORM_CHECKM2_REPORT(ch_checkm2_report, extension_ch)
}

