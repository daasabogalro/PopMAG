#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CHECKM2_DATABASEDOWNLOAD } from './modules/local/checkm2_database'
//include { CHECKM2_DATABASEDOWNLOAD } from './modules/nf-core/checkm2/databasedownload/main'
include { CHECKM2_PREDICT } from './modules/nf-core/checkm2/predict/main'
include { TRANSFORM_CHECKM2_REPORT } from './modules/local/transform_checkm2_report'
include { FILTER_BINS } from './modules/local/filter_bins'
include { DREP } from './modules/local/drep'
include { COVERM } from './modules/local/coverm'
include { GENERATE_HEATMAP } from './modules/local/generate_heatmap'
include { BOWTIE2_INSTRAIN_BUILD } from './modules/local/bowtie2_instrain_build'
include { BOWTIE2_INSTRAIN_ALIGN } from './modules/local/bowtie2_instrain_align'
include { EXTRACT_CONTIG_NAMES } from './modules/local/mag_to_contig'
include { PRODIGAL } from './modules/nf-core/prodigal/main'
include { METACERBERUS_DATABASEDOWNLOAD } from './modules/local/metacerberus_database'
include { INSTRAIN_PROFILE } from './modules/local/instrain_profile'
include { INSTRAIN_COMPARE } from './modules/local/instrain_compare'
include { METACERBERUS_ANNOTATION } from './modules/local/metacerberus_annotation'
include { SNVS_TO_VCF } from './modules/local/snvs_to_vcf'
include { SUBSET_VCF_BY_GENOME } from './modules/local/subset_vcf_by_genome'
include { POGENOM } from './modules/local/pogenom'
include { EXTRACT_BIN_METRICS } from './modules/local/extract_bin_metrics'
include { MERGE_REPORTS_METACERBERUS_INSTRAIN } from './modules/local/merge_reports_metacerberus_instrain'
include { SINGLEM_METAPACKAGE } from './modules/local/singleM_metapackage'
include { SINGLEM } from './modules/local/singleM'

workflow {
    // MAGs channel from input file
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
            CHECKM2_DATABASEDOWNLOAD()//params.checkm2_db_version) This was needed for nf-core module
            ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
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
            ch_dereplicated_genomes = DREP.out.dereplicated_genomes

            EXTRACT_CONTIG_NAMES(ch_dereplicated_genomes)
            ch_combined_contig_names = EXTRACT_CONTIG_NAMES.out.contigs_to_bin_combined
            ch_individual_contig_names = EXTRACT_CONTIG_NAMES.out.contigs_to_bin_individual

    }

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
    ch_bams = BOWTIE2_INSTRAIN_ALIGN.out.bams.groupTuple()

    ch_coverm_input = ch_individual_contig_names
    .combine(ch_bams, by: 0
    )

    COVERM(ch_coverm_input)
    ch_coverage_file = COVERM.out.coverage

    GENERATE_HEATMAP(ch_coverage_file)

    // TODO: Find a way to simplify the usage of channels that manipulate params.mag_paths 
    ch_prodigal_mags = ch_dereplicated_genomes
    .map { meta, files ->
        // Ensure files is always a list
        def fileList = files instanceof List ? files : [files]
        tuple(meta, fileList)
    }
    .flatMap { meta, files ->
        files.collect { file ->
            def mag_id = file.name.replaceFirst(/\.fa$/, '')
            def sample_id = meta.id
            tuple([id: mag_id, sample: sample_id], file)
        }
    }

    PRODIGAL(ch_prodigal_mags,'gff')
    ch_prodigal_faa = PRODIGAL.out.amino_acid_fasta
    ch_prodigal_fna = PRODIGAL.out.nucleotide_fasta
    
    if (!params.skip_metacerberus) {
        if (params.metacerberus_db) {
            ch_metacerberus_db = [[:], file(params.metacerberus_db, checkIfExists: true)]
        }
        else {
            METACERBERUS_DATABASEDOWNLOAD()
            ch_metacerberus_db = METACERBERUS_DATABASEDOWNLOAD.out.database

            METACERBERUS_ANNOTATION(ch_prodigal_faa, ch_metacerberus_db)
            ch_metacerberus_annotations = METACERBERUS_ANNOTATION.out.final_annotation
        }
    }

    ch_instrain_input = ch_bowtie2_mapping
    .combine(ch_prodigal_fna
        .map { metadata, path ->
        tuple([id: metadata.sample], path)
        }
        .groupTuple(
        by: [0],
        ),by: 0)
    .combine(ch_combined_contig_names, by: 0)

    INSTRAIN_PROFILE(ch_instrain_input)
    ch_profiles = INSTRAIN_PROFILE.out.profile
    ch_snvs = INSTRAIN_PROFILE.out.snvs
    ch_gene_info = INSTRAIN_PROFILE.out.gene_info

    // Convert SNVs to VCF format
    ch_snv_vcf_input = ch_snvs
        .groupTuple(by: 0)
        .combine(ch_combined_contig_names, by: 0)

    SNVS_TO_VCF(ch_snv_vcf_input)
    ch_vcf_files = SNVS_TO_VCF.out.vcf_files

    ch_subset_vcf_input = ch_vcf_files.combine(ch_combined_contig_names, by: 0)

    SUBSET_VCF_BY_GENOME(ch_subset_vcf_input)
    ch_genome_vcfs = SUBSET_VCF_BY_GENOME.out.genome_vcfs

    //TODO: Check if we can access the mag.id from the vcf file. Since we are combining the channels
    //by sample we are getting all the sample - MAG combinations
    ch_pogenom_input = ch_dereplicated_genomes
    .map { meta, files ->
        // Ensure files is always a list
        def fileList = files instanceof List ? files : [files]
        tuple(meta, fileList)
    }
    .flatMap { meta, files ->
        files.collect { file ->
            def mag_id = file.name.replaceFirst(/\.fa$/, '')
            def sample_id = meta.id
            tuple([sample: sample_id], file)
        }
    }.combine( ch_genome_vcfs
    .map { meta, files ->

        def fileList = files instanceof List ? files : [files]
        tuple(meta, fileList)
    }
    .flatMap { meta, files ->
        files.collect { file ->
            def mag_id = file.name.replaceFirst(/\.fa$/, '')
            def sample_id = meta.id
            tuple([sample: sample_id], file)
        }}, by: 0)

    POGENOM(ch_pogenom_input)
    ch_pogenom_results = POGENOM.out.pogenom_results

    // Extract bin metrics from inStrain gene_info files
    ch_extract_bin_metrics_input = ch_gene_info
        .groupTuple(by: 0)  
        .combine(ch_combined_contig_names, by: 0) 

    EXTRACT_BIN_METRICS(ch_extract_bin_metrics_input)
    ch_bin_metrics = EXTRACT_BIN_METRICS.out.bin_metrics

    // Merge inStrain metrics with MetaCerberus annotations
    if (!params.skip_metacerberus) {

        ch_bin_metrics_for_merge = ch_bin_metrics
            .map { meta, files ->
                def fileList = files instanceof List ? files : [files]
                tuple(meta, fileList)
            }
           .flatMap { meta, files ->
            files.collect { file ->
                def sample_id = meta.id
                def bin_name = file.name.replaceFirst(/_metrics\.tsv$/, '')
                tuple([id: bin_name, sample: sample_id], file)
            }
        }
       

        ch_merge_input = ch_bin_metrics_for_merge.combine(ch_metacerberus_annotations, by: 0)

        MERGE_REPORTS_METACERBERUS_INSTRAIN(ch_merge_input)
        ch_merged_reports = MERGE_REPORTS_METACERBERUS_INSTRAIN.out.merged_reports
    }

    SINGLEM_METAPACKAGE()
    ch_singleM_metapackage = SINGLEM_METAPACKAGE.out.singleM_metapackage

    SINGLEM(ch_dereplicated_genomes, ch_singleM_metapackage)
    ch_singleM_results = SINGLEM.out.singleM_profile

    //TODO: Add a way to merge the singleM results with the other results

    INSTRAIN_COMPARE(ch_profiles.groupTuple())
}