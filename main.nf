#!/usr/bin/env nextflow

/*
========================================================================================
    Immune Repertoire + Neoantigen Pipeline
========================================================================================
    A comprehensive pipeline for immune repertoire and neoantigen prediction
    in lung cancer patients using WES, RNA-seq, and TCR-seq data.

    Author: lilei
    Version: 1.0.0
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES AND WORKFLOWS
========================================================================================
*/

// Include utility functions
def validateParameters() { return true }
def paramsHelp() {
    log.info """
    Usage: nextflow run main.nf --input samplesheet.csv --outdir results
    """
}
def paramsSummaryLog(workflow, params) {
    return "Pipeline: ${workflow.manifest.name} v${workflow.manifest.version}"
}
include { WES_WORKFLOW }        from './workflows/wes_workflow'
include { RNASEQ_WORKFLOW }     from './workflows/rnaseq_workflow'
include { TCR_WORKFLOW }        from './workflows/tcr_workflow'
include { NEOANTIGEN_WORKFLOW } from './workflows/neoantigen_workflow'

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow IMMUNE_NEOANTIGEN_PIPELINE {

    take:
    input_samples // channel: [ meta, [files] ]

    main:

    // Initialize channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Separate samples by data type
    input_samples
        .branch { meta, files ->
            wes: meta.data_type == 'wes'
            rna: meta.data_type == 'rna'
            tcr: meta.data_type == 'tcr'
        }
        .set { ch_samples_by_type }

    // WES workflow - variant calling and HLA typing
    if (params.run_wes) {
        WES_WORKFLOW (
            ch_samples_by_type.wes
        )
        ch_variants = WES_WORKFLOW.out.variants
        ch_hla_types = WES_WORKFLOW.out.hla_types
        ch_versions = ch_versions.mix(WES_WORKFLOW.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(WES_WORKFLOW.out.multiqc_files)
    }

    // RNA-seq workflow - transcript quantification
    if (params.run_rnaseq) {
        RNASEQ_WORKFLOW (
            ch_samples_by_type.rna
        )
        ch_transcripts = RNASEQ_WORKFLOW.out.transcripts
        ch_versions = ch_versions.mix(RNASEQ_WORKFLOW.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(RNASEQ_WORKFLOW.out.multiqc_files)
    }

    // TCR workflow - clonotype extraction
    if (params.run_tcr) {
        TCR_WORKFLOW (
            ch_samples_by_type.tcr
        )
        ch_clonotypes = TCR_WORKFLOW.out.clonotypes
        ch_versions = ch_versions.mix(TCR_WORKFLOW.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TCR_WORKFLOW.out.multiqc_files)
    }

    // Neoantigen prediction workflow
    if (params.run_neoantigen && params.run_wes && params.run_rnaseq) {
        // Combine variants and transcripts by sample
        ch_variants
            .join(ch_transcripts, by: [0])
            .join(ch_hla_types, by: [0])
            .set { ch_neoantigen_input }

        NEOANTIGEN_WORKFLOW (
            ch_neoantigen_input
        )
        ch_neoantigens = NEOANTIGEN_WORKFLOW.out.prioritized
        ch_versions = ch_versions.mix(NEOANTIGEN_WORKFLOW.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(NEOANTIGEN_WORKFLOW.out.multiqc_files)
    }

    emit:
    variants     = params.run_wes ? ch_variants : Channel.empty()
    transcripts  = params.run_rnaseq ? ch_transcripts : Channel.empty()
    clonotypes   = params.run_tcr ? ch_clonotypes : Channel.empty()
    neoantigens  = (params.run_neoantigen && params.run_wes && params.run_rnaseq) ? ch_neoantigens : Channel.empty()
    hla_types    = params.run_wes ? ch_hla_types : Channel.empty()
    versions     = ch_versions
    multiqc_files = ch_multiqc_files
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Print help message if requested
    if (params.help) {
        paramsHelp()
        exit 0
    }

    // Validate parameters
    validateParameters()

    // Print parameter summary
    log.info paramsSummaryLog(workflow, params)

    // Create input channel from samplesheet
    if (params.input) {
        ch_input = Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .flatMap { row ->
                def meta = [:]
                meta.id = row.sample_id
                meta.patient_id = row.patient_id
                meta.timepoint = row.timepoint
                meta.sample_type = row.sample_type
                meta.hla_alleles = row.hla_alleles ?: null

                // Check which data types are available
                meta.has_wes = row.wes_r1 && row.wes_r2
                meta.has_rna = row.rna_r1 && row.rna_r2
                meta.has_tcr = row.tcr_r1 && row.tcr_r2

                // Create separate entries for each data type
                def entries = []

                if (meta.has_wes) {
                    def wes_meta = meta.clone()
                    wes_meta.data_type = 'wes'
                    entries.add([wes_meta, [file(row.wes_r1), file(row.wes_r2)]])
                }
                if (meta.has_rna) {
                    def rna_meta = meta.clone()
                    rna_meta.data_type = 'rna'
                    entries.add([rna_meta, [file(row.rna_r1), file(row.rna_r2)]])
                }
                if (meta.has_tcr) {
                    def tcr_meta = meta.clone()
                    tcr_meta.data_type = 'tcr'
                    entries.add([tcr_meta, [file(row.tcr_r1), file(row.tcr_r2)]])
                }

                return entries
            }
    } else {
        exit 1, "ERROR: Please provide input samplesheet with --input"
    }

    // Run main workflow
    IMMUNE_NEOANTIGEN_PIPELINE (
        ch_input
    )
}

/*
========================================================================================
    THE END
========================================================================================
*/