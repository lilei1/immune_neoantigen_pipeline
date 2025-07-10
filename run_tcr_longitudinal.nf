#!/usr/bin/env nextflow

/*
========================================================================================
    TCR LONGITUDINAL ANALYSIS PIPELINE
========================================================================================
    A standalone pipeline for comprehensive TCR longitudinal analysis
    
    Usage:
    nextflow run run_tcr_longitudinal.nf -profile docker,test_tcr_longitudinal
    
    Or with custom data:
    nextflow run run_tcr_longitudinal.nf -profile docker --input samplesheet.csv
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

log.info """\
         TCR LONGITUDINAL ANALYSIS PIPELINE
         ===================================
         input    : ${params.input}
         outdir   : ${params.outdir}
         """
         .stripIndent()

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { TCR_LONGITUDINAL } from './workflows/tcr_longitudinal'

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_TCRLONGITUDINAL {

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    Channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:',')
        .map { create_tcr_channels_from_samplesheet(it) }
        .set { ch_input }

    //
    // WORKFLOW: Run pipeline
    //
    TCR_LONGITUDINAL (
        ch_input
    )

    emit:
    multiqc_report = Channel.empty() // TODO: Add MultiQC support
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    NFCORE_TCRLONGITUDINAL ()
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_tcr_channels_from_samplesheet(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.patient      = row.patient
    meta.type         = row.type
    meta.timepoint    = row.timepoint ?: 'unknown'
    meta.single_end   = false

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []

    // Check if TCR files exist
    def tcr1 = file(row.tcr_1)
    def tcr2 = file(row.tcr_2)

    if (!tcr1.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> TCR FastQ file does not exist!\n${row.tcr_1}"
    }
    if (!tcr2.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> TCR FastQ file does not exist!\n${row.tcr_2}"
    }

    fastq_meta = [ meta, [ tcr1, tcr2 ] ]
    return fastq_meta
}

/*
========================================================================================
    THE END
========================================================================================
*/
