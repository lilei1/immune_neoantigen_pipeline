/*
========================================================================================
    TCR WORKFLOW
========================================================================================
    T-cell receptor sequencing workflow for clonotype extraction and tracking
========================================================================================
*/

include { FASTQC                    } from '../modules/qc'
include { MIXCR_ANALYZE             } from '../modules/tcr_analysis'
include { MIXCR_EXPORTCLONES        } from '../modules/tcr_analysis'
include { TRACK_CLONES              } from '../modules/tcr_analysis'

workflow TCR_WORKFLOW {
    
    take:
    ch_reads // channel: [ meta, [ fastq1, fastq2 ] ]
    
    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // MODULE: Run FastQC
    //
    if (!params.skip_qc) {
        FASTQC (
            ch_reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }
    
    //
    // MODULE: Analyze TCR sequences with MiXCR
    //
    MIXCR_ANALYZE (
        ch_reads
    )
    ch_versions = ch_versions.mix(MIXCR_ANALYZE.out.versions.first())
    
    //
    // MODULE: Export clonotypes
    //
    MIXCR_EXPORTCLONES (
        MIXCR_ANALYZE.out.clns
    )
    ch_versions = ch_versions.mix(MIXCR_EXPORTCLONES.out.versions.first())
    
    //
    // MODULE: Track clonotypes across timepoints
    //
    // Group clonotypes by patient for longitudinal tracking
    ch_patient_clonotypes = MIXCR_EXPORTCLONES.out.clonotypes
        .map { meta, clonotypes ->
            def patient_id = meta.patient_id
            return [patient_id, meta, clonotypes]
        }
        .groupTuple(by: 0)
        .map { patient_id, metas, clonotype_files ->
            def meta = [:]
            meta.id = patient_id
            meta.patient_id = patient_id
            meta.sample_count = metas.size()
            meta.timepoints = metas.collect { it.timepoint }.unique().sort()
            meta.sample_types = metas.collect { it.sample_type }.unique()
            return [meta, clonotype_files, metas]
        }
    
    TRACK_CLONES (
        ch_patient_clonotypes
    )
    ch_versions = ch_versions.mix(TRACK_CLONES.out.versions.first())
    
    emit:
    clonotypes   = MIXCR_EXPORTCLONES.out.clonotypes // channel: [ meta, clonotypes.txt ]
    tracking     = TRACK_CLONES.out.tracking         // channel: [ meta, tracking_results ]
    reports      = TRACK_CLONES.out.plots            // channel: [ meta, plots ]
    raw_data     = MIXCR_ANALYZE.out.clns            // channel: [ meta, clns ]
    versions     = ch_versions                        // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files                 // channel: [ files ]
}
