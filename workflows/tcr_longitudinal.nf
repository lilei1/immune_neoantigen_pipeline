/*
========================================================================================
    TCR LONGITUDINAL ANALYSIS WORKFLOW
========================================================================================
*/

include { MIXCR_ALIGN                } from '../modules/tcr_analysis'
include { MIXCR_ASSEMBLE             } from '../modules/tcr_analysis'
include { MIXCR_FILTER               } from '../modules/tcr_analysis'
include { MIXCR_COLLAPSE             } from '../modules/tcr_analysis'
include { MIXCR_EXPORTCLONES         } from '../modules/tcr_analysis'
include { IMMUNARCH_ANALYSIS         } from '../modules/tcr_analysis'

workflow TCR_LONGITUDINAL {
    take:
    ch_reads    // channel: [ val(meta), [ reads ] ]
    
    main:
    ch_versions = Channel.empty()
    
    //
    // MODULE: MiXCR alignment
    //
    MIXCR_ALIGN (
        ch_reads
    )
    ch_versions = ch_versions.mix(MIXCR_ALIGN.out.versions.first())
    
    //
    // MODULE: MiXCR assembly
    //
    MIXCR_ASSEMBLE (
        MIXCR_ALIGN.out.vdjca
    )
    ch_versions = ch_versions.mix(MIXCR_ASSEMBLE.out.versions.first())
    
    //
    // MODULE: MiXCR filtering
    //
    MIXCR_FILTER (
        MIXCR_ASSEMBLE.out.clns
    )
    ch_versions = ch_versions.mix(MIXCR_FILTER.out.versions.first())

    //
    // MODULE: MiXCR collapsing
    //
    MIXCR_COLLAPSE (
        MIXCR_FILTER.out.clns_filtered
    )
    ch_versions = ch_versions.mix(MIXCR_COLLAPSE.out.versions.first())

    //
    // MODULE: MiXCR export clones
    //
    MIXCR_EXPORTCLONES (
        MIXCR_COLLAPSE.out.clns_collapsed
    )
    ch_versions = ch_versions.mix(MIXCR_EXPORTCLONES.out.versions.first())
    
    //
    // GROUP: Collect clonotype files by patient for longitudinal analysis
    //
    ch_clonotypes_grouped = MIXCR_EXPORTCLONES.out.clonotypes_tsv
        .map { meta, clonotypes ->
            def patient = meta.patient ?: meta.id.split('_')[0]
            return [ patient, meta, clonotypes ]
        }
        .groupTuple(by: 0)
        .map { patient, metas, clonotype_files ->
            def combined_meta = [
                id: patient,
                patient: patient,
                timepoints: metas.collect { it.timepoint ?: 'unknown' }.unique(),
                samples: metas.collect { it.id }
            ]
            return [ combined_meta, clonotype_files.flatten() ]
        }
    
    //
    // MODULE: Immunarch analysis
    //
    IMMUNARCH_ANALYSIS (
        ch_clonotypes_grouped
    )
    ch_versions = ch_versions.mix(IMMUNARCH_ANALYSIS.out.versions.first())
    
    emit:
    vdjca           = MIXCR_ALIGN.out.vdjca           // channel: [ val(meta), [ vdjca ] ]
    clns            = MIXCR_ASSEMBLE.out.clns         // channel: [ val(meta), [ clns ] ]
    clns_filtered   = MIXCR_FILTER.out.clns_filtered  // channel: [ val(meta), [ clns_filtered ] ]
    clns_collapsed  = MIXCR_COLLAPSE.out.clns_collapsed // channel: [ val(meta), [ clns_collapsed ] ]
    clonotypes      = MIXCR_EXPORTCLONES.out.clonotypes // channel: [ val(meta), [ clonotypes ] ]
    clonotypes_tsv  = MIXCR_EXPORTCLONES.out.clonotypes_tsv // channel: [ val(meta), [ clonotypes_tsv ] ]
    
    // Immunarch outputs
    diversity       = IMMUNARCH_ANALYSIS.out.diversity    // channel: [ val(meta), [ diversity_metrics ] ]
    expansion       = IMMUNARCH_ANALYSIS.out.expansion    // channel: [ val(meta), [ clonal_expansion ] ]
    tracking        = IMMUNARCH_ANALYSIS.out.tracking     // channel: [ val(meta), [ longitudinal_tracking ] ]
    plots           = IMMUNARCH_ANALYSIS.out.plots        // channel: [ val(meta), [ plots ] ]
    report          = IMMUNARCH_ANALYSIS.out.report       // channel: [ val(meta), [ report ] ]
    
    // Reports
    align_reports   = MIXCR_ALIGN.out.report         // channel: [ val(meta), [ report ] ]
    assemble_reports = MIXCR_ASSEMBLE.out.report     // channel: [ val(meta), [ report ] ]
    filter_reports  = MIXCR_FILTER.out.report        // channel: [ val(meta), [ report ] ]
    
    versions        = ch_versions                     // channel: [ versions.yml ]
}
