/*
========================================================================================
    RNA-SEQ WORKFLOW
========================================================================================
    RNA sequencing workflow for transcript quantification
========================================================================================
*/

include { FASTQC                    } from '../modules/qc'
include { SALMON_INDEX              } from '../modules/rna_quantification'
include { SALMON_QUANT              } from '../modules/rna_quantification'
include { MERGE_TRANSCRIPTS         } from '../modules/rna_quantification'

workflow RNASEQ_WORKFLOW {
    
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
    // MODULE: Build Salmon index if not provided
    //
    ch_salmon_index = Channel.empty()
    if (params.salmon_index) {
        ch_salmon_index = Channel.fromPath(params.salmon_index)
    } else {
        // Build index from reference files
        if (params.fasta && params.gtf) {
            SALMON_INDEX (
                Channel.fromPath(params.fasta),
                Channel.fromPath(params.gtf)
            )
            ch_salmon_index = SALMON_INDEX.out.index
            ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())
        } else {
            error "ERROR: Either provide --salmon_index or both --fasta and --gtf for index building"
        }
    }
    
    //
    // MODULE: Quantify transcripts with Salmon
    //
    SALMON_QUANT (
        ch_reads,
        ch_salmon_index.collect()
    )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]}.ifEmpty([]))
    
    //
    // MODULE: Merge transcript quantifications across timepoints
    //
    // Group transcripts by patient for longitudinal analysis
    ch_patient_transcripts = SALMON_QUANT.out.results
        .map { meta, quant_dir ->
            def patient_id = meta.patient_id
            return [patient_id, meta, quant_dir]
        }
        .groupTuple(by: 0)
        .map { patient_id, metas, quant_dirs ->
            def meta = [:]
            meta.id = patient_id
            meta.patient_id = patient_id
            meta.sample_count = metas.size()
            meta.timepoints = metas.collect { it.timepoint }.unique()
            return [meta, quant_dirs]
        }
    
    MERGE_TRANSCRIPTS (
        ch_patient_transcripts
    )
    ch_versions = ch_versions.mix(MERGE_TRANSCRIPTS.out.versions.first())
    
    emit:
    transcripts = SALMON_QUANT.out.results      // channel: [ meta, quant_dir ]
    merged_transcripts = MERGE_TRANSCRIPTS.out.merged // channel: [ meta, merged_tsv ]
    index       = ch_salmon_index               // channel: [ index ]
    versions    = ch_versions                   // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files           // channel: [ files ]
}
