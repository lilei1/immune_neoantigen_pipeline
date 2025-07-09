/*
========================================================================================
    WES WORKFLOW
========================================================================================
    Whole Exome Sequencing workflow for variant calling and HLA typing
========================================================================================
*/

include { FASTQC                    } from '../modules/qc'
include { MUTECT2                   } from '../modules/variant_calling'
include { FILTERMUTECTCALLS         } from '../modules/variant_calling'
include { OPTITYPE                  } from '../modules/hla_typing'
include { MERGE_VARIANTS            } from '../modules/variant_calling'

workflow WES_WORKFLOW {
    
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
    // MODULE: Variant calling with Mutect2
    //
    // Separate tumor and normal samples for paired analysis
    ch_reads
        .branch { meta, reads ->
            tumor: meta.sample_type == 'tumor' || meta.sample_type == 'cfDNA'
            normal: meta.sample_type == 'normal'
        }
        .set { ch_samples_by_type }
    
    // Create tumor-normal pairs
    ch_tumor_normal_pairs = ch_samples_by_type.tumor
        .map { meta, reads -> 
            def patient_id = meta.patient_id
            def tumor_meta = meta.clone()
            tumor_meta.id = "${patient_id}_tumor"
            return [patient_id, tumor_meta, reads]
        }
        .combine(
            ch_samples_by_type.normal
                .map { meta, reads -> 
                    def patient_id = meta.patient_id
                    def normal_meta = meta.clone()
                    normal_meta.id = "${patient_id}_normal"
                    return [patient_id, normal_meta, reads]
                },
            by: 0
        )
        .map { patient_id, tumor_meta, tumor_reads, normal_meta, normal_reads ->
            def meta = [:]
            meta.id = patient_id
            meta.patient_id = patient_id
            meta.tumor_id = tumor_meta.id
            meta.normal_id = normal_meta.id
            return [meta, tumor_reads, normal_reads]
        }
    
    MUTECT2 (
        ch_tumor_normal_pairs
    )
    ch_versions = ch_versions.mix(MUTECT2.out.versions.first())
    
    //
    // MODULE: Filter Mutect2 calls
    //
    FILTERMUTECTCALLS (
        MUTECT2.out.vcf
    )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions.first())
    
    //
    // MODULE: HLA typing with OptiType
    //
    // Run HLA typing on all samples (tumor and normal)
    OPTITYPE (
        ch_reads
    )
    ch_versions = ch_versions.mix(OPTITYPE.out.versions.first())
    
    //
    // MODULE: Merge variants across timepoints for each patient
    //
    // Group variants by patient for longitudinal analysis
    ch_patient_variants = FILTERMUTECTCALLS.out.vcf
        .map { meta, vcf, tbi ->
            def patient_id = meta.patient_id
            return [patient_id, meta, vcf, tbi]
        }
        .groupTuple(by: 0)
        .map { patient_id, metas, vcfs, tbis ->
            def meta = [:]
            meta.id = patient_id
            meta.patient_id = patient_id
            meta.sample_count = metas.size()
            return [meta, vcfs, tbis]
        }
    
    MERGE_VARIANTS (
        ch_patient_variants
    )
    ch_versions = ch_versions.mix(MERGE_VARIANTS.out.versions.first())
    
    emit:
    variants     = FILTERMUTECTCALLS.out.vcf    // channel: [ meta, vcf, tbi ]
    merged_variants = MERGE_VARIANTS.out.vcf    // channel: [ meta, vcf, tbi ]
    hla_types    = OPTITYPE.out.result          // channel: [ meta, tsv ]
    versions     = ch_versions                  // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files           // channel: [ files ]
}
