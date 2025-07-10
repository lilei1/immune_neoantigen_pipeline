/*
========================================================================================
    WES WORKFLOW
========================================================================================
    Whole Exome Sequencing workflow for variant calling and HLA typing
========================================================================================
*/

include { FASTQC                    } from '../modules/qc'
include { BWA_INDEX                 } from '../modules/variant_calling'
include { BWA_MEM                   } from '../modules/variant_calling'
include { SAMTOOLS_INDEX            } from '../modules/variant_calling'
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

    // Prepare reference genome files
    ch_fasta = params.fasta ? Channel.fromPath(params.fasta).collect() : Channel.empty()
    ch_fai   = params.fasta ? Channel.fromPath(params.fasta + '.fai').collect() : Channel.empty()
    ch_dict  = params.fasta ? Channel.fromPath(params.fasta.replaceAll(/\.fa(sta)?$/, '.dict')).collect() : Channel.empty()
    
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
    // MODULE: BWA Index
    //
    BWA_INDEX (
        ch_fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    //
    // MODULE: BWA-MEM Alignment
    //
    BWA_MEM (
        ch_reads,
        BWA_INDEX.out.index
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Index BAM files
    //
    SAMTOOLS_INDEX (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Variant calling with Mutect2
    //
    // Separate tumor and normal samples for paired analysis
    SAMTOOLS_INDEX.out.bam_bai
        .branch { meta, bam, bai ->
            tumor: meta.sample_type == 'tumor' || meta.sample_type == 'cfDNA'
            normal: meta.sample_type == 'normal'
        }
        .set { ch_samples_by_type }
    
    // Create tumor-normal pairs
    ch_tumor_normal_pairs = ch_samples_by_type.tumor
        .map { meta, bam, bai ->
            def patient_id = meta.patient_id
            def tumor_meta = meta.clone()
            tumor_meta.id = "${patient_id}_tumor"
            return [patient_id, tumor_meta, bam]
        }
        .combine(
            ch_samples_by_type.normal
                .map { meta, bam, bai ->
                    def patient_id = meta.patient_id
                    def normal_meta = meta.clone()
                    normal_meta.id = "${patient_id}_normal"
                    return [patient_id, normal_meta, bam]
                },
            by: 0
        )
        .map { patient_id, tumor_meta, tumor_bam, normal_meta, normal_bam ->
            def meta = [:]
            meta.id = patient_id
            meta.patient_id = patient_id
            meta.tumor_id = tumor_meta.id
            meta.normal_id = normal_meta.id
            return [meta, tumor_bam, normal_bam]
        }
    
    MUTECT2 (
        ch_tumor_normal_pairs,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(MUTECT2.out.versions.first())
    
    //
    // MODULE: Filter Mutect2 calls
    //
    // Combine VCF outputs with stats for filtering
    ch_mutect2_for_filtering = MUTECT2.out.vcf
        .join(MUTECT2.out.tbi, by: 0)
        .join(MUTECT2.out.stats, by: 0)

    FILTERMUTECTCALLS (
        ch_mutect2_for_filtering,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions.first())
    
    //
    // MODULE: HLA typing with OptiType (conditional)
    //
    // Run HLA typing on all samples (tumor and normal)
    if (!params.skip_hla_typing) {
        OPTITYPE (
            ch_reads
        )
        ch_versions = ch_versions.mix(OPTITYPE.out.versions.first())
        ch_hla_types = OPTITYPE.out.result
    } else {
        ch_hla_types = Channel.empty()
    }
    
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
    hla_types    = ch_hla_types                 // channel: [ meta, tsv ]
    versions     = ch_versions                  // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files           // channel: [ files ]
}
