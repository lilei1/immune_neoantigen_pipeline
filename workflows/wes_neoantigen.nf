/*
========================================================================================
    WES TO NEOANTIGEN PREDICTION WORKFLOW
========================================================================================
*/

include { MUTECT2                    } from '../modules/variant_calling'
include { VEP_ANNOTATION             } from '../modules/variant_annotation'
include { FILTER_NONSYNONYMOUS       } from '../modules/variant_annotation'
include { CROSS_CHECK_EXPRESSION     } from '../modules/variant_annotation'
include { GENERATE_PROTEIN_FASTA     } from '../modules/neoantigen_prediction'
include { NETMHCPAN_PREDICT          } from '../modules/neoantigen_prediction'

workflow WES_NEOANTIGEN {
    take:
    ch_tumor_normal_bam     // channel: [ val(meta), [ tumor_bam, tumor_bai, normal_bam, normal_bai ] ]
    ch_expression           // channel: [ val(meta), path(expression_file) ]
    ch_hla_alleles          // channel: [ val(meta), path(hla_alleles) ]
    fasta                   // reference genome FASTA
    fasta_fai               // reference genome FASTA index
    dict                    // reference genome dictionary
    vep_cache               // VEP cache directory
    transcript_fasta        // transcript sequences FASTA
    gtf                     // gene annotation GTF
    
    main:
    ch_versions = Channel.empty()
    
    //
    // MODULE: Mutect2 variant calling
    //
    MUTECT2 (
        ch_tumor_normal_bam,
        fasta,
        fasta_fai,
        dict
    )
    ch_versions = ch_versions.mix(MUTECT2.out.versions.first())
    
    //
    // MODULE: VEP annotation
    //
    VEP_ANNOTATION (
        MUTECT2.out.vcf,
        vep_cache,
        fasta
    )
    ch_versions = ch_versions.mix(VEP_ANNOTATION.out.versions.first())
    
    //
    // MODULE: Filter for nonsynonymous variants
    //
    FILTER_NONSYNONYMOUS (
        VEP_ANNOTATION.out.vcf.join(VEP_ANNOTATION.out.tbi)
    )
    ch_versions = ch_versions.mix(FILTER_NONSYNONYMOUS.out.versions.first())
    
    //
    // MODULE: Cross-check with expression data
    //
    // Join variant data with expression data by patient ID
    ch_variants_with_expression = FILTER_NONSYNONYMOUS.out.vcf
        .join(FILTER_NONSYNONYMOUS.out.tbi)
        .map { meta, vcf, tbi ->
            def patient_id = meta.patient ?: meta.id.split('_')[0]
            return [ patient_id, meta, vcf, tbi ]
        }
        .combine(
            ch_expression.map { meta, expr ->
                def patient_id = meta.patient ?: meta.id.split('_')[0]
                return [ patient_id, meta, expr ]
            },
            by: 0
        )
        .map { patient_id, meta1, vcf, tbi, meta2, expr ->
            return [ meta1, vcf, tbi, meta2, expr ]
        }
    
    CROSS_CHECK_EXPRESSION (
        ch_variants_with_expression.map { meta1, vcf, tbi, meta2, expr ->
            return [ meta1, vcf, tbi ]
        },
        ch_variants_with_expression.map { meta1, vcf, tbi, meta2, expr ->
            return [ meta2, expr ]
        }
    )
    ch_versions = ch_versions.mix(CROSS_CHECK_EXPRESSION.out.versions.first())
    
    //
    // MODULE: Generate protein FASTA sequences
    //
    GENERATE_PROTEIN_FASTA (
        CROSS_CHECK_EXPRESSION.out.vcf.join(CROSS_CHECK_EXPRESSION.out.tbi),
        transcript_fasta,
        gtf
    )
    ch_versions = ch_versions.mix(GENERATE_PROTEIN_FASTA.out.versions.first())
    
    //
    // MODULE: NetMHCpan prediction
    //
    // Join mutant peptides with HLA alleles by patient ID
    ch_peptides_with_hla = GENERATE_PROTEIN_FASTA.out.mutant_fasta
        .map { meta, fasta ->
            def patient_id = meta.patient ?: meta.id.split('_')[0]
            return [ patient_id, meta, fasta ]
        }
        .combine(
            ch_hla_alleles.map { meta, hla ->
                def patient_id = meta.patient ?: meta.id.split('_')[0]
                return [ patient_id, hla ]
            },
            by: 0
        )
        .map { patient_id, meta, fasta, hla ->
            return [ meta, fasta, hla ]
        }
    
    NETMHCPAN_PREDICT (
        ch_peptides_with_hla.map { meta, fasta, hla -> [ meta, fasta ] },
        ch_peptides_with_hla.map { meta, fasta, hla -> hla }.first()
    )
    ch_versions = ch_versions.mix(NETMHCPAN_PREDICT.out.versions.first())
    
    emit:
    // Raw variant calling outputs
    raw_vcf             = MUTECT2.out.vcf                    // channel: [ val(meta), path(vcf) ]
    raw_stats           = MUTECT2.out.stats                  // channel: [ val(meta), path(stats) ]
    
    // Annotation outputs
    annotated_vcf       = VEP_ANNOTATION.out.vcf            // channel: [ val(meta), path(vcf) ]
    vep_report          = VEP_ANNOTATION.out.report         // channel: [ val(meta), path(report) ]
    
    // Filtering outputs
    nonsynonymous_vcf   = FILTER_NONSYNONYMOUS.out.vcf      // channel: [ val(meta), path(vcf) ]
    filtering_stats     = FILTER_NONSYNONYMOUS.out.stats    // channel: [ val(meta), path(stats) ]
    
    // Expression filtering outputs
    expressed_vcf       = CROSS_CHECK_EXPRESSION.out.vcf    // channel: [ val(meta), path(vcf) ]
    expression_stats    = CROSS_CHECK_EXPRESSION.out.stats  // channel: [ val(meta), path(stats) ]
    
    // Protein generation outputs
    mutant_proteins     = GENERATE_PROTEIN_FASTA.out.mutant_fasta    // channel: [ val(meta), path(fasta) ]
    wildtype_proteins   = GENERATE_PROTEIN_FASTA.out.wildtype_fasta  // channel: [ val(meta), path(fasta) ]
    mutation_info       = GENERATE_PROTEIN_FASTA.out.mutation_info   // channel: [ val(meta), path(tsv) ]
    
    // Neoantigen prediction outputs
    predictions         = NETMHCPAN_PREDICT.out.predictions  // channel: [ val(meta), path(tsv) ]
    binding_summary     = NETMHCPAN_PREDICT.out.summary      // channel: [ val(meta), path(tsv) ]
    
    versions            = ch_versions                         // channel: [ versions.yml ]
}
