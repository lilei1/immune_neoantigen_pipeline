/*
========================================================================================
    NEOANTIGEN WORKFLOW
========================================================================================
    Neoantigen prediction workflow combining variants, transcripts, and HLA types
========================================================================================
*/

include { GENERATE_PEPTIDES         } from '../modules/neoantigen_prediction'
include { NETMHCPAN                 } from '../modules/neoantigen_prediction'
include { FILTER_NEOANTIGENS        } from '../modules/neoantigen_prediction'
include { PRIORITIZE_NEOANTIGENS    } from '../modules/neoantigen_prediction'

workflow NEOANTIGEN_WORKFLOW {
    
    take:
    ch_input // channel: [ meta, variants, transcripts, hla_types ]
    
    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // MODULE: Generate peptide sequences from variants
    //
    GENERATE_PEPTIDES (
        ch_input
    )
    ch_versions = ch_versions.mix(GENERATE_PEPTIDES.out.versions.first())
    
    //
    // MODULE: Predict MHC binding with NetMHCpan
    //
    NETMHCPAN (
        GENERATE_PEPTIDES.out.peptides
    )
    ch_versions = ch_versions.mix(NETMHCPAN.out.versions.first())
    
    //
    // MODULE: Filter neoantigens by binding affinity
    //
    FILTER_NEOANTIGENS (
        NETMHCPAN.out.predictions
    )
    ch_versions = ch_versions.mix(FILTER_NEOANTIGENS.out.versions.first())
    
    //
    // MODULE: Prioritize neoantigens using expression and other criteria
    //
    // Combine filtered neoantigens with transcript expression data
    ch_prioritization_input = FILTER_NEOANTIGENS.out.filtered
        .join(
            ch_input.map { meta, variants, transcripts, hla_types ->
                [meta, transcripts]
            },
            by: 0
        )
    
    PRIORITIZE_NEOANTIGENS (
        ch_prioritization_input
    )
    ch_versions = ch_versions.mix(PRIORITIZE_NEOANTIGENS.out.versions.first())
    
    emit:
    peptides     = GENERATE_PEPTIDES.out.peptides       // channel: [ meta, peptides.fasta ]
    predictions  = NETMHCPAN.out.predictions            // channel: [ meta, predictions.txt ]
    filtered     = FILTER_NEOANTIGENS.out.filtered      // channel: [ meta, filtered.txt ]
    prioritized  = PRIORITIZE_NEOANTIGENS.out.prioritized // channel: [ meta, prioritized.txt ]
    summary      = PRIORITIZE_NEOANTIGENS.out.summary   // channel: [ meta, summary.txt ]
    versions     = ch_versions                           // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files                    // channel: [ files ]
}
