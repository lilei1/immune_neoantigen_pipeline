/*
========================================================================================
    Utility functions for the immune neoantigen pipeline
========================================================================================
*/

import groovy.json.JsonSlurper

/*
 * Validate parameters
 */
def validateParameters() {
    
    // Check required parameters
    if (!params.input) {
        log.error "ERROR: Please provide input samplesheet with --input"
        System.exit(1)
    }
    
    // Check output directory
    if (!params.outdir) {
        log.error "ERROR: Please provide output directory with --outdir"
        System.exit(1)
    }
    
    // Check workflow flags
    if (!params.run_wes && !params.run_rnaseq && !params.run_tcr) {
        log.error "ERROR: At least one workflow must be enabled (run_wes, run_rnaseq, or run_tcr)"
        System.exit(1)
    }
    
    // Check neoantigen prediction requirements
    if (params.run_neoantigen && (!params.run_wes || !params.run_rnaseq)) {
        log.error "ERROR: Neoantigen prediction requires both WES and RNA-seq workflows to be enabled"
        System.exit(1)
    }
    
    // Validate peptide lengths
    if (params.peptide_lengths) {
        def lengths = params.peptide_lengths.split(',').collect { it.trim() as Integer }
        if (lengths.any { it < 8 || it > 15 }) {
            log.error "ERROR: Peptide lengths must be between 8 and 15 amino acids"
            System.exit(1)
        }
    }
    
    // Validate binding threshold
    if (params.binding_threshold && (params.binding_threshold < 0 || params.binding_threshold > 50000)) {
        log.error "ERROR: Binding threshold must be between 0 and 50000 nM"
        System.exit(1)
    }
}

/*
 * Print help message
 */
def paramsHelp() {
    log.info"""
    ============================================================================
     Immune Repertoire + Neoantigen Pipeline v${workflow.manifest.version}
    ============================================================================
    
    Usage:
        nextflow run main.nf --input samplesheet.csv --outdir results [options]
    
    Required arguments:
        --input                 Path to input samplesheet (CSV format)
        --outdir                Output directory for results
    
    Workflow control:
        --run_wes               Run WES workflow (variant calling, HLA typing) [default: true]
        --run_rnaseq            Run RNA-seq workflow (transcript quantification) [default: true]
        --run_tcr               Run TCR-seq workflow (clonotype extraction) [default: true]
        --run_neoantigen        Run neoantigen prediction workflow [default: true]
    
    Reference genome:
        --genome                Reference genome [default: GRCh38]
        --fasta                 Path to reference FASTA file
        --gtf                   Path to gene annotation GTF file
    
    Variant calling options:
        --mutect2_extra_args    Additional arguments for Mutect2 [default: '']
        --min_base_quality_score Minimum base quality score [default: 20]
    
    RNA-seq options:
        --salmon_index          Path to pre-built Salmon index
        --salmon_extra_args     Additional arguments for Salmon [default: '']
    
    TCR-seq options:
        --mixcr_species         Species for MiXCR analysis [default: 'hsa']
        --mixcr_extra_args      Additional arguments for MiXCR [default: '']
    
    Neoantigen prediction:
        --netmhcpan_version     NetMHCpan version [default: '4.1']
        --peptide_lengths       Peptide lengths to predict (comma-separated) [default: '8,9,10,11']
        --binding_threshold     Binding affinity threshold in nM [default: 500]
    
    HLA typing:
        --optitype_extra_args   Additional arguments for OptiType [default: '']
    
    Quality control:
        --skip_qc               Skip quality control steps [default: false]
        --skip_multiqc          Skip MultiQC report generation [default: false]
    
    Resource limits:
        --max_cpus              Maximum number of CPUs [default: 16]
        --max_memory            Maximum memory [default: 128.GB]
        --max_time              Maximum time [default: 240.h]
    
    Other options:
        --help                  Show this help message and exit
        --version               Show pipeline version and exit
    
    Profiles:
        -profile test           Run with test data
        -profile docker         Run with Docker containers
        -profile singularity    Run with Singularity containers
        -profile slurm          Run on SLURM cluster
        -profile aws            Run on AWS Batch
    
    Examples:
        # Run with test data
        nextflow run main.nf -profile test,docker
        
        # Run on SLURM cluster
        nextflow run main.nf -profile slurm --input samples.csv --outdir results
        
        # Run only WES and neoantigen prediction
        nextflow run main.nf --input samples.csv --outdir results --run_rnaseq false --run_tcr false
        
        # Run with custom parameters
        nextflow run main.nf --input samples.csv --outdir results --peptide_lengths "9,10" --binding_threshold 1000
    
    ============================================================================
    """.stripIndent()
}

/*
 * Print parameter summary
 */
def paramsSummaryLog(workflow, params) {
    def summary = [:]
    
    // Core parameters
    summary['Pipeline Name']     = workflow.manifest.name
    summary['Pipeline Version']  = workflow.manifest.version
    summary['Run Name']          = workflow.runName
    summary['Input']             = params.input
    summary['Output directory']  = params.outdir
    
    // Workflow control
    summary['Run WES']           = params.run_wes
    summary['Run RNA-seq']       = params.run_rnaseq
    summary['Run TCR-seq']       = params.run_tcr
    summary['Run Neoantigen']    = params.run_neoantigen
    
    // Reference genome
    summary['Genome']            = params.genome
    if (params.fasta)            summary['FASTA'] = params.fasta
    if (params.gtf)              summary['GTF'] = params.gtf
    
    // Tool-specific parameters
    if (params.run_wes) {
        summary['Mutect2 args']      = params.mutect2_extra_args ?: 'None'
        summary['Min base quality']  = params.min_base_quality_score
    }
    
    if (params.run_rnaseq) {
        summary['Salmon args']       = params.salmon_extra_args ?: 'None'
        if (params.salmon_index)     summary['Salmon index'] = params.salmon_index
    }
    
    if (params.run_tcr) {
        summary['MiXCR species']     = params.mixcr_species
        summary['MiXCR args']        = params.mixcr_extra_args ?: 'None'
    }
    
    if (params.run_neoantigen) {
        summary['NetMHCpan version'] = params.netmhcpan_version
        summary['Peptide lengths']  = params.peptide_lengths
        summary['Binding threshold'] = "${params.binding_threshold} nM"
    }
    
    // Resource limits
    summary['Max CPUs']          = params.max_cpus
    summary['Max Memory']        = params.max_memory
    summary['Max Time']          = params.max_time
    
    // Configuration
    summary['Config Profile']    = workflow.profile
    summary['Container Engine']  = workflow.containerEngine
    if (workflow.containerEngine) {
        summary['Container']     = workflow.container
    }
    summary['Work dir']          = workflow.workDir
    summary['Launch dir']        = workflow.launchDir
    
    return summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
}

/*
 * Check if a file exists
 */
def checkPathParamList = { it ->
    def list = []
    for (path in it) {
        if (path) {
            def resolved = file(path, checkIfExists: true)
            list.add(resolved)
        }
    }
    return list
}

/*
 * Check mandatory parameters
 */
def checkMandatoryParams = { params_list ->
    def missing_params = []
    for (param in params_list) {
        if (!params[param]) {
            missing_params.add(param)
        }
    }
    if (missing_params.size() > 0) {
        log.error "ERROR: Missing required parameters: ${missing_params.join(', ')}"
        System.exit(1)
    }
}
