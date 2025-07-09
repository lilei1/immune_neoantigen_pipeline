/*
========================================================================================
    VARIANT CALLING MODULES
========================================================================================
*/

process MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py39hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py39hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_reads), path(normal_reads)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor_command = tumor_reads.collect{"--input $it"}.join(' ')
    def normal_command = normal_reads.collect{"--input $it"}.join(' ')
    
    """
    gatk Mutect2 \\
        $tumor_command \\
        $normal_command \\
        --reference $fasta \\
        --output ${prefix}.vcf.gz \\
        --tumor-sample ${meta.tumor_id} \\
        --normal-sample ${meta.normal_id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process FILTERMUTECTCALLS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py39hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py39hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.filtered.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*.filteringStats.tsv")  , emit: stats
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    gatk FilterMutectCalls \\
        --variant $vcf \\
        --reference $fasta \\
        --output ${prefix}.filtered.vcf.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

process MERGE_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0':
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.merged.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.merged.vcf.gz.tbi") , emit: tbi
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_list = vcfs.collect{it}.join(' ')
    
    """
    # Create list of VCF files
    echo "${vcf_list.split(' ').join('\\n')}" > vcf_list.txt
    
    # Merge VCF files
    bcftools merge \\
        --file-list vcf_list.txt \\
        --output-type z \\
        --output ${prefix}.merged.vcf.gz \\
        $args
    
    # Index the merged VCF
    bcftools index --tbi ${prefix}.merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
