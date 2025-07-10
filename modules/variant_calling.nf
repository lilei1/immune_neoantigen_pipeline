/*
========================================================================================
    VARIANT CALLING MODULES
========================================================================================
*/

process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7':
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    path fasta

    output:
    tuple path(fasta), path("*.{amb,ann,bwt,pac,sa}"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bwa index $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0':
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0' }"

    input:
    tuple val(meta), path(reads)
    tuple path(fasta), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        -R "$read_group" \\
        \$INDEX \\
        ${reads[0]} \\
        ${reads[1]} \\
        | samtools sort -@ $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1':
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools index $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py39hdfd78af_0':
        'broadinstitute/gatk:4.3.0.0' }"

    input:
    tuple val(meta), path(tumor_bam), path(normal_bam)
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

    """
    gatk Mutect2 \\
        --input $tumor_bam \\
        --input $normal_bam \\
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
        'broadinstitute/gatk:4.3.0.0' }"

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
