/*
========================================================================================
    RNA QUANTIFICATION MODULES
========================================================================================
*/

process SALMON_INDEX {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::salmon=1.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    input:
    path fasta
    path gtf

    output:
    path "salmon"       , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Extract transcript sequences from genome
    grep "^>" $fasta | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat $fasta > gentrome.fa
    
    # Build Salmon index
    salmon index \\
        --threads $task.cpus \\
        --transcripts gentrome.fa \\
        --decoys decoys.txt \\
        --index salmon \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}

process SALMON_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::salmon=1.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = meta.strandedness ?: 'A'
    
    if (meta.single_end) {
        """
        salmon quant \\
            --threads $task.cpus \\
            --libType=${strandedness} \\
            --index $index \\
            --reads $reads \\
            --output $prefix \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
    } else {
        """
        salmon quant \\
            --threads $task.cpus \\
            --libType=${strandedness} \\
            --index $index \\
            --mates1 ${reads[0]} \\
            --mates2 ${reads[1]} \\
            --output $prefix \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
    }
}

process MERGE_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'python:3.9-slim' }"

    input:
    tuple val(meta), path(quant_dirs)

    output:
    tuple val(meta), path("*.merged_transcripts.tsv"), emit: merged
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import sys
    
    # List of quantification directories
    quant_dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    
    merged_data = None
    
    for quant_dir in quant_dirs:
        quant_file = os.path.join(quant_dir, 'quant.sf')
        if os.path.exists(quant_file):
            df = pd.read_csv(quant_file, sep='\\t')
            df = df[['Name', 'TPM', 'NumReads']]
            df.columns = ['transcript_id', f'{quant_dir}_TPM', f'{quant_dir}_NumReads']
            
            if merged_data is None:
                merged_data = df
            else:
                merged_data = merged_data.merge(df, on='transcript_id', how='outer')
    
    # Fill NaN values with 0
    merged_data = merged_data.fillna(0)
    
    # Save merged data
    merged_data.to_csv('${prefix}.merged_transcripts.tsv', sep='\\t', index=False)
    
    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}
