/*
========================================================================================
    HLA TYPING MODULES
========================================================================================
*/

process OPTITYPE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::optitype=1.3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/optitype:1.3.5--py39h5371cbf_3':
        'umccr/optitype:latest' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_result.tsv")  , emit: result
    tuple val(meta), path("*_coverage_plot.pdf"), optional: true, emit: pdf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sequencing_type = meta.single_end ? "--rna" : "--dna"
    
    if (meta.single_end) {
        """
        OptiTypePipeline.py \\
            --input $reads \\
            $sequencing_type \\
            --prefix $prefix \\
            --outdir . \\
            --verbose $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            optitype: \$(OptiTypePipeline.py --version 2>&1 | grep -o 'OptiType [0-9.]*' | cut -d' ' -f2 || echo "unknown")
        END_VERSIONS
        """
    } else {
        """
        OptiTypePipeline.py \\
            --input ${reads[0]} ${reads[1]} \\
            $sequencing_type \\
            --prefix $prefix \\
            --outdir . \\
            --verbose $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            optitype: \$(OptiTypePipeline.py --version 2>&1 | grep -o 'OptiType [0-9.]*' | cut -d' ' -f2 || echo "unknown")
        END_VERSIONS
        """
    }
}

process PARSE_HLA_TYPES {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'python:3.9-slim' }"

    input:
    tuple val(meta), path(optitype_result)

    output:
    tuple val(meta), path("*.hla_alleles.txt"), emit: alleles
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import sys
    
    # Read OptiType result
    df = pd.read_csv('$optitype_result', sep='\\t')
    
    # Extract HLA alleles
    alleles = []
    for col in ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']:
        if col in df.columns:
            allele = df[col].iloc[0]
            if pd.notna(allele):
                # Format as HLA-A*02:01 style
                locus = col[0]
                alleles.append(f"HLA-{locus}*{allele}")
    
    # Write to output file
    with open('${prefix}.hla_alleles.txt', 'w') as f:
        for allele in alleles:
            f.write(f"{allele}\\n")
    
    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}
