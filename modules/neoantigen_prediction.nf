/*
========================================================================================
    NEOANTIGEN PREDICTION MODULES
========================================================================================
*/

process GENERATE_PEPTIDES {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ensembl-vep=107.0 conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:107.0--pl5321h4a94de4_0' :
        'quay.io/biocontainers/ensembl-vep:107.0--pl5321h4a94de4_0' }"

    input:
    tuple val(meta), path(vcf), path(transcripts), path(hla_alleles)

    output:
    tuple val(meta), path("*.peptides.fasta"), emit: peptides
    tuple val(meta), path("*.variant_info.tsv"), emit: variant_info
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def peptide_lengths = params.peptide_lengths ?: "8,9,10,11"
    
    """
    #!/usr/bin/env python3
    
    import sys
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import vcf
    
    def generate_mutant_peptides(variant, transcript_seq, peptide_lengths):
        \"\"\"Generate mutant peptides from variant\"\"\"
        peptides = []
        
        # Simple implementation - in practice would use more sophisticated tools
        # like pVACtools or similar
        
        for length in peptide_lengths:
            # Generate peptides around the mutation site
            # This is a simplified version
            start_pos = max(0, variant.POS - length)
            end_pos = min(len(transcript_seq), variant.POS + length)
            
            if end_pos - start_pos >= length:
                peptide_seq = transcript_seq[start_pos:start_pos + length]
                peptides.append({
                    'peptide_id': f"{variant.CHROM}_{variant.POS}_{length}mer",
                    'sequence': str(peptide_seq),
                    'length': length,
                    'variant_id': f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"
                })
        
        return peptides
    
    # Parse peptide lengths
    lengths = [int(x.strip()) for x in "$peptide_lengths".split(',')]
    
    # Read VCF file
    vcf_reader = vcf.Reader(open('$vcf', 'r'))
    
    # Read transcript data (simplified - would normally use transcript sequences)
    # For demo purposes, create dummy sequences
    all_peptides = []
    variant_info = []
    
    for i, variant in enumerate(vcf_reader):
        # Create dummy transcript sequence for demonstration
        dummy_seq = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG" * 10
        
        # Generate peptides
        peptides = generate_mutant_peptides(variant, dummy_seq, lengths)
        all_peptides.extend(peptides)
        
        # Store variant info
        variant_info.append({
            'variant_id': f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}",
            'chromosome': variant.CHROM,
            'position': variant.POS,
            'ref_allele': variant.REF,
            'alt_allele': str(variant.ALT[0]),
            'gene': 'DEMO_GENE'  # Would extract from VCF annotations
        })
    
    # Write peptides to FASTA
    peptide_records = []
    for peptide in all_peptides:
        record = SeqRecord(
            Seq(peptide['sequence']),
            id=peptide['peptide_id'],
            description=f"length={peptide['length']} variant={peptide['variant_id']}"
        )
        peptide_records.append(record)
    
    with open('${prefix}.peptides.fasta', 'w') as f:
        SeqIO.write(peptide_records, f, 'fasta')
    
    # Write variant info
    variant_df = pd.DataFrame(variant_info)
    variant_df.to_csv('${prefix}.variant_info.tsv', sep='\\t', index=False)
    
    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write('    biopython: 1.79\\n')
        f.write('    pandas: 1.4.3\\n')
    """
}

process NETMHCPAN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::netmhcpan=4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/netmhcpan:4.1--hdfd78af_0' :
        'quay.io/biocontainers/netmhcpan:4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(peptides), path(hla_alleles)

    output:
    tuple val(meta), path("*.netmhcpan_predictions.txt"), emit: predictions
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    # Read HLA alleles
    HLA_ALLELES=\$(cat $hla_alleles | tr '\\n' ',' | sed 's/,\$//')
    
    # Run NetMHCpan prediction
    netMHCpan \\
        -f $peptides \\
        -a \$HLA_ALLELES \\
        -BA \\
        -xls \\
        -xlsfile ${prefix}.netmhcpan_predictions.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhcpan: \$(netMHCpan -h 2>&1 | grep -o 'NetMHCpan [0-9.]*' | head -1 | sed 's/NetMHCpan //')
    END_VERSIONS
    """
}

process FILTER_NEOANTIGENS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'python:3.9-slim' }"

    input:
    tuple val(meta), path(predictions)

    output:
    tuple val(meta), path("*.filtered_neoantigens.tsv"), emit: filtered
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def binding_threshold = params.binding_threshold ?: 500
    
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import sys
    
    # Read NetMHCpan predictions
    df = pd.read_csv('$predictions', sep='\\t', skiprows=1)
    
    # Filter by binding affinity threshold
    strong_binders = df[df['Rank'] <= 0.5]  # Strong binders (top 0.5%)
    weak_binders = df[(df['Rank'] > 0.5) & (df['Rank'] <= 2.0)]  # Weak binders (0.5-2%)
    
    # Combine and add binding category
    strong_binders['binding_category'] = 'strong'
    weak_binders['binding_category'] = 'weak'
    
    filtered_neoantigens = pd.concat([strong_binders, weak_binders])
    
    # Additional filtering by binding affinity if specified
    if $binding_threshold > 0:
        filtered_neoantigens = filtered_neoantigens[
            filtered_neoantigens['Aff(nM)'] <= $binding_threshold
        ]
    
    # Sort by binding affinity
    filtered_neoantigens = filtered_neoantigens.sort_values('Aff(nM)')
    
    # Save filtered results
    filtered_neoantigens.to_csv('${prefix}.filtered_neoantigens.tsv', sep='\\t', index=False)
    
    print(f"Filtered {len(filtered_neoantigens)} neoantigens from {len(df)} predictions")
    
    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}

process PRIORITIZE_NEOANTIGENS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'python:3.9-slim' }"

    input:
    tuple val(meta), path(filtered_neoantigens), path(transcripts)

    output:
    tuple val(meta), path("*.prioritized_neoantigens.tsv"), emit: prioritized
    tuple val(meta), path("*.neoantigen_summary.txt")     , emit: summary
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import sys
    
    # Read filtered neoantigens
    neoantigens = pd.read_csv('$filtered_neoantigens', sep='\\t')
    
    # Read transcript expression data
    # This is simplified - would normally match transcript IDs to genes
    transcripts = pd.read_csv('$transcripts', sep='\\t')
    
    # Add expression scores (simplified)
    # In practice, would match genes to transcripts properly
    neoantigens['expression_score'] = 1.0  # Placeholder
    
    # Calculate priority score
    # Combine binding affinity, expression, and other factors
    neoantigens['priority_score'] = (
        (1000 - neoantigens['Aff(nM)']) / 1000 * 0.4 +  # Binding affinity (40%)
        neoantigens['expression_score'] * 0.3 +          # Expression (30%)
        (2.0 - neoantigens['Rank']) / 2.0 * 0.3          # Rank score (30%)
    )
    
    # Sort by priority score
    prioritized = neoantigens.sort_values('priority_score', ascending=False)
    
    # Save prioritized results
    prioritized.to_csv('${prefix}.prioritized_neoantigens.tsv', sep='\\t', index=False)
    
    # Generate summary
    summary_stats = {
        'total_neoantigens': len(prioritized),
        'strong_binders': len(prioritized[prioritized['binding_category'] == 'strong']),
        'weak_binders': len(prioritized[prioritized['binding_category'] == 'weak']),
        'top_10_avg_affinity': prioritized.head(10)['Aff(nM)'].mean(),
        'top_candidate': prioritized.iloc[0]['Peptide'] if len(prioritized) > 0 else 'None'
    }
    
    with open('${prefix}.neoantigen_summary.txt', 'w') as f:
        f.write("Neoantigen Prediction Summary\\n")
        f.write("=" * 30 + "\\n")
        for key, value in summary_stats.items():
            f.write(f"{key}: {value}\\n")
    
    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: {sys.version.split()[0]}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """
}
