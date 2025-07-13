/*
========================================================================================
    NEOANTIGEN PREDICTION MODULES
========================================================================================
*/

process GENERATE_PROTEIN_FASTA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.9 bioconda::pysam=0.21.0 bioconda::biopython=1.81" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'python:3.9-slim' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(transcript_fasta)
    path(gtf)

    output:
    tuple val(meta), path("*.mutant_proteins.fasta")   , emit: mutant_fasta
    tuple val(meta), path("*.wildtype_proteins.fasta") , emit: wildtype_fasta
    tuple val(meta), path("*.mutation_info.tsv")       , emit: mutation_info
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3

    import pysam
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import re
    import sys

    def parse_vep_csq(csq_field):
        \"\"\"Parse VEP CSQ field to extract relevant information.\"\"\"
        if not csq_field or csq_field == '.':
            return []
        
        consequences = []
        for entry in csq_field.split(','):
            fields = entry.split('|')
            if len(fields) >= 15:  # Ensure we have enough fields
                consequence = {
                    'consequence': fields[1],
                    'gene': fields[4],
                    'transcript': fields[6],
                    'protein_position': fields[14] if len(fields) > 14 else '',
                    'amino_acids': fields[15] if len(fields) > 15 else '',
                    'hgvsp': fields[11] if len(fields) > 11 else ''
                }
                consequences.append(consequence)
        return consequences

    def extract_mutation_position(protein_pos_str):
        \"\"\"Extract numeric position from protein position string.\"\"\"
        if not protein_pos_str or protein_pos_str == '':
            return None
        
        # Handle formats like "123", "123-125", "123/456"
        match = re.search(r'(\\d+)', protein_pos_str)
        if match:
            return int(match.group(1))
        return None

    def generate_peptides_around_mutation(protein_seq, mut_pos, window=10):
        \"\"\"Generate peptides around mutation site.\"\"\"
        if mut_pos is None or mut_pos < 1:
            return []
        
        # Convert to 0-based indexing
        mut_pos_0 = mut_pos - 1
        
        peptides = []
        
        # MHC Class I peptides (8-11 mers)
        for length in range(8, 12):
            # Generate peptides centered on mutation
            for offset in range(-window, window + 1):
                start = mut_pos_0 + offset
                end = start + length
                
                if start >= 0 and end <= len(protein_seq):
                    peptide = str(protein_seq[start:end])
                    if len(peptide) == length and mut_pos_0 >= start and mut_pos_0 < end:
                        peptides.append({
                            'peptide': peptide,
                            'length': length,
                            'start': start + 1,  # Convert back to 1-based
                            'end': end,
                            'class': 'I',
                            'mut_position_in_peptide': mut_pos_0 - start + 1
                        })
        
        # MHC Class II peptides (13-25 mers)
        for length in range(13, 26):
            # Generate overlapping peptides
            for offset in range(-window, window + 1):
                start = mut_pos_0 + offset
                end = start + length
                
                if start >= 0 and end <= len(protein_seq):
                    peptide = str(protein_seq[start:end])
                    if len(peptide) == length and mut_pos_0 >= start and mut_pos_0 < end:
                        peptides.append({
                            'peptide': peptide,
                            'length': length,
                            'start': start + 1,  # Convert back to 1-based
                            'end': end,
                            'class': 'II',
                            'mut_position_in_peptide': mut_pos_0 - start + 1
                        })
        
        return peptides

    def apply_mutation(protein_seq, mut_pos, ref_aa, alt_aa):
        \"\"\"Apply amino acid mutation to protein sequence.\"\"\"
        if mut_pos is None or mut_pos < 1 or mut_pos > len(protein_seq):
            return protein_seq
        
        # Convert to 0-based indexing
        mut_pos_0 = mut_pos - 1
        
        # Check if reference amino acid matches
        if protein_seq[mut_pos_0] != ref_aa:
            print(f"Warning: Reference AA mismatch at position {mut_pos}: expected {ref_aa}, found {protein_seq[mut_pos_0]}")
        
        # Apply mutation
        mutant_seq = protein_seq[:mut_pos_0] + alt_aa + protein_seq[mut_pos_0 + 1:]
        return mutant_seq

    def main():
        # Load transcript sequences
        print("Loading transcript sequences...")
        transcripts = {}
        try:
            for record in SeqIO.parse("${transcript_fasta}", "fasta"):
                # Extract transcript ID (remove version if present)
                transcript_id = record.id.split('.')[0]
                transcripts[transcript_id] = record.seq
            print(f"Loaded {len(transcripts)} transcript sequences")
        except Exception as e:
            print(f"Error loading transcript FASTA: {e}")
            sys.exit(1)

        # Process VCF file
        print("Processing variants...")
        vcf_file = pysam.VariantFile("${vcf}")
        
        mutant_proteins = []
        wildtype_proteins = []
        mutation_info = []
        
        variant_count = 0
        processed_count = 0
        
        for record in vcf_file:
            variant_count += 1
            
            # Extract VEP consequences
            csq_field = record.info.get('CSQ', '')
            consequences = parse_vep_csq(csq_field)
            
            for csq in consequences:
                # Only process protein-coding consequences
                if 'missense_variant' not in csq['consequence']:
                    continue
                
                transcript_id = csq['transcript'].split('.')[0]
                if transcript_id not in transcripts:
                    continue
                
                # Extract mutation information
                mut_pos = extract_mutation_position(csq['protein_position'])
                if mut_pos is None:
                    continue
                
                # Parse amino acid change
                aa_change = csq['amino_acids']
                if '/' not in aa_change:
                    continue
                
                ref_aa, alt_aa = aa_change.split('/')
                
                # Get protein sequence (translate transcript)
                transcript_seq = transcripts[transcript_id]
                protein_seq = transcript_seq.translate()
                
                # Generate wildtype peptides
                wt_peptides = generate_peptides_around_mutation(protein_seq, mut_pos)
                
                # Apply mutation and generate mutant peptides
                mutant_protein_seq = apply_mutation(protein_seq, mut_pos, ref_aa, alt_aa)
                mut_peptides = generate_peptides_around_mutation(mutant_protein_seq, mut_pos)
                
                # Store results
                variant_id = f"{record.chrom}_{record.pos}_{record.ref}_{record.alts[0]}"
                
                for i, (wt_pep, mut_pep) in enumerate(zip(wt_peptides, mut_peptides)):
                    peptide_id = f"{variant_id}_{transcript_id}_{i}"
                    
                    # Wildtype protein
                    wt_record = SeqRecord(
                        Seq(wt_pep['peptide']),
                        id=f"WT_{peptide_id}",
                        description=f"Wildtype peptide from {transcript_id} position {wt_pep['start']}-{wt_pep['end']}"
                    )
                    wildtype_proteins.append(wt_record)
                    
                    # Mutant protein
                    mut_record = SeqRecord(
                        Seq(mut_pep['peptide']),
                        id=f"MUT_{peptide_id}",
                        description=f"Mutant peptide from {transcript_id} position {mut_pep['start']}-{mut_pep['end']}"
                    )
                    mutant_proteins.append(mut_record)
                    
                    # Mutation info
                    mutation_info.append({
                        'variant_id': variant_id,
                        'transcript_id': transcript_id,
                        'peptide_id': peptide_id,
                        'chromosome': record.chrom,
                        'position': record.pos,
                        'ref_allele': record.ref,
                        'alt_allele': str(record.alts[0]),
                        'gene': csq['gene'],
                        'protein_position': mut_pos,
                        'ref_aa': ref_aa,
                        'alt_aa': alt_aa,
                        'wildtype_peptide': wt_pep['peptide'],
                        'mutant_peptide': mut_pep['peptide'],
                        'peptide_length': wt_pep['length'],
                        'mhc_class': wt_pep['class'],
                        'peptide_start': wt_pep['start'],
                        'peptide_end': wt_pep['end'],
                        'mutation_position_in_peptide': wt_pep['mut_position_in_peptide']
                    })
                
                processed_count += 1

        print(f"Processed {processed_count} variants out of {variant_count} total variants")
        
        # Write FASTA files
        if wildtype_proteins:
            SeqIO.write(wildtype_proteins, "${prefix}.wildtype_proteins.fasta", "fasta")
            print(f"Wrote {len(wildtype_proteins)} wildtype peptides")
        else:
            # Create empty file
            open("${prefix}.wildtype_proteins.fasta", 'w').close()
        
        if mutant_proteins:
            SeqIO.write(mutant_proteins, "${prefix}.mutant_proteins.fasta", "fasta")
            print(f"Wrote {len(mutant_proteins)} mutant peptides")
        else:
            # Create empty file
            open("${prefix}.mutant_proteins.fasta", 'w').close()
        
        # Write mutation info
        if mutation_info:
            df = pd.DataFrame(mutation_info)
            df.to_csv("${prefix}.mutation_info.tsv", sep='\\t', index=False)
            print(f"Wrote mutation info for {len(mutation_info)} peptides")
        else:
            # Create empty file
            pd.DataFrame().to_csv("${prefix}.mutation_info.tsv", sep='\\t', index=False)

    if __name__ == "__main__":
        main()

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}

process NETMHCPAN_PREDICT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::netmhcpan=4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/netmhcpan:4.1--hdfd78af_1' :
        'quay.io/biocontainers/netmhcpan:4.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(peptides_fasta)
    path(hla_alleles)

    output:
    tuple val(meta), path("*.netmhcpan_predictions.tsv"), emit: predictions
    tuple val(meta), path("*.binding_summary.tsv")     , emit: summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ic50_threshold = params.netmhcpan_ic50_threshold ?: 500
    def percentile_threshold = params.netmhcpan_percentile_threshold ?: 2.0

    """
    # Read HLA alleles
    HLA_ALLELES=\$(cat ${hla_alleles} | tr '\\n' ',' | sed 's/,\$//')

    # Run NetMHCpan for each peptide length
    for LENGTH in 8 9 10 11; do
        # Extract peptides of specific length
        awk -v len=\$LENGTH '
        /^>/ { header=\$0; getline; if(length(\$0)==len) print header "\\n" \$0 }
        ' ${peptides_fasta} > peptides_\${LENGTH}mer.fasta

        if [ -s peptides_\${LENGTH}mer.fasta ]; then
            echo "Running NetMHCpan for \${LENGTH}-mer peptides..."
            netMHCpan -f peptides_\${LENGTH}mer.fasta \\
                -a \$HLA_ALLELES \\
                -l \$LENGTH \\
                -BA \\
                -xls \\
                -xlsfile ${prefix}_\${LENGTH}mer.xls \\
                > ${prefix}_\${LENGTH}mer.out
        fi
    done

    # Process NetMHCpan output
    python3 << 'EOF'
import pandas as pd
import glob
import os

def process_netmhcpan_output(xls_file):
    \"\"\"Process NetMHCpan XLS output file.\"\"\"
    try:
        # Read the XLS file (actually tab-separated)
        df = pd.read_csv(xls_file, sep='\\t', skiprows=1)

        # Standardize column names
        df.columns = df.columns.str.strip()

        # Extract relevant columns
        result_df = df[['Peptide', 'HLA', 'IC50(nM)', '%Rank']].copy()
        result_df.columns = ['peptide', 'hla_allele', 'ic50_nm', 'percentile_rank']

        # Add binding classification
        result_df['strong_binder'] = result_df['ic50_nm'] < 50
        result_df['weak_binder'] = (result_df['ic50_nm'] >= 50) & (result_df['ic50_nm'] < 500)
        result_df['binding_level'] = 'non_binder'
        result_df.loc[result_df['weak_binder'], 'binding_level'] = 'weak_binder'
        result_df.loc[result_df['strong_binder'], 'binding_level'] = 'strong_binder'

        # Add peptide length
        result_df['peptide_length'] = result_df['peptide'].str.len()

        return result_df
    except Exception as e:
        print(f"Error processing {xls_file}: {e}")
        return pd.DataFrame()

# Process all NetMHCpan output files
all_predictions = []
xls_files = glob.glob("${prefix}_*mer.xls")

for xls_file in xls_files:
    if os.path.exists(xls_file) and os.path.getsize(xls_file) > 0:
        df = process_netmhcpan_output(xls_file)
        if not df.empty:
            all_predictions.append(df)

if all_predictions:
    # Combine all predictions
    combined_df = pd.concat(all_predictions, ignore_index=True)

    # Save all predictions
    combined_df.to_csv("${prefix}.netmhcpan_predictions.tsv", sep='\\t', index=False)

    # Create binding summary
    ic50_threshold = ${ic50_threshold}
    percentile_threshold = ${percentile_threshold}

    # Filter for significant binders
    significant_binders = combined_df[
        (combined_df['ic50_nm'] < ic50_threshold) |
        (combined_df['percentile_rank'] < percentile_threshold)
    ].copy()

    # Add additional annotations
    significant_binders['meets_ic50_threshold'] = significant_binders['ic50_nm'] < ic50_threshold
    significant_binders['meets_percentile_threshold'] = significant_binders['percentile_rank'] < percentile_threshold

    # Sort by binding strength
    significant_binders = significant_binders.sort_values(['ic50_nm', 'percentile_rank'])

    # Save binding summary
    significant_binders.to_csv("${prefix}.binding_summary.tsv", sep='\\t', index=False)

    print(f"Total predictions: {len(combined_df)}")
    print(f"Significant binders (IC50 < {ic50_threshold} nM or %Rank < {percentile_threshold}%): {len(significant_binders)}")
    print(f"Strong binders (IC50 < 50 nM): {len(combined_df[combined_df['strong_binder']])}")
    print(f"Weak binders (50-500 nM): {len(combined_df[combined_df['weak_binder']])}")
else:
    # Create empty files
    pd.DataFrame().to_csv("${prefix}.netmhcpan_predictions.tsv", sep='\\t', index=False)
    pd.DataFrame().to_csv("${prefix}.binding_summary.tsv", sep='\\t', index=False)
    print("No NetMHCpan predictions generated")

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhcpan: \$(netMHCpan -h 2>&1 | grep "NetMHCpan version" | sed 's/.*version //; s/ .*//')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
