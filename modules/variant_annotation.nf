/*
========================================================================================
    VARIANT ANNOTATION MODULES
========================================================================================
*/

process VEP_ANNOTATION {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::ensembl-vep=110.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:110.1--pl5321h2a3209d_0' :
        'ensemblorg/ensembl-vep:release_110.1' }"

    input:
    tuple val(meta), path(vcf)
    path(vep_cache)
    path(fasta)

    output:
    tuple val(meta), path("*.vep.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.vep.vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*.vep.txt")        , emit: txt
    tuple val(meta), path("*.vep_summary.html"), emit: report
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = vep_cache ? "--cache --dir_cache ${vep_cache}" : "--database"
    def fasta_command = fasta ? "--fasta ${fasta}" : ""
    
    """
    # Run VEP annotation
    vep \\
        --input_file ${vcf} \\
        --output_file ${prefix}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --stats_file ${prefix}.vep_summary.html \\
        --tab \\
        --output_file ${prefix}.vep.txt \\
        ${cache_command} \\
        ${fasta_command} \\
        --species homo_sapiens \\
        --assembly GRCh38 \\
        --offline \\
        --everything \\
        --canonical \\
        --protein \\
        --symbol \\
        --numbers \\
        --domains \\
        --regulatory \\
        --variant_class \\
        --sift b \\
        --polyphen b \\
        --humdiv \\
        --af \\
        --af_1kg \\
        --af_gnomad \\
        --max_af \\
        --pubmed \\
        --uniprot \\
        --mane \\
        --tsl \\
        --appris \\
        --ccds \\
        --hgvs \\
        --hgvsp \\
        --hgvsc \\
        --var_synonyms \\
        --gene_phenotype \\
        --fork ${task.cpus} \\
        $args

    # Index the output VCF
    tabix -p vcf ${prefix}.vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}

process FILTER_NONSYNONYMOUS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0' :
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.nonsynonymous.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.nonsynonymous.vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*.filtering_stats.txt")      , emit: stats
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Count total variants before filtering
    echo "=== Variant Filtering Statistics ===" > ${prefix}.filtering_stats.txt
    echo "Total variants before filtering: \$(bcftools view -H ${vcf} | wc -l)" >> ${prefix}.filtering_stats.txt

    # Filter for nonsynonymous variants
    # Include: missense_variant, stop_gained, stop_lost, start_lost, frameshift_variant, 
    #          inframe_insertion, inframe_deletion, protein_altering_variant
    bcftools view \\
        --include 'CSQ~"missense_variant" || CSQ~"stop_gained" || CSQ~"stop_lost" || CSQ~"start_lost" || CSQ~"frameshift_variant" || CSQ~"inframe_insertion" || CSQ~"inframe_deletion" || CSQ~"protein_altering_variant"' \\
        --output-type z \\
        --output ${prefix}.nonsynonymous.vcf.gz \\
        $args \\
        ${vcf}

    # Index the filtered VCF
    tabix -p vcf ${prefix}.nonsynonymous.vcf.gz

    # Count variants after filtering
    echo "Nonsynonymous variants after filtering: \$(bcftools view -H ${prefix}.nonsynonymous.vcf.gz | wc -l)" >> ${prefix}.filtering_stats.txt

    # Extract variant consequences for summary
    echo "" >> ${prefix}.filtering_stats.txt
    echo "=== Variant Consequences Summary ===" >> ${prefix}.filtering_stats.txt
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CSQ\\n' ${prefix}.nonsynonymous.vcf.gz | \\
        cut -f5 | tr '|' '\\t' | cut -f2 | sort | uniq -c | sort -nr >> ${prefix}.filtering_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

process CROSS_CHECK_EXPRESSION {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.5.3 bioconda::pysam=0.21.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'python:3.9-slim' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(expression_file)

    output:
    tuple val(meta), path("*.expressed.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.expressed.vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*.expression_filter.txt"), emit: stats
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tpm_threshold = params.expression_tpm_threshold ?: 0.2
    def count_threshold = params.expression_count_threshold ?: 10
    
    """
    #!/usr/bin/env python3

    import pandas as pd
    import pysam
    import gzip
    import sys
    import os

    def load_expression_data(expression_file):
        \"\"\"Load expression data from Salmon or similar output.\"\"\"
        try:
            # Try to read as Salmon quant.sf format
            df = pd.read_csv(expression_file, sep='\\t')
            
            # Check if it's Salmon format (has TPM column)
            if 'TPM' in df.columns and 'NumReads' in df.columns:
                # Salmon format
                df['gene_id'] = df['Name'].str.split('.').str[0]  # Remove version
                return df[['gene_id', 'TPM', 'NumReads']].copy()
            elif 'tpm' in df.columns and 'expected_count' in df.columns:
                # Alternative format
                df['gene_id'] = df.index
                return df[['gene_id', 'tpm', 'expected_count']].rename(
                    columns={'tpm': 'TPM', 'expected_count': 'NumReads'})
            else:
                print(f"Warning: Unknown expression file format for {expression_file}")
                return pd.DataFrame(columns=['gene_id', 'TPM', 'NumReads'])
        except Exception as e:
            print(f"Error reading expression file {expression_file}: {e}")
            return pd.DataFrame(columns=['gene_id', 'TPM', 'NumReads'])

    def extract_gene_from_csq(csq_field):
        \"\"\"Extract gene ID from VEP CSQ field.\"\"\"
        if not csq_field or csq_field == '.':
            return None
        
        # CSQ format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...
        try:
            csq_entries = csq_field.split(',')
            genes = set()
            for entry in csq_entries:
                fields = entry.split('|')
                if len(fields) > 4:
                    gene_id = fields[4]  # Gene field
                    if gene_id and gene_id != '':
                        genes.add(gene_id.split('.')[0])  # Remove version
            return list(genes)
        except:
            return None

    def main():
        # Load expression data
        print(f"Loading expression data from ${expression_file}")
        expr_df = load_expression_data("${expression_file}")
        
        if expr_df.empty:
            print("Warning: No expression data loaded, keeping all variants")
            # Copy input to output
            os.system(f"cp ${vcf} ${prefix}.expressed.vcf.gz")
            os.system(f"cp ${tbi} ${prefix}.expressed.vcf.gz.tbi")
            with open("${prefix}.expression_filter.txt", 'w') as f:
                f.write("No expression data available - no filtering applied\\n")
            return

        # Filter expression data
        tpm_threshold = ${tpm_threshold}
        count_threshold = ${count_threshold}
        
        expressed_genes = set(expr_df[
            (expr_df['TPM'] > tpm_threshold) | 
            (expr_df['NumReads'] >= count_threshold)
        ]['gene_id'].tolist())
        
        print(f"Found {len(expressed_genes)} expressed genes (TPM > {tpm_threshold} or counts >= {count_threshold})")
        
        # Process VCF
        vcf_in = pysam.VariantFile("${vcf}")
        vcf_out = pysam.VariantFile("${prefix}.expressed.vcf.gz", 'w', header=vcf_in.header)
        
        total_variants = 0
        expressed_variants = 0
        
        for record in vcf_in:
            total_variants += 1
            
            # Extract gene IDs from CSQ field
            csq_info = record.info.get('CSQ', '')
            genes = extract_gene_from_csq(csq_info)
            
            # Check if any gene is expressed
            if genes and any(gene in expressed_genes for gene in genes):
                vcf_out.write(record)
                expressed_variants += 1
            elif not genes:
                # If no gene annotation, keep the variant
                vcf_out.write(record)
                expressed_variants += 1
        
        vcf_in.close()
        vcf_out.close()
        
        # Index output VCF
        pysam.tabix_index("${prefix}.expressed.vcf.gz", preset="vcf")
        
        # Write statistics
        with open("${prefix}.expression_filter.txt", 'w') as f:
            f.write("=== Expression Filtering Statistics ===\\n")
            f.write(f"Total variants before expression filtering: {total_variants}\\n")
            f.write(f"Variants in expressed genes: {expressed_variants}\\n")
            f.write(f"Filtering rate: {(total_variants - expressed_variants) / total_variants * 100:.2f}%\\n")
            f.write(f"TPM threshold: {tpm_threshold}\\n")
            f.write(f"Count threshold: {count_threshold}\\n")
            f.write(f"Total expressed genes: {len(expressed_genes)}\\n")

    if __name__ == "__main__":
        main()

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
