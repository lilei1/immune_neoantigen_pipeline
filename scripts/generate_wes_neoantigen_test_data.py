#!/usr/bin/env python3
"""
Generate test data for WES-to-neoantigen workflow demonstration.
Creates synthetic VCF files, expression data, and HLA alleles for testing.
"""

import os
import pandas as pd
import random
from pathlib import Path

def create_synthetic_vcf_with_vep(output_file, num_variants=50):
    """Create a synthetic VCF file with VEP annotations."""
    
    # VCF header
    header = """##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP">
##VEP="v110" time="2024-01-01 12:00:00" cache="/path/to/cache" ensembl-version=110.1 ensembl-variation=110.1 ensembl-funcgen=110.1 ensembl-io=110.1 reference="GRCh38.fa" command="vep --everything --cache --offline"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
"""
    
    # Common genes and their transcript IDs
    genes = [
        ("TP53", "ENST00000269305", "ENSG00000141510"),
        ("KRAS", "ENST00000256078", "ENSG00000133703"),
        ("PIK3CA", "ENST00000263967", "ENSG00000121879"),
        ("EGFR", "ENST00000275493", "ENSG00000146648"),
        ("BRAF", "ENST00000288602", "ENSG00000157764"),
        ("APC", "ENST00000257430", "ENSG00000134982"),
        ("PTEN", "ENST00000371953", "ENSG00000171862"),
        ("MYC", "ENST00000377970", "ENSG00000136997"),
        ("RB1", "ENST00000267163", "ENSG00000139687"),
        ("BRCA1", "ENST00000357654", "ENSG00000012048")
    ]
    
    # Amino acid substitutions
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    
    variants = []
    
    for i in range(num_variants):
        # Random chromosome and position
        chrom = random.choice(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'])
        pos = random.randint(1000000, 200000000)
        
        # Random nucleotide change
        ref = random.choice(['A', 'T', 'G', 'C'])
        alt = random.choice([x for x in ['A', 'T', 'G', 'C'] if x != ref])
        
        # Select a gene
        gene_symbol, transcript_id, gene_id = random.choice(genes)
        
        # Generate amino acid change
        protein_pos = random.randint(1, 500)
        ref_aa = random.choice(amino_acids)
        alt_aa = random.choice([x for x in amino_acids if x != ref_aa])
        
        # Create VEP CSQ annotation
        csq_fields = [
            alt,  # Allele
            "missense_variant",  # Consequence
            "MODERATE",  # IMPACT
            gene_symbol,  # SYMBOL
            gene_id,  # Gene
            "Transcript",  # Feature_type
            transcript_id,  # Feature
            "protein_coding",  # BIOTYPE
            "5/10",  # EXON
            "",  # INTRON
            f"c.{protein_pos*3}G>A",  # HGVSc
            f"p.{ref_aa}{protein_pos}{alt_aa}",  # HGVSp
            str(protein_pos*3),  # cDNA_position
            str(protein_pos*3),  # CDS_position
            str(protein_pos),  # Protein_position
            f"{ref_aa}/{alt_aa}",  # Amino_acids
            f"{ref}gt/{alt}gt",  # Codons
            "",  # Existing_variation
            "",  # DISTANCE
            "1",  # STRAND
            "",  # FLAGS
            "SNV",  # VARIANT_CLASS
            "HGNC",  # SYMBOL_SOURCE
            "12345",  # HGNC_ID
            "YES",  # CANONICAL
        ]
        
        # Pad with empty fields to match VEP format
        while len(csq_fields) < 50:
            csq_fields.append("")
        
        csq_annotation = "|".join(csq_fields)
        
        # Create variant line
        variant_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\tCSQ={csq_annotation}\tGT\t0/1"
        variants.append(variant_line)
    
    # Write VCF file
    with open(output_file, 'w') as f:
        f.write(header)
        for variant in variants:
            f.write(variant + "\n")
    
    print(f"Created synthetic VCF with {num_variants} variants: {output_file}")

def create_synthetic_expression_data(output_file, num_genes=1000):
    """Create synthetic RNA-seq expression data (Salmon-like format)."""
    
    # Gene list with some overlap with VCF genes
    genes = [
        "ENST00000269305", "ENST00000256078", "ENST00000263967", "ENST00000275493", "ENST00000288602",
        "ENST00000257430", "ENST00000371953", "ENST00000377970", "ENST00000267163", "ENST00000357654"
    ]
    
    # Add random genes
    for i in range(num_genes - len(genes)):
        genes.append(f"ENST{random.randint(10000000, 99999999)}")
    
    # Generate expression data
    expression_data = []
    for gene in genes:
        # Some genes highly expressed, some lowly expressed
        if random.random() < 0.1:  # 10% highly expressed
            tpm = random.uniform(10, 1000)
            num_reads = random.uniform(100, 10000)
        elif random.random() < 0.3:  # 30% moderately expressed
            tpm = random.uniform(1, 10)
            num_reads = random.uniform(10, 100)
        else:  # 60% lowly expressed or not expressed
            tpm = random.uniform(0, 1)
            num_reads = random.uniform(0, 10)
        
        expression_data.append({
            'Name': gene,
            'Length': random.randint(500, 5000),
            'EffectiveLength': random.randint(400, 4500),
            'TPM': tpm,
            'NumReads': num_reads
        })
    
    # Create DataFrame and save
    df = pd.DataFrame(expression_data)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Created synthetic expression data with {num_genes} genes: {output_file}")

def create_hla_alleles_file(output_file):
    """Create HLA alleles file for NetMHCpan."""
    
    # Common HLA alleles
    hla_alleles = [
        "HLA-A*02:01",
        "HLA-A*01:01", 
        "HLA-A*03:01",
        "HLA-B*07:02",
        "HLA-B*08:01",
        "HLA-B*44:02",
        "HLA-C*07:01",
        "HLA-C*07:02",
        "HLA-C*16:01"
    ]
    
    with open(output_file, 'w') as f:
        for allele in hla_alleles:
            f.write(allele + "\n")
    
    print(f"Created HLA alleles file: {output_file}")

def create_transcript_fasta(output_file):
    """Create a simple transcript FASTA file."""
    
    # Sample transcript sequences (simplified)
    transcripts = {
        "ENST00000269305": "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCT",
        "ENST00000256078": "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCTCTTGGATATTCTCGACACAGCAGGTCAAGAGGAGTACAGTGCAATGAGGGACCAGTACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAAATGTGATTTGCCTTCTAGAACAGTAGACACAAAACAGGCTCAGGACTTAGCAAGAAGTTATGGAATTCCTTTTATTGAAACATCAGCAAAGACAAGACAGGGTGTTGATGATGGCGTGGTGGTGGGCGCTGTGGTGAAGGACATGCACATGCCCCGGTGCGAGCTGCCCAACACAAGCTGCATGACCCACTACCGCTGGGAGCTGGCCAAGAACCGCAAGAAGGCCATCGAGGGCTACAGTGACATGACCCTGGGCAGCGGCAAGGTGGTGACCTACGCCAAGACCGTGCTGCTGGTGGGCAACAAGAACGCCATGAACGCCTACGTGATGGCCCTGGGCAACATGCTGGACATGCAGAACGAGGTGAAGTTCCTGCGCTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGCCCCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACCGGGAGACACGGAACATGAAGGCCCAGTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCG"
    }
    
    with open(output_file, 'w') as f:
        for transcript_id, sequence in transcripts.items():
            f.write(f">{transcript_id}\n{sequence}\n")
    
    print(f"Created transcript FASTA file: {output_file}")

def create_simple_gtf(output_file):
    """Create a simple GTF file."""
    
    gtf_content = """#!genome-build GRCh38
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession NCBI:GCA_000001405.15
#!genebuild-last-updated 2013-09
1	ensembl	gene	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1	ensembl	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript";
17	ensembl	gene	7565097	7590856	.	-	.	gene_id "ENSG00000141510"; gene_version "17"; gene_name "TP53"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
17	ensembl	transcript	7565097	7590856	.	-	.	gene_id "ENSG00000141510"; gene_version "17"; transcript_id "ENST00000269305"; transcript_version "8"; gene_name "TP53"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TP53-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding";
"""
    
    with open(output_file, 'w') as f:
        f.write(gtf_content)
    
    print(f"Created GTF file: {output_file}")

def main():
    """Generate all test data for WES-to-neoantigen workflow."""
    
    print("🧬 GENERATING WES-TO-NEOANTIGEN TEST DATA")
    print("=" * 50)
    
    # Create output directory
    output_dir = Path("test_data_wes_neoantigen")
    output_dir.mkdir(exist_ok=True)
    
    # Generate test data for two patients
    patients = ["PATIENT_01", "PATIENT_02"]
    
    for patient in patients:
        print(f"\nGenerating data for {patient}...")
        
        # Create VCF with VEP annotations
        vcf_file = output_dir / f"{patient}_T.annotated.vcf"
        create_synthetic_vcf_with_vep(vcf_file, num_variants=30)
        
        # Create expression data
        expr_file = output_dir / f"{patient}_T.quant.sf"
        create_synthetic_expression_data(expr_file, num_genes=500)
        
        # Create HLA alleles
        hla_file = output_dir / f"{patient}.hla_alleles.txt"
        create_hla_alleles_file(hla_file)
    
    # Create reference files (shared)
    print("\nGenerating reference files...")
    
    # Transcript FASTA
    transcript_fasta = output_dir / "transcripts.fasta"
    create_transcript_fasta(transcript_fasta)
    
    # GTF file
    gtf_file = output_dir / "genes.gtf"
    create_simple_gtf(gtf_file)
    
    print(f"\n✅ Test data generation complete!")
    print(f"📁 Output directory: {output_dir}")
    print(f"📊 Files generated:")
    for file in sorted(output_dir.glob("*")):
        size = file.stat().st_size
        if size > 1024:
            size_str = f"{size/1024:.1f}KB"
        else:
            size_str = f"{size}B"
        print(f"   {file.name} ({size_str})")
    
    print(f"\n🔬 To test the workflow:")
    print(f"   1. Use the generated VCF files as input to VEP annotation")
    print(f"   2. Use the expression files for cross-checking")
    print(f"   3. Use the HLA alleles for NetMHCpan prediction")

if __name__ == "__main__":
    main()
