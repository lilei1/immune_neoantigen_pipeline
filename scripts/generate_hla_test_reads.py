#!/usr/bin/env python3
"""
Generate synthetic HLA reads for OptiType testing.
This creates minimal FASTQ files with HLA-like sequences that OptiType can process.
"""

import random
import gzip
from pathlib import Path

# Common HLA-A*02:01 exon 2 sequence (partial, for testing)
HLA_A_EXON2 = "GCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACG"

# Common HLA-B*07:02 exon 2 sequence (partial, for testing)  
HLA_B_EXON2 = "GCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGTCCGAGGATGGCGCCCCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACCGGGAGACACGGAACATGAAGGCCCAGTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCG"

# Common HLA-C*07:01 exon 2 sequence (partial, for testing)
HLA_C_EXON2 = "GCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACG"

def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def generate_reads_from_sequence(sequence, num_reads=100, read_length=150):
    """Generate synthetic reads from a reference sequence."""
    reads = []
    seq_len = len(sequence)
    
    for i in range(num_reads):
        # Random start position
        if seq_len > read_length:
            start = random.randint(0, seq_len - read_length)
            read_seq = sequence[start:start + read_length]
        else:
            # If sequence is shorter than read length, pad with random bases
            read_seq = sequence + ''.join(random.choices('ATGC', k=read_length - seq_len))
        
        # Add some random mutations (1% error rate)
        read_list = list(read_seq)
        for j in range(len(read_list)):
            if random.random() < 0.01:  # 1% mutation rate
                read_list[j] = random.choice('ATGC')
        
        read_seq = ''.join(read_list)
        
        # Generate quality scores (mostly high quality)
        quality = ''.join(chr(33 + min(40, max(20, int(random.gauss(35, 5))))) for _ in range(read_length))
        
        reads.append((read_seq, quality))
    
    return reads

def write_fastq_pair(reads1, reads2, prefix, output_dir):
    """Write paired FASTQ files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Write R1
    with gzip.open(output_dir / f"{prefix}_1.fastq.gz", 'wt') as f1:
        for i, (seq, qual) in enumerate(reads1):
            f1.write(f"@read_{i+1}/1\n{seq}\n+\n{qual}\n")
    
    # Write R2  
    with gzip.open(output_dir / f"{prefix}_2.fastq.gz", 'wt') as f2:
        for i, (seq, qual) in enumerate(reads2):
            f2.write(f"@read_{i+1}/2\n{seq}\n+\n{qual}\n")

def generate_hla_test_data(output_dir="test_data_hla", num_reads_per_gene=200):
    """Generate synthetic HLA test data."""
    print(f"Generating HLA test data in {output_dir}/")
    
    # Generate reads for each HLA gene
    hla_sequences = {
        'HLA_A': HLA_A_EXON2,
        'HLA_B': HLA_B_EXON2, 
        'HLA_C': HLA_C_EXON2
    }
    
    all_reads_r1 = []
    all_reads_r2 = []
    
    for gene_name, sequence in hla_sequences.items():
        print(f"Generating reads for {gene_name}...")
        
        # Generate forward reads
        reads_r1 = generate_reads_from_sequence(sequence, num_reads_per_gene, 150)
        
        # Generate reverse reads (reverse complement)
        reads_r2 = []
        for seq, qual in reads_r1:
            rev_seq = reverse_complement(seq)
            reads_r2.append((rev_seq, qual))
        
        all_reads_r1.extend(reads_r1)
        all_reads_r2.extend(reads_r2)
    
    # Add some random non-HLA reads to simulate real data
    print("Adding background reads...")
    random_seq = ''.join(random.choices('ATGC', k=1000))
    bg_reads_r1 = generate_reads_from_sequence(random_seq, 100, 150)
    bg_reads_r2 = [(reverse_complement(seq), qual) for seq, qual in bg_reads_r1]
    
    all_reads_r1.extend(bg_reads_r1)
    all_reads_r2.extend(bg_reads_r2)
    
    # Shuffle reads
    combined = list(zip(all_reads_r1, all_reads_r2))
    random.shuffle(combined)
    all_reads_r1, all_reads_r2 = zip(*combined)
    
    # Write FASTQ files for each sample
    samples = ['SAMPLE_01_T', 'SAMPLE_01_N', 'SAMPLE_02_T', 'SAMPLE_02_N']
    
    for sample in samples:
        print(f"Writing {sample} FASTQ files...")
        # Use subset of reads for each sample
        sample_size = len(all_reads_r1) // 2
        start_idx = hash(sample) % (len(all_reads_r1) - sample_size)
        
        sample_r1 = all_reads_r1[start_idx:start_idx + sample_size]
        sample_r2 = all_reads_r2[start_idx:start_idx + sample_size]
        
        write_fastq_pair(sample_r1, sample_r2, f"{sample}_wes", output_dir)
    
    print(f"HLA test data generated successfully in {output_dir}/")
    print(f"Total reads per sample: {sample_size}")
    print(f"HLA reads per gene: {num_reads_per_gene}")

if __name__ == "__main__":
    generate_hla_test_data()
