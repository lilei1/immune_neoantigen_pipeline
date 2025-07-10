#!/usr/bin/env python3
"""
Generate synthetic TCR-seq data for longitudinal analysis with 5 timepoints.
Creates realistic V(D)J sequences with clonal expansion patterns and temporal dynamics.
"""

import random
import gzip
import json
from pathlib import Path
from collections import defaultdict
import numpy as np

# Real TCR V, D, J gene segments (human)
TCR_V_GENES = [
    "TRBV2", "TRBV3-1", "TRBV4-1", "TRBV4-2", "TRBV4-3", "TRBV5-1", "TRBV5-4", "TRBV5-5", "TRBV5-6", "TRBV5-8",
    "TRBV6-1", "TRBV6-2", "TRBV6-3", "TRBV6-4", "TRBV6-5", "TRBV6-6", "TRBV6-8", "TRBV6-9", "TRBV7-2", "TRBV7-3",
    "TRBV7-4", "TRBV7-6", "TRBV7-7", "TRBV7-8", "TRBV7-9", "TRBV9", "TRBV10-1", "TRBV10-2", "TRBV10-3", "TRBV11-1",
    "TRBV11-2", "TRBV11-3", "TRBV12-3", "TRBV12-4", "TRBV12-5", "TRBV13", "TRBV14", "TRBV15", "TRBV16", "TRBV18",
    "TRBV19", "TRBV20-1", "TRBV24-1", "TRBV25-1", "TRBV27", "TRBV28", "TRBV29-1", "TRBV30"
]

TCR_D_GENES = [
    "TRBD1", "TRBD2"
]

TCR_J_GENES = [
    "TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6",
    "TRBJ2-1", "TRBJ2-2", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-6", "TRBJ2-7"
]

# CDR3 amino acid sequences (realistic examples)
CDR3_TEMPLATES = [
    "CASSQETQYF", "CASSLGQGAYEQYF", "CASSLDRGQPQHF", "CASSLGQGAYEQYF", "CASSQETQYF",
    "CASSLGQGAYEQYF", "CASSLDRGQPQHF", "CASSLGQGAYEQYF", "CASSQETQYF", "CASSLGQGAYEQYF",
    "CASSLDRGQPQHF", "CASSLGQGAYEQYF", "CASSQETQYF", "CASSLGQGAYEQYF", "CASSLDRGQPQHF",
    "CASSLGQGAYEQYF", "CASSQETQYF", "CASSLGQGAYEQYF", "CASSLDRGQPQHF", "CASSLGQGAYEQYF"
]

def generate_cdr3_sequence(template=None, length_range=(8, 20)):
    """Generate a realistic CDR3 amino acid sequence."""
    if template and random.random() < 0.3:  # 30% chance to use template
        return template
    
    # Generate random CDR3
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    length = random.randint(*length_range)
    
    # CDR3 typically starts with C and ends with F
    if length < 4:
        length = 8
    
    sequence = "C"  # Start with C
    for _ in range(length - 3):
        sequence += random.choice(amino_acids)
    sequence += "F"  # End with F
    
    return sequence

def generate_nucleotide_sequence(cdr3_aa, v_gene, j_gene):
    """Generate a realistic nucleotide sequence for the CDR3 region."""
    # Simple codon mapping (not exhaustive, but realistic)
    codon_map = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'], 'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 'K': ['AAA', 'AAG'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'], 'N': ['AAT', 'AAC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC']
    }
    
    # Convert amino acid sequence to nucleotides
    nucleotides = ""
    for aa in cdr3_aa:
        if aa in codon_map:
            nucleotides += random.choice(codon_map[aa])
        else:
            nucleotides += "NNN"  # Unknown amino acid
    
    # Add some flanking sequences (simplified)
    v_end = "TGTGCC"  # Typical V gene ending
    j_start = "TTCGGG"  # Typical J gene starting
    
    full_sequence = v_end + nucleotides + j_start
    return full_sequence

class TCRClone:
    """Represents a TCR clone with its properties and temporal dynamics."""
    
    def __init__(self, clone_id, v_gene, d_gene, j_gene, cdr3_aa, cdr3_nt):
        self.clone_id = clone_id
        self.v_gene = v_gene
        self.d_gene = d_gene
        self.j_gene = j_gene
        self.cdr3_aa = cdr3_aa
        self.cdr3_nt = cdr3_nt
        self.frequencies = {}  # timepoint -> frequency
        
    def set_temporal_pattern(self, pattern_type="stable"):
        """Set temporal dynamics for this clone."""
        base_freq = random.uniform(0.0001, 0.01)  # Base frequency
        
        if pattern_type == "expanding":
            # Clone expands over time
            multipliers = [1, 2, 5, 8, 12]
        elif pattern_type == "contracting":
            # Clone contracts over time
            multipliers = [10, 8, 5, 2, 1]
        elif pattern_type == "transient":
            # Clone appears and disappears
            multipliers = [1, 5, 10, 3, 0.5]
        else:  # stable
            # Clone remains relatively stable
            multipliers = [1, 1.2, 0.8, 1.1, 0.9]
        
        timepoints = ["T0_baseline", "T1_cycle1", "T2_cycle3", "T3_progression", "T4_posttreatment"]
        
        for i, tp in enumerate(timepoints):
            freq = base_freq * multipliers[i]
            # Add some noise
            freq *= random.uniform(0.8, 1.2)
            self.frequencies[tp] = max(freq, 0.00001)  # Minimum frequency

def generate_tcr_clones(num_clones=1000):
    """Generate a population of TCR clones with diverse properties."""
    clones = []
    
    # Define clone expansion patterns
    patterns = ["stable"] * 600 + ["expanding"] * 200 + ["contracting"] * 150 + ["transient"] * 50
    random.shuffle(patterns)
    
    for i in range(num_clones):
        clone_id = f"clone_{i+1:04d}"
        
        # Select V, D, J genes
        v_gene = random.choice(TCR_V_GENES)
        d_gene = random.choice(TCR_D_GENES)
        j_gene = random.choice(TCR_J_GENES)
        
        # Generate CDR3 sequence
        template = random.choice(CDR3_TEMPLATES) if random.random() < 0.1 else None
        cdr3_aa = generate_cdr3_sequence(template)
        cdr3_nt = generate_nucleotide_sequence(cdr3_aa, v_gene, j_gene)
        
        # Create clone
        clone = TCRClone(clone_id, v_gene, d_gene, j_gene, cdr3_aa, cdr3_nt)
        clone.set_temporal_pattern(patterns[i])
        
        clones.append(clone)
    
    return clones

def generate_reads_for_clone(clone, timepoint, total_reads=100000):
    """Generate sequencing reads for a clone at a specific timepoint."""
    frequency = clone.frequencies.get(timepoint, 0)
    num_reads = int(total_reads * frequency)
    
    if num_reads == 0:
        return []
    
    reads = []
    sequence = clone.cdr3_nt
    
    # Generate reads with some sequencing errors
    for i in range(num_reads):
        # Add sequencing errors (1% error rate)
        read_seq = ""
        for base in sequence:
            if random.random() < 0.01:  # 1% error rate
                read_seq += random.choice("ATGC")
            else:
                read_seq += base
        
        # Generate quality scores
        quality = ''.join(chr(33 + min(40, max(20, int(random.gauss(35, 5))))) for _ in range(len(read_seq)))
        
        reads.append((read_seq, quality, clone))
    
    return reads

def write_fastq_files(all_reads, sample_name, output_dir):
    """Write reads to paired FASTQ files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Shuffle reads to simulate random sequencing
    random.shuffle(all_reads)
    
    # Write R1 and R2 files
    with gzip.open(output_dir / f"{sample_name}_tcr_1.fastq.gz", 'wt') as f1, \
         gzip.open(output_dir / f"{sample_name}_tcr_2.fastq.gz", 'wt') as f2:
        
        for i, (seq, qual, clone) in enumerate(all_reads):
            # R1 (forward read)
            f1.write(f"@{clone.clone_id}_{i+1}/1\n{seq}\n+\n{qual}\n")
            
            # R2 (reverse read - reverse complement)
            rev_seq = seq[::-1].translate(str.maketrans("ATGC", "TACG"))
            f2.write(f"@{clone.clone_id}_{i+1}/2\n{rev_seq}\n+\n{qual}\n")

def generate_metadata_file(clones, output_dir):
    """Generate metadata file with clone information."""
    output_dir = Path(output_dir)
    
    metadata = {
        "clones": [],
        "timepoints": ["T0_baseline", "T1_cycle1", "T2_cycle3", "T3_progression", "T4_posttreatment"],
        "generation_info": {
            "total_clones": len(clones),
            "v_genes": len(set(c.v_gene for c in clones)),
            "j_genes": len(set(c.j_gene for c in clones)),
            "description": "Synthetic TCR data with longitudinal dynamics"
        }
    }
    
    for clone in clones:
        clone_data = {
            "clone_id": clone.clone_id,
            "v_gene": clone.v_gene,
            "d_gene": clone.d_gene,
            "j_gene": clone.j_gene,
            "cdr3_aa": clone.cdr3_aa,
            "cdr3_nt": clone.cdr3_nt,
            "frequencies": clone.frequencies
        }
        metadata["clones"].append(clone_data)
    
    with open(output_dir / "tcr_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

def generate_tcr_longitudinal_data(output_dir="test_data_tcr", num_clones=1000, reads_per_timepoint=50000):
    """Generate complete TCR longitudinal dataset."""
    print(f"Generating TCR longitudinal data in {output_dir}/")
    
    # Generate clone population
    print(f"Creating {num_clones} TCR clones with temporal dynamics...")
    clones = generate_tcr_clones(num_clones)
    
    # Generate data for each timepoint
    timepoints = ["T0_baseline", "T1_cycle1", "T2_cycle3", "T3_progression", "T4_posttreatment"]
    patients = ["PATIENT_01", "PATIENT_02"]
    
    for patient in patients:
        print(f"Generating data for {patient}...")
        
        for timepoint in timepoints:
            print(f"  Processing {timepoint}...")
            
            # Generate reads for all clones at this timepoint
            all_reads = []
            for clone in clones:
                reads = generate_reads_for_clone(clone, timepoint, reads_per_timepoint)
                all_reads.extend(reads)
            
            # Write FASTQ files
            sample_name = f"{patient}_{timepoint}"
            write_fastq_files(all_reads, sample_name, output_dir)
            
            print(f"    Generated {len(all_reads)} reads for {sample_name}")
    
    # Generate metadata
    print("Writing metadata...")
    generate_metadata_file(clones, output_dir)
    
    print(f"TCR longitudinal data generated successfully!")
    print(f"Total timepoints: {len(timepoints)}")
    print(f"Total patients: {len(patients)}")
    print(f"Total clones: {num_clones}")
    print(f"Reads per timepoint: ~{reads_per_timepoint}")

if __name__ == "__main__":
    generate_tcr_longitudinal_data()
