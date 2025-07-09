#!/usr/bin/env python3

"""
Merge variant files across timepoints for longitudinal analysis
Author: lilei
"""

import argparse
import pandas as pd
import vcf
import sys
from pathlib import Path

def parse_vcf_to_dataframe(vcf_file):
    """Parse VCF file and convert to DataFrame"""
    variants = []
    
    with open(vcf_file, 'r') as f:
        vcf_reader = vcf.Reader(f)
        
        for record in vcf_reader:
            variant_data = {
                'chromosome': record.CHROM,
                'position': record.POS,
                'ref_allele': record.REF,
                'alt_allele': str(record.ALT[0]) if record.ALT else '',
                'quality': record.QUAL,
                'filter': ';'.join(record.FILTER) if record.FILTER else 'PASS',
                'variant_id': f"{record.CHROM}_{record.POS}_{record.REF}_{record.ALT[0] if record.ALT else ''}"
            }
            
            # Extract sample-specific information
            for sample in record.samples:
                sample_data = variant_data.copy()
                sample_data['sample'] = sample.sample
                sample_data['genotype'] = sample['GT'] if 'GT' in sample.data._fields else './.'
                sample_data['depth'] = sample['DP'] if 'DP' in sample.data._fields else 0
                sample_data['allele_freq'] = sample['AF'] if 'AF' in sample.data._fields else 0.0
                
                variants.append(sample_data)
    
    return pd.DataFrame(variants)

def merge_variant_dataframes(variant_dfs, sample_names):
    """Merge variant DataFrames from multiple timepoints"""
    
    if not variant_dfs:
        return pd.DataFrame()
    
    # Add sample information
    for i, df in enumerate(variant_dfs):
        df['timepoint'] = sample_names[i]
    
    # Combine all variants
    combined_df = pd.concat(variant_dfs, ignore_index=True)
    
    # Create summary statistics
    variant_summary = combined_df.groupby('variant_id').agg({
        'timepoint': lambda x: ';'.join(sorted(set(x))),
        'sample': lambda x: ';'.join(sorted(set(x))),
        'quality': 'mean',
        'depth': 'mean',
        'allele_freq': 'mean',
        'chromosome': 'first',
        'position': 'first',
        'ref_allele': 'first',
        'alt_allele': 'first'
    }).reset_index()
    
    # Add persistence information
    variant_summary['timepoint_count'] = variant_summary['timepoint'].str.split(';').str.len()
    variant_summary['is_persistent'] = variant_summary['timepoint_count'] > 1
    
    return variant_summary

def main():
    parser = argparse.ArgumentParser(description='Merge variant files across timepoints')
    parser.add_argument('--input', '-i', nargs='+', required=True,
                       help='Input VCF files')
    parser.add_argument('--output', '-o', required=True,
                       help='Output merged variants file')
    parser.add_argument('--sample-names', '-s', nargs='+',
                       help='Sample names corresponding to input files')
    parser.add_argument('--format', '-f', choices=['tsv', 'csv'], default='tsv',
                       help='Output format (default: tsv)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if len(args.input) == 0:
        print("Error: No input files provided", file=sys.stderr)
        sys.exit(1)
    
    # Use file basenames as sample names if not provided
    if not args.sample_names:
        args.sample_names = [Path(f).stem for f in args.input]
    
    if len(args.input) != len(args.sample_names):
        print("Error: Number of input files must match number of sample names", file=sys.stderr)
        sys.exit(1)
    
    # Parse VCF files
    print(f"Processing {len(args.input)} VCF files...")
    variant_dfs = []
    
    for vcf_file, sample_name in zip(args.input, args.sample_names):
        print(f"  Processing {vcf_file} (sample: {sample_name})")
        try:
            df = parse_vcf_to_dataframe(vcf_file)
            if not df.empty:
                variant_dfs.append(df)
            else:
                print(f"    Warning: No variants found in {vcf_file}")
        except Exception as e:
            print(f"    Error processing {vcf_file}: {e}", file=sys.stderr)
            continue
    
    if not variant_dfs:
        print("Error: No valid variant data found", file=sys.stderr)
        sys.exit(1)
    
    # Merge variants
    print("Merging variant data...")
    merged_variants = merge_variant_dataframes(variant_dfs, args.sample_names)
    
    # Save results
    separator = '\t' if args.format == 'tsv' else ','
    merged_variants.to_csv(args.output, sep=separator, index=False)
    
    # Print summary
    total_variants = len(merged_variants)
    persistent_variants = len(merged_variants[merged_variants['is_persistent']])
    
    print(f"\nMerging completed:")
    print(f"  Total unique variants: {total_variants}")
    print(f"  Persistent variants: {persistent_variants}")
    print(f"  Results saved to: {args.output}")

if __name__ == '__main__':
    main()
