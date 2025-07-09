#!/usr/bin/env python3

"""
Validate sample metadata and samplesheet format
Author: lilei
"""

import argparse
import pandas as pd
import sys
from pathlib import Path
import re

def validate_sample_id(sample_id):
    """Validate sample ID format"""
    if not sample_id or pd.isna(sample_id):
        return False, "Sample ID cannot be empty"
    
    # Check for valid characters (alphanumeric, underscore, hyphen)
    if not re.match(r'^[A-Za-z0-9_-]+$', str(sample_id)):
        return False, "Sample ID contains invalid characters"
    
    return True, ""

def validate_patient_id(patient_id):
    """Validate patient ID format"""
    if not patient_id or pd.isna(patient_id):
        return False, "Patient ID cannot be empty"
    
    if not re.match(r'^[A-Za-z0-9_-]+$', str(patient_id)):
        return False, "Patient ID contains invalid characters"
    
    return True, ""

def validate_timepoint(timepoint):
    """Validate timepoint format"""
    if not timepoint or pd.isna(timepoint):
        return False, "Timepoint cannot be empty"
    
    valid_timepoints = ['baseline', 'cycle1', 'cycle2', 'cycle3', 'progression', 'response']
    if str(timepoint).lower() not in valid_timepoints:
        return False, f"Invalid timepoint. Must be one of: {', '.join(valid_timepoints)}"
    
    return True, ""

def validate_sample_type(sample_type):
    """Validate sample type"""
    if not sample_type or pd.isna(sample_type):
        return False, "Sample type cannot be empty"
    
    valid_types = ['tumor', 'normal', 'cfDNA', 'PBMC']
    if str(sample_type) not in valid_types:
        return False, f"Invalid sample type. Must be one of: {', '.join(valid_types)}"
    
    return True, ""

def validate_file_path(file_path):
    """Validate file path"""
    if pd.isna(file_path) or file_path == '':
        return True, ""  # Optional files are allowed to be empty
    
    # Check if it's a URL or local path
    if str(file_path).startswith(('http://', 'https://', 'ftp://', 's3://')):
        return True, ""  # Assume URLs are valid
    
    # Check if local file exists
    if not Path(file_path).exists():
        return False, f"File does not exist: {file_path}"
    
    return True, ""

def validate_hla_alleles(hla_alleles):
    """Validate HLA allele format"""
    if pd.isna(hla_alleles) or hla_alleles == '':
        return True, ""  # HLA alleles are optional
    
    # Check HLA format (e.g., HLA-A*02:01;HLA-B*07:02)
    alleles = str(hla_alleles).split(';')
    for allele in alleles:
        if not re.match(r'^HLA-[ABC]\*\d{2}:\d{2}$', allele.strip()):
            return False, f"Invalid HLA allele format: {allele}. Expected format: HLA-A*02:01"
    
    return True, ""

def validate_samplesheet(df):
    """Validate entire samplesheet"""
    errors = []
    warnings = []
    
    # Check required columns
    required_columns = ['sample_id', 'patient_id', 'timepoint', 'sample_type']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        errors.append(f"Missing required columns: {', '.join(missing_columns)}")
        return errors, warnings
    
    # Validate each row
    for idx, row in df.iterrows():
        row_errors = []
        
        # Validate sample_id
        valid, msg = validate_sample_id(row['sample_id'])
        if not valid:
            row_errors.append(f"Row {idx + 1}: {msg}")
        
        # Validate patient_id
        valid, msg = validate_patient_id(row['patient_id'])
        if not valid:
            row_errors.append(f"Row {idx + 1}: {msg}")
        
        # Validate timepoint
        valid, msg = validate_timepoint(row['timepoint'])
        if not valid:
            row_errors.append(f"Row {idx + 1}: {msg}")
        
        # Validate sample_type
        valid, msg = validate_sample_type(row['sample_type'])
        if not valid:
            row_errors.append(f"Row {idx + 1}: {msg}")
        
        # Validate file paths
        file_columns = ['wes_r1', 'wes_r2', 'rna_r1', 'rna_r2', 'tcr_r1', 'tcr_r2']
        for col in file_columns:
            if col in df.columns:
                valid, msg = validate_file_path(row[col])
                if not valid:
                    row_errors.append(f"Row {idx + 1}: {msg}")
        
        # Validate HLA alleles
        if 'hla_alleles' in df.columns:
            valid, msg = validate_hla_alleles(row['hla_alleles'])
            if not valid:
                row_errors.append(f"Row {idx + 1}: {msg}")
        
        # Check data consistency
        sample_type = str(row['sample_type'])
        has_wes = not pd.isna(row.get('wes_r1', '')) and not pd.isna(row.get('wes_r2', ''))
        has_rna = not pd.isna(row.get('rna_r1', '')) and not pd.isna(row.get('rna_r2', ''))
        has_tcr = not pd.isna(row.get('tcr_r1', '')) and not pd.isna(row.get('tcr_r2', ''))
        
        if sample_type == 'normal' and (has_rna or has_tcr):
            warnings.append(f"Row {idx + 1}: Normal samples typically don't have RNA-seq or TCR-seq data")
        
        if sample_type == 'cfDNA' and has_rna:
            warnings.append(f"Row {idx + 1}: cfDNA samples typically don't have RNA-seq data")
        
        if not (has_wes or has_rna or has_tcr):
            row_errors.append(f"Row {idx + 1}: Sample must have at least one data type (WES, RNA-seq, or TCR-seq)")
        
        errors.extend(row_errors)
    
    # Check for duplicate sample IDs
    duplicate_samples = df[df.duplicated('sample_id', keep=False)]['sample_id'].tolist()
    if duplicate_samples:
        errors.append(f"Duplicate sample IDs found: {', '.join(set(duplicate_samples))}")
    
    # Check patient-timepoint combinations
    patient_timepoints = df.groupby(['patient_id', 'timepoint']).size()
    multiple_samples = patient_timepoints[patient_timepoints > 1]
    if not multiple_samples.empty:
        warnings.append(f"Multiple samples per patient-timepoint combination found: {multiple_samples.to_dict()}")
    
    return errors, warnings

def main():
    parser = argparse.ArgumentParser(description='Validate sample metadata and samplesheet')
    parser.add_argument('--input', '-i', required=True,
                       help='Input samplesheet CSV file')
    parser.add_argument('--output', '-o',
                       help='Output validated samplesheet (optional)')
    parser.add_argument('--strict', action='store_true',
                       help='Treat warnings as errors')
    
    args = parser.parse_args()
    
    # Read samplesheet
    try:
        df = pd.read_csv(args.input)
        print(f"Read samplesheet with {len(df)} samples")
    except Exception as e:
        print(f"Error reading samplesheet: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate samplesheet
    errors, warnings = validate_samplesheet(df)
    
    # Print results
    if errors:
        print(f"\nFound {len(errors)} error(s):")
        for error in errors:
            print(f"  ERROR: {error}")
    
    if warnings:
        print(f"\nFound {len(warnings)} warning(s):")
        for warning in warnings:
            print(f"  WARNING: {warning}")
    
    # Determine exit status
    exit_code = 0
    if errors:
        exit_code = 1
    elif warnings and args.strict:
        exit_code = 1
    
    if exit_code == 0:
        print("\n✓ Samplesheet validation passed!")
        
        # Save validated samplesheet if output specified
        if args.output:
            df.to_csv(args.output, index=False)
            print(f"Validated samplesheet saved to: {args.output}")
    else:
        print("\n✗ Samplesheet validation failed!")
    
    sys.exit(exit_code)

if __name__ == '__main__':
    main()
