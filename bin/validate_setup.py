#!/usr/bin/env python3

"""
Validate pipeline setup without requiring Nextflow
Author: lilei
"""

import os
import sys
import subprocess
from pathlib import Path

def check_file_exists(file_path, description=""):
    """Check if a file exists"""
    if Path(file_path).exists():
        print(f"✓ {description or file_path}")
        return True
    else:
        print(f"✗ {description or file_path} - NOT FOUND")
        return False

def check_command_exists(command):
    """Check if a command is available"""
    try:
        subprocess.run([command, '--version'], capture_output=True, check=True)
        print(f"✓ {command} is available")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        try:
            subprocess.run([command, '-v'], capture_output=True, check=True)
            print(f"✓ {command} is available")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"✗ {command} is NOT available")
            return False

def validate_python_packages():
    """Validate required Python packages"""
    required_packages = [
        ('pandas', 'pandas'),
        ('numpy', 'numpy'),
        ('biopython', 'Bio')  # biopython package imports as Bio
    ]
    missing_packages = []

    for display_name, import_name in required_packages:
        try:
            __import__(import_name)
            print(f"✓ Python package: {display_name}")
        except ImportError:
            print(f"✗ Python package: {display_name} - NOT FOUND")
            missing_packages.append(display_name)

    return len(missing_packages) == 0

def validate_r_packages():
    """Validate required R packages"""
    required_packages = ['dplyr', 'ggplot2', 'readr']
    
    try:
        # Create R script to check packages
        r_script = """
        packages <- c('dplyr', 'ggplot2', 'readr')
        missing <- packages[!packages %in% installed.packages()[,'Package']]
        if(length(missing) > 0) {
            cat('MISSING:', paste(missing, collapse=','), '\\n')
            quit(status=1)
        } else {
            cat('ALL_FOUND\\n')
        }
        """
        
        result = subprocess.run(['R', '--slave', '-e', r_script], 
                              capture_output=True, text=True)
        
        if result.returncode == 0 and 'ALL_FOUND' in result.stdout:
            print("✓ Required R packages are available")
            return True
        else:
            print("✗ Some R packages are missing")
            if 'MISSING:' in result.stdout:
                missing = result.stdout.split('MISSING:')[1].strip()
                print(f"  Missing packages: {missing}")
            return False
            
    except FileNotFoundError:
        print("✗ R is not available")
        return False

def validate_file_structure():
    """Validate pipeline file structure"""
    required_files = [
        ('main.nf', 'Main pipeline file'),
        ('conf/nextflow.config', 'Main configuration'),
        ('conf/base.config', 'Base configuration'),
        ('conf/modules.config', 'Module configuration'),
        ('conf/test.config', 'Test configuration'),
        ('assets/samplesheet_test.csv', 'Test samplesheet'),
        ('workflows/wes_workflow.nf', 'WES workflow'),
        ('workflows/rnaseq_workflow.nf', 'RNA-seq workflow'),
        ('workflows/tcr_workflow.nf', 'TCR workflow'),
        ('workflows/neoantigen_workflow.nf', 'Neoantigen workflow'),
        ('modules/qc.nf', 'QC module'),
        ('modules/variant_calling.nf', 'Variant calling module'),
        ('modules/hla_typing.nf', 'HLA typing module'),
        ('modules/rna_quantification.nf', 'RNA quantification module'),
        ('modules/tcr_analysis.nf', 'TCR analysis module'),
        ('modules/neoantigen_prediction.nf', 'Neoantigen prediction module'),
        ('bin/track_clones.R', 'Clonotype tracking script'),
        ('bin/merge_variants.py', 'Variant merging script'),
        ('bin/metadata_validator.py', 'Metadata validator'),
    ]
    
    print("\nValidating file structure:")
    all_found = True
    
    for file_path, description in required_files:
        if not check_file_exists(file_path, description):
            all_found = False
    
    return all_found

def validate_samplesheet():
    """Validate test samplesheet"""
    print("\nValidating test samplesheet:")
    
    samplesheet_path = "assets/samplesheet_test.csv"
    if not Path(samplesheet_path).exists():
        print(f"✗ Test samplesheet not found: {samplesheet_path}")
        return False
    
    try:
        # Run metadata validator
        result = subprocess.run([
            'python3', 'bin/metadata_validator.py', 
            '--input', samplesheet_path
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Test samplesheet validation passed")
            return True
        else:
            print("✗ Test samplesheet validation failed")
            print(result.stdout)
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"✗ Error validating samplesheet: {e}")
        return False

def main():
    print("Immune Repertoire + Neoantigen Pipeline Setup Validation")
    print("=" * 60)
    
    all_checks_passed = True
    
    # Check basic commands
    print("\nChecking basic commands:")
    commands = ['python3', 'R']
    for cmd in commands:
        if not check_command_exists(cmd):
            all_checks_passed = False
    
    # Check optional commands
    print("\nChecking optional commands:")
    optional_commands = ['nextflow', 'docker', 'singularity']
    for cmd in optional_commands:
        check_command_exists(cmd)  # Don't fail if these are missing
    
    # Validate Python packages
    print("\nValidating Python packages:")
    if not validate_python_packages():
        all_checks_passed = False
        print("  Install missing packages with: pip3 install pandas numpy biopython")
    
    # Validate R packages
    print("\nValidating R packages:")
    if not validate_r_packages():
        all_checks_passed = False
        print("  Install missing packages with: R -e \"install.packages(c('dplyr', 'ggplot2', 'readr'))\"")
    
    # Validate file structure
    if not validate_file_structure():
        all_checks_passed = False
    
    # Validate samplesheet
    if not validate_samplesheet():
        all_checks_passed = False
    
    # Print summary
    print("\n" + "=" * 60)
    if all_checks_passed:
        print("✓ ALL CHECKS PASSED!")
        print("The pipeline setup is valid and ready to use.")
        print("\nNext steps:")
        print("1. Install Nextflow if not already installed")
        print("2. Install Docker or Singularity for containerization")
        print("3. Run: nextflow run main.nf -profile test,docker")
    else:
        print("✗ SOME CHECKS FAILED!")
        print("Please address the issues above before running the pipeline.")
        sys.exit(1)

if __name__ == '__main__':
    main()
