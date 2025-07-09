#!/bin/bash

# Test script for immune neoantigen pipeline
# Author: lilei

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local status=$1
    local message=$2
    case $status in
        "INFO")
            echo -e "${GREEN}[INFO]${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} $message"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} $message"
            ;;
    esac
}

# Function to check if command exists
check_command() {
    local cmd=$1
    if command -v "$cmd" &> /dev/null; then
        print_status "INFO" "$cmd is available"
        return 0
    else
        print_status "ERROR" "$cmd is not available"
        return 1
    fi
}

# Function to check file exists
check_file() {
    local file=$1
    if [[ -f "$file" ]]; then
        print_status "INFO" "File exists: $file"
        return 0
    else
        print_status "ERROR" "File missing: $file"
        return 1
    fi
}

# Function to validate Nextflow syntax
validate_nextflow_syntax() {
    print_status "INFO" "Validating Nextflow syntax..."
    
    if nextflow run main.nf --help &> /dev/null; then
        print_status "INFO" "Nextflow syntax validation passed"
        return 0
    else
        print_status "ERROR" "Nextflow syntax validation failed"
        return 1
    fi
}

# Function to test samplesheet validation
test_samplesheet_validation() {
    print_status "INFO" "Testing samplesheet validation..."
    
    if [[ -f "assets/samplesheet_test.csv" ]]; then
        if python3 bin/metadata_validator.py --input assets/samplesheet_test.csv; then
            print_status "INFO" "Samplesheet validation passed"
            return 0
        else
            print_status "ERROR" "Samplesheet validation failed"
            return 1
        fi
    else
        print_status "WARN" "Test samplesheet not found, skipping validation test"
        return 0
    fi
}

# Function to check Docker/Singularity availability
check_containers() {
    print_status "INFO" "Checking container engines..."
    
    local has_docker=false
    local has_singularity=false
    
    if command -v docker &> /dev/null; then
        print_status "INFO" "Docker is available"
        has_docker=true
    fi
    
    if command -v singularity &> /dev/null; then
        print_status "INFO" "Singularity is available"
        has_singularity=true
    fi
    
    if [[ "$has_docker" == false && "$has_singularity" == false ]]; then
        print_status "WARN" "Neither Docker nor Singularity found. Pipeline will require Conda or manual tool installation."
    fi
}

# Function to run dry run test
run_dry_run() {
    print_status "INFO" "Running pipeline dry run..."
    
    if nextflow run main.nf -profile test --outdir test_results -resume --dry-run &> dry_run.log; then
        print_status "INFO" "Dry run completed successfully"
        return 0
    else
        print_status "ERROR" "Dry run failed. Check dry_run.log for details"
        return 1
    fi
}

# Main test function
main() {
    print_status "INFO" "Starting pipeline tests..."
    
    local test_passed=0
    local test_failed=0
    
    # Check required commands
    print_status "INFO" "Checking required commands..."
    for cmd in nextflow python3; do
        if check_command "$cmd"; then
            ((test_passed++))
        else
            ((test_failed++))
        fi
    done
    
    # Check essential files
    print_status "INFO" "Checking essential files..."
    essential_files=(
        "main.nf"
        "conf/nextflow.config"
        "conf/base.config"
        "lib/utils.nf"
    )
    
    for file in "${essential_files[@]}"; do
        if check_file "$file"; then
            ((test_passed++))
        else
            ((test_failed++))
        fi
    done
    
    # Check workflow files
    print_status "INFO" "Checking workflow files..."
    workflow_files=(
        "workflows/wes_workflow.nf"
        "workflows/rnaseq_workflow.nf"
        "workflows/tcr_workflow.nf"
        "workflows/neoantigen_workflow.nf"
    )
    
    for file in "${workflow_files[@]}"; do
        if check_file "$file"; then
            ((test_passed++))
        else
            ((test_failed++))
        fi
    done
    
    # Check module files
    print_status "INFO" "Checking module files..."
    module_files=(
        "modules/qc.nf"
        "modules/variant_calling.nf"
        "modules/hla_typing.nf"
        "modules/rna_quantification.nf"
        "modules/tcr_analysis.nf"
        "modules/neoantigen_prediction.nf"
    )
    
    for file in "${module_files[@]}"; do
        if check_file "$file"; then
            ((test_passed++))
        else
            ((test_failed++))
        fi
    done
    
    # Validate Nextflow syntax
    if validate_nextflow_syntax; then
        ((test_passed++))
    else
        ((test_failed++))
    fi
    
    # Test samplesheet validation
    if test_samplesheet_validation; then
        ((test_passed++))
    else
        ((test_failed++))
    fi
    
    # Check container engines
    check_containers
    
    # Run dry run if basic tests pass
    if [[ $test_failed -eq 0 ]]; then
        if run_dry_run; then
            ((test_passed++))
        else
            ((test_failed++))
        fi
    else
        print_status "WARN" "Skipping dry run due to previous failures"
    fi
    
    # Print summary
    print_status "INFO" "Test Summary:"
    print_status "INFO" "  Tests passed: $test_passed"
    if [[ $test_failed -gt 0 ]]; then
        print_status "ERROR" "  Tests failed: $test_failed"
        print_status "ERROR" "Pipeline setup has issues that need to be resolved"
        exit 1
    else
        print_status "INFO" "  Tests failed: $test_failed"
        print_status "INFO" "✓ All tests passed! Pipeline is ready to use."
        exit 0
    fi
}

# Run main function
main "$@"
