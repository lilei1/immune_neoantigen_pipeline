# Changelog

All notable changes to the Immune Neoantigen Pipeline will be documented in this file.

## [2.0.0] - 2025-07-10

### 🚀 Major Features
- **Full ARM64 Support**: Complete compatibility with Apple Silicon (M1/M2/M3) Macs
- **HLA-Enriched Test Data**: Added synthetic HLA reads for comprehensive OptiType testing
- **Enhanced Container Architecture**: All containers verified for ARM64 compatibility

### 🔧 Fixed
- **MiXCR Container Architecture**: Resolved "cannot execute binary file" errors on ARM64
  - Changed from `milaboratory/mixcr:latest` to `mgibio/mixcr:latest`
  - Updated command syntax from v4.x to v3.0.3 format
  - Fixed duplicate `--species` parameter conflicts
- **OptiType Nextflow Execution**: Resolved command line parsing issues
  - Changed from `fred2/optitype:latest` to `umccr/optitype:latest`
  - Fixed environment variable corruption in Docker execution
  - Simplified command structure to avoid argument parsing conflicts
- **GATK4 Compatibility**: Updated to `broadinstitute/gatk:4.3.0.0` for ARM64 support

### 🆕 Added
- **HLA Test Data Generation Script**: `scripts/generate_hla_test_reads.py`
  - Generates synthetic HLA reads based on real gene sequences
  - Supports HLA-A*02:01, HLA-B*07:02, HLA-C*07:01 alleles
  - Configurable read counts and mutation rates
- **New Test Configuration**: `assets/samplesheet_hla.csv`
  - HLA-enriched test samples for OptiType validation
  - 600+ HLA reads per sample with realistic quality scores
- **Comprehensive Documentation**: Updated README with troubleshooting and ARM64 notes

### 🔄 Changed
- **Container Updates**:
  - MiXCR: `milaboratory/mixcr:latest` → `mgibio/mixcr:latest`
  - OptiType: `fred2/optitype:latest` → `umccr/optitype:latest`
  - GATK4: Updated to `broadinstitute/gatk:4.3.0.0`
- **Command Syntax**: Updated MiXCR commands for v3.0.3 compatibility
- **Configuration**: Cleaned up duplicate parameters in modules.config

### 🧪 Testing
- **Verified ARM64 Compatibility**: All containers tested on Apple Silicon
- **Enhanced Test Coverage**: HLA-enriched data provides realistic OptiType testing
- **Improved Error Handling**: Better version detection and error messages

### 📊 Performance
- **Reduced Container Size**: More efficient container selection
- **Faster Startup**: Eliminated container architecture conflicts
- **Better Resource Usage**: Optimized for ARM64 execution

## [1.0.0] - 2025-07-09

### 🎉 Initial Release
- **Core Pipeline**: Nextflow DSL2 implementation
- **Multi-omics Support**: WES, RNA-seq, and TCR-seq workflows
- **Modular Architecture**: Independent, reusable workflow components
- **Container Integration**: Docker/Singularity support
- **Quality Control**: FastQC and MultiQC integration
- **Test Framework**: Basic test data and validation

### 🔧 Core Tools
- **GATK4**: Mutect2 variant calling
- **Salmon**: RNA-seq transcript quantification
- **MiXCR**: TCR clonotype extraction
- **OptiType**: HLA typing
- **BWA**: Read alignment
- **Samtools**: BAM processing

### 📋 Features
- **Sample Management**: CSV-based samplesheet input
- **Multi-environment Support**: Local, SLURM, AWS execution
- **Failure Recovery**: Nextflow checkpointing
- **Metadata Tracking**: Patient-timepoint associations
- **Modular Design**: Independent workflow components

---

## Legend
- 🚀 Major Features
- 🆕 Added
- 🔧 Fixed
- 🔄 Changed
- 🧪 Testing
- 📊 Performance
- 📋 Features
- 🎉 Release
