# Setup Guide for Immune Repertoire + Neoantigen Pipeline

This guide will help you set up the immune repertoire and neoantigen prediction pipeline on your system.

## Prerequisites

### 1. Install Nextflow

```bash
# Install Java (required for Nextflow)
sudo apt-get update
sudo apt-get install openjdk-11-jdk

# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Install Container Engine

Choose one of the following:

#### Option A: Docker
```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER
# Log out and back in for group changes to take effect
```

#### Option B: Singularity
```bash
# Install Singularity (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y singularity-container
```

### 3. Install Python Dependencies (for utility scripts)

```bash
pip3 install pandas numpy biopython pysam pyvcf matplotlib seaborn
```

### 4. Install R Dependencies (for clonotype tracking)

```bash
R -e "install.packages(c('dplyr', 'ggplot2', 'readr', 'tidyr', 'stringr', 'RColorBrewer'), repos='https://cran.rstudio.com/')"
```

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/lilei1/immune_neoantigen_pipeline.git
cd immune_neoantigen_pipeline
```

### 2. Test the Pipeline Setup

```bash
# Run the test script
./bin/test_pipeline.sh
```

### 3. Run with Test Data

```bash
# Run with test profile and Docker
nextflow run main.nf -profile test,docker

# Run with test profile and Singularity
nextflow run main.nf -profile test,singularity
```

### 4. Run with Your Data

```bash
# Prepare your samplesheet (see example in assets/samplesheet_test.csv)
# Validate your samplesheet
python3 bin/metadata_validator.py --input your_samplesheet.csv

# Run the pipeline
nextflow run main.nf --input your_samplesheet.csv --outdir results -profile docker
```

## Configuration Profiles

### Local Execution
- `docker`: Use Docker containers locally
- `singularity`: Use Singularity containers locally

### HPC Execution
- `slurm`: Run on SLURM cluster (requires configuration)

### Cloud Execution
- `aws`: Run on AWS Batch (requires AWS setup)

### Testing
- `test`: Run with small test dataset

## Customizing the Pipeline

### 1. Modify Parameters

Create a custom parameters file:

```yaml
# custom_params.yaml
run_wes: true
run_rnaseq: true
run_tcr: false
peptide_lengths: "9,10"
binding_threshold: 1000
```

Run with custom parameters:
```bash
nextflow run main.nf -params-file custom_params.yaml --input samplesheet.csv --outdir results
```

### 2. Modify Resource Requirements

Edit `conf/base.config` to adjust CPU, memory, and time requirements for each process.

### 3. Add Custom Modules

1. Create new module in `modules/` directory
2. Include in appropriate workflow in `workflows/` directory
3. Update configuration in `conf/modules.config`

## Troubleshooting

### Common Issues

1. **Nextflow not found**
   - Ensure Nextflow is installed and in your PATH
   - Check Java installation: `java -version`

2. **Container engine issues**
   - Verify Docker/Singularity installation
   - Check permissions for Docker

3. **Memory/CPU issues**
   - Adjust resource requirements in `conf/base.config`
   - Use appropriate profile for your system

4. **File permission errors**
   - Ensure input files are readable
   - Check output directory permissions

### Getting Help

1. Check the pipeline logs in the `work/` directory
2. Use Nextflow's built-in debugging: `nextflow run main.nf -with-trace -with-timeline -with-report`
3. Open an issue on GitHub with error logs

## Pipeline Architecture

```
main.nf
├── workflows/
│   ├── wes_workflow.nf      # WES processing
│   ├── rnaseq_workflow.nf   # RNA-seq processing
│   ├── tcr_workflow.nf      # TCR-seq processing
│   └── neoantigen_workflow.nf # Neoantigen prediction
├── modules/
│   ├── qc.nf                # Quality control
│   ├── variant_calling.nf   # Variant calling
│   ├── hla_typing.nf        # HLA typing
│   ├── rna_quantification.nf # RNA quantification
│   ├── tcr_analysis.nf      # TCR analysis
│   └── neoantigen_prediction.nf # Neoantigen prediction
├── conf/
│   ├── nextflow.config      # Main configuration
│   ├── base.config          # Base process configuration
│   ├── modules.config       # Module-specific configuration
│   ├── test.config          # Test configuration
│   ├── slurm.config         # SLURM configuration
│   └── aws.config           # AWS configuration
└── bin/
    ├── track_clones.R       # Clonotype tracking script
    ├── merge_variants.py    # Variant merging script
    ├── metadata_validator.py # Samplesheet validation
    └── test_pipeline.sh     # Pipeline testing script
```

## Data Requirements

### Input Samplesheet Format

The pipeline expects a CSV file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| sample_id | Yes | Unique sample identifier |
| patient_id | Yes | Patient identifier |
| timepoint | Yes | Collection timepoint |
| sample_type | Yes | Sample type (tumor/normal/cfDNA/PBMC) |
| wes_r1, wes_r2 | Optional | WES FASTQ files |
| rna_r1, rna_r2 | Optional | RNA-seq FASTQ files |
| tcr_r1, tcr_r2 | Optional | TCR-seq FASTQ files |
| hla_alleles | Optional | Known HLA alleles |

### Reference Files

The pipeline can use:
- iGenomes references (automatic download)
- Custom reference files (specify with --fasta and --gtf)
- Pre-built indices (specify with --salmon_index)

## Performance Optimization

### For Large Datasets
- Use appropriate resource profiles (slurm, aws)
- Increase memory and CPU allocations
- Use pre-built indices when possible

### For Small Datasets
- Use test profile for quick validation
- Reduce resource requirements
- Skip time-consuming QC steps if needed

## Security Considerations

- Ensure input data is properly secured
- Use appropriate access controls for cloud execution
- Validate all input files before processing
- Monitor resource usage to prevent abuse
