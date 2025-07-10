#!/bin/bash -euo pipefail
fastqc --quiet --threads 2 SAMPLE_01_T_rna_1.fastq.gz SAMPLE_01_T_rna_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:RNASEQ_WORKFLOW:FASTQC":
    fastqc: $(fastqc --version | sed 's/FastQC v//')
END_VERSIONS
