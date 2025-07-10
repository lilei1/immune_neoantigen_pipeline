#!/bin/bash -euo pipefail
fastqc --quiet --threads 2 SAMPLE_02_T_wes_1.fastq.gz SAMPLE_02_T_wes_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:WES_WORKFLOW:FASTQC":
    fastqc: $(fastqc --version | sed 's/FastQC v//')
END_VERSIONS
