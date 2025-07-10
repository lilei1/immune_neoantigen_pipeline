#!/bin/bash -euo pipefail
mixcr analyze shotgun \
    -s hsa \
    --starting-material rna \
    --only-productive \
    --report SAMPLE_01_T.report \
     \
    SAMPLE_01_T_tcr_1.fastq.gz SAMPLE_01_T_tcr_2.fastq.gz \
    SAMPLE_01_T

cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:TCR_WORKFLOW:MIXCR_ANALYZE":
    mixcr: $(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
END_VERSIONS
