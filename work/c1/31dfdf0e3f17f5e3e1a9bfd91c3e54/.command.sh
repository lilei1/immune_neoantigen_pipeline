#!/bin/bash -euo pipefail
mixcr align \
    -s hsa \
    --report PATIENT_01_T4_posttreatment_align.report \
     \
    PATIENT_01_T4_posttreatment_tcr_1.fastq.gz PATIENT_01_T4_posttreatment_tcr_2.fastq.gz \
    PATIENT_01_T4_posttreatment.vdjca

cat <<-END_VERSIONS > versions.yml
"NFCORE_TCRLONGITUDINAL:TCR_LONGITUDINAL:MIXCR_ALIGN":
    mixcr: $(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
END_VERSIONS
