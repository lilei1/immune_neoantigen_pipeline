#!/bin/bash -euo pipefail
# Filter clones by count and frequency
mixcr filter \
    --min-count 2 \
    --min-fraction 0.00001 \
    --report PATIENT_01_T3_progression_filter.report \
     \
    PATIENT_01_T3_progression.clns \
    PATIENT_01_T3_progression_filtered.clns

cat <<-END_VERSIONS > versions.yml
"NFCORE_TCRLONGITUDINAL:TCR_LONGITUDINAL:MIXCR_FILTER":
    mixcr: $(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
END_VERSIONS
