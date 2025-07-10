#!/bin/bash -euo pipefail
mixcr assemble \
    --report PATIENT_01_T3_progression_assemble.report \
     \
    PATIENT_01_T3_progression.vdjca \
    PATIENT_01_T3_progression.clns

cat <<-END_VERSIONS > versions.yml
"NFCORE_TCRLONGITUDINAL:TCR_LONGITUDINAL:MIXCR_ASSEMBLE":
    mixcr: $(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
END_VERSIONS
