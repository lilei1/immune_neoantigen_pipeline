#!/bin/bash -euo pipefail
mixcr assemble \
    --report PATIENT_02_T1_cycle1_assemble.report \
     \
    PATIENT_02_T1_cycle1.vdjca \
    PATIENT_02_T1_cycle1.clns

cat <<-END_VERSIONS > versions.yml
"NFCORE_TCRLONGITUDINAL:TCR_LONGITUDINAL:MIXCR_ASSEMBLE":
    mixcr: $(mixcr --version 2>&1 | grep -o 'MiXCR v[0-9.]*' | sed 's/MiXCR v//')
END_VERSIONS
