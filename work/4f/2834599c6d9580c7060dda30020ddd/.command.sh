#!/bin/bash -euo pipefail
OptiTypePipeline.py \
    --input SAMPLE_02_T_wes_1.fastq.gz SAMPLE_02_T_wes_2.fastq.gz \
    --dna \
    --prefix SAMPLE_02_T \
    --outdir . \
    --verbose 

cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:WES_WORKFLOW:OPTITYPE":
    optitype: $(OptiTypePipeline.py --version 2>&1 | grep -o 'OptiType [0-9.]*' | cut -d' ' -f2 || echo "unknown")
END_VERSIONS
