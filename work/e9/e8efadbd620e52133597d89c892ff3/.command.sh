#!/bin/bash -euo pipefail
OptiTypePipeline.py \
    --input SAMPLE_01_cfDNA_wes_1.fastq.gz SAMPLE_01_cfDNA_wes_2.fastq.gz \
    --dna \
    --prefix SAMPLE_01_cfDNA \
    --outdir . \
    --verbose 

cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:WES_WORKFLOW:OPTITYPE":
    optitype: $(OptiTypePipeline.py --version 2>&1 | grep -o 'OptiType [0-9.]*' | cut -d' ' -f2 || echo "unknown")
END_VERSIONS
