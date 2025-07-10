#!/bin/bash -euo pipefail
bwa index genome.fasta

cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:WES_WORKFLOW:BWA_INDEX":
    bwa: $(echo $(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*$//')
END_VERSIONS
