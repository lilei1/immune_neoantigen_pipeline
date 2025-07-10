#!/bin/bash -euo pipefail
# Extract transcript sequences from genome
grep "^>" genome.fasta | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat genome.fasta > gentrome.fa

# Build Salmon index
salmon index \
    --threads 2 \
    --transcripts gentrome.fa \
    --decoys decoys.txt \
    --index salmon \


cat <<-END_VERSIONS > versions.yml
"IMMUNE_NEOANTIGEN_PIPELINE:RNASEQ_WORKFLOW:SALMON_INDEX":
    salmon: $(echo $(salmon --version) | sed -e "s/salmon //g")
END_VERSIONS
