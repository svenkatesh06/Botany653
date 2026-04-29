#!/bin/bash

## pwd = scripts/coalescent_astral

mkdir -p ../../results/coalescent_astral/gene_trees_nt
mkdir -p ../../results/coalescent_astral/astral4_nt

indir="../../results/msa_macse/reordered_cleaned_macse_nt"
outdir="../../results/coalescent_astral/gene_trees_nt"

logfile="$outdir/_output_iqtree.log"

{
	echo "IQ-TREE  maximum likelihood analysis started at $(date)"

for f in "$indir"/*.fna; do
    base=$(basename "$f" .fna)

    iqtree \
        -s "$f" \
        --seqtype DNA \
        -m MFP \
        -B 1000 \
        -alrt 1000 \
        -T 2 \
        --seed 12345 \
		--verbose \
        --prefix "$outdir/$base"
done
	echo ""
	echo "IQ-TREE partitioned maximum likelihood analysis finished at $(date)"

} 2>&1 | tee -a "$logfile"