#!/bin/bash

## pwd = Botany653/scripts/ml_iqtree

indir="../../results/msa_macse/reordered_cleaned_macse_nt"
outdir="../../results/ml_iqtree/nt_partitioned/"
mkdir -p "$outdir"

logfile="$outdir/_output_iqtree_partitioned.log"

{
	echo "IQ-TREE partitioned maximum likelihood analysis started at $(date)"
	echo "Input directory: $indir"
	echo "Output directory: $outdir"
	echo ""
	echo "Running NT partitioned ML analysis"
	echo ""

	iqtree \
		-p "$indir" \
		--seqtype DNA \
		-m MFP+MERGE \
		-B 1000 \
		-alrt 1000 \
		-T 2 \
		-o Ip \
		--seed 12345 \
		--verbose \
		--prefix "$outdir/nt_partitioned"

	echo ""
	echo "NT partitioned ML analysis completed"
	echo ""

	echo "IQ-TREE partitioned maximum likelihood analysis completed at $(date)"

} 2>&1 | tee -a "$logfile"