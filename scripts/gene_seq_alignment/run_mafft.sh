#!/bin/bash

## pwd = Botany563/scripts/gene_seq_alignment
outdir="../../results/msa_mafft"
mkdir -p $outdir
logfile="$outdir/_output_mafft_alignment.log"

{
    echo "MAFFT alignment started at $(date)"
    echo "Output directory: $outdir"
    echo ""

    for f in ../../data/input_to_MSA_perGene_fa/*.fna; do
        base=$(basename "$f" .fna)

		nseq=$(grep -c '^>' "$f")
    	if [ "$nseq" -lt 2 ]; then
        	echo "SKIP: $f has $nseq sequence(s)" >&2
        	continue
    	fi
        echo "################################ Processing: $f"

        mafft --localpair --maxiterate 1000 "$f" > "$outdir/${base}.mafft.fna"

        echo "alignment done for $f"
		echo ""
    done

    echo ""
    echo "MAFFT alignment completed at $(date)"
} 2>&1 | tee -a "$logfile"