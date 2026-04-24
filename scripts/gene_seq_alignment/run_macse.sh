#!/bin/bash

## pwd = Botany563/scripts/gene_seq_alignment
outdir_nt="../../results/msa_macse/nt"
outdir_aa="../../results/msa_macse/aa"
mkdir -p "$outdir_nt" "$outdir_aa"

logfile="../../results/msa_macse/_output_macse_alignment.log"

{
    echo "MACSE alignment started at $(date)"
    echo "NT output directory: $outdir_nt"
    echo "AA output directory: $outdir_aa"
    echo ""

    for f in ../../data/input_to_MSA_perGene_fa/*.fna; do
        base=$(basename "$f" .fna)

        nseq=$(grep -c '^>' "$f")
        if [ "$nseq" -lt 2 ]; then
            echo "SKIP: $f has $nseq sequence(s)" >&2
            continue
        fi

        echo "################################ Processing: $f"

        macse -prog alignSequences \
            -seq "$f" \
            -replace_FS_by_gaps \
            -out_NT "$outdir_nt/${base}.macse.fna" \
            -out_AA "$outdir_aa/${base}.macse.faa"

        echo "alignment done for $f"
        echo ""
    done

    echo ""
    echo "MACSE alignment completed at $(date)"
} 2>&1 | tee -a "$logfile"