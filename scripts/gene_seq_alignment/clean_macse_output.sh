#!/bin/bash

## pwd = Botany563/scripts/gene_seq_alignment

indir_nt="../../results/msa_macse/nt"
outdir_nt="../../results/msa_macse/cleaned_nt"
outdir_aa="../../results/msa_macse/cleaned_aa"

mkdir -p "$outdir_nt" "$outdir_aa"

logfile="../../results/msa_macse/_output_cleanup_macse.log"

{
    echo "MACSE export/cleanup started at $(date)"
    echo "Input NT directory: $indir_nt"
    echo "Cleaned NT output directory: $outdir_nt"
    echo "Cleaned AA output directory: $outdir_aa"
    echo ""

    for f in "$indir_nt"/*.fna; do
        [ -e "$f" ] || continue

        base=$(basename "$f" .macse.fna)

        echo "################################ Processing: $f"

        macse -prog exportAlignment \
            -align "$f" \
            -codonForInternalStop NNN \
            -codonForInternalFS --- \
            -codonForExternalFS --- \
            -charForRemainingFS - \
            -out_NT "$outdir_nt/${base}.macse.cleaned.fna" \
            -out_AA "$outdir_aa/${base}.macse.cleaned.faa"

        echo "cleanup/export done for $f"
        echo ""
    done

    echo ""
    echo "MACSE export/cleanup completed at $(date)"
} 2>&1 | tee -a "$logfile"