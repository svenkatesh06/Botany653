#!/bin/bash

## pwd = scripts/coalescent_astral

cd ../../results/coalescent_astral/astral4_nt || exit 1

logfile="_output_astral4.log"

{
    echo "ASTRAL-IV analysis started at $(date)"
    echo "Input gene trees: all_gene_trees_nt.tre"
    echo "Output species tree: astral4_species_tree_nt.tre"
    echo ""

    astral4 \
        -r 16 \
        -s 16 \
        -t 2 \
        --root Ip \
        -i all_gene_trees_nt.tre \
        -o astral4_species_tree_nt.tre

    echo ""
    echo "ASTRAL-IV analysis finished at $(date)"
} 2>&1 | tee -a "$logfile"