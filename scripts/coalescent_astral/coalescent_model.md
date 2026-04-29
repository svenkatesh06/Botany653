# THE COALESCENT MODEL

Useful Links:
1. https://github.com/chaoszhang/ASTER
2. https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md
3. https://academic.oup.com/mbe/article/42/8/msaf172/8203509?login=true

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :--- | :--- | :--- | :--- | :--- |
| ASTRAL-IV / ASTER | Coalescent-based species-tree method that estimates a species tree from a collection of gene-tree topologies. | Accounts for gene-tree discordance under the multi-species coalescent; appropriate for incomplete lineage sorting; fast for this dataset. | Does not use the original sequence alignments directly; depends on the quality of the input gene trees. | Gene-tree discordance is modeled mainly under incomplete lineage sorting, and the input trees represent single-copy orthologous markers. | I used ASTRAL-IV from ASTER with 50 per-gene IQ-TREE trees, rooted the output with `Ip`, and used `-r 16 -s 16` for a more thorough search. |


## PREPARING THE GENE TREES

For the coalescent analysis, I first inferred a **maximum-likelihood tree per MACSE-cleaned nucleotide gene alignment using IQ-TREE with automatic model selection**. This is different from the concatenated IQ-TREE analysis, where all genes were analyzed together as one partitioned dataset. ASTRAL-IV requires a collection of gene trees as input, so each gene alignment was analyzed separately.

The resulting per-gene trees were concatenated into a single gene-tree file and used as input to ASTRAL-IV. The gene names did not need to be included inside the Newick trees because ASTRAL treats each tree as 1 gene tree. However, the tip labels must correspond to the same species/taxon names across trees which is the case here.

```bash
# pwd = scripts/coalescent_astral
bash run_iqtree_per_gene.sh

cat ../../results/coalescent_astral/gene_trees_nt/*.treefile \
    > ../../results/coalescent_astral/astral4_nt/all_gene_trees_nt.tre

conda install -c bioconda aster -y
```


## RUN ASTRAL-IV

I used ASTRAL-IV from the ASTER package for the coalescent species-tree analysis. ASTRAL-IV estimates a species tree from a set of gene-tree topologies and is designed to account for gene-tree discordance under the multi-species coalescent model.

The species tree was rooted with Ip, the non-electric outgroup. ASTRAL branch support values were interpreted as local posterior probabilities, not bootstrap values.

Because this dataset contains only 5 taxa and a modest number of genes, I used the more exhaustive ASTRAL-IV search setting -r 16 -s 16, as recommended for small datasets. This increases the search effort while remaining computationally inexpensive for this dataset.

```bash
# pwd = scripts/coalescent_astral

bash run_astral4.sh
```
Important log output:

```text
Accurate Species TRee ALgorithm IV (ASTRAL-IV)
*** NOW with integrated CASTLES-2 ***
Version: v1.25.4.8
#Genetrees: 50
#Species: 5
#Rounds: 16
#Samples: 16
#Threads: 2

Current score: 232
Current tree: (((Ee,Bg),(Bb,Pk)),Ip);
Final Tree: (((Ee,Bg),(Bb,Pk)),Ip);
Score: 232
```

ASTRAL-IV recovered the same major topology as the concatenated IQ-TREE and MrBayes analyses:

`(((Ee,Bg),(Bb,Pk)),Ip);`

This places `Ip` as the outgroup and separates the electric fishes into two major lineages:

Ee + Bg: gymnotiform electric fishes
Bb + Pk: mormyroid electric fishes

Both major ingroup clades had local posterior probability support of 1.0.

### VISUALIZE THE ASTRAL-IV SPECIES TREE

I visualized the ASTRAL-IV species tree using the `visualize_trees_astral.ipynb` notebook. The tree was rooted with `Ip` as the outgroup and plotted with tip labels colored by lineage. The resulting plot is saved in `Botany653/results/coalescent_astral/astral4_nt/tree_plots/astral4_nt_rooted.pdf`.


The ASTRAL-IV tree recovered the same topology as the concatenated IQ-TREE and the codon-partitioned MrBayes analyses, but the branch lengths are different across methods which makes sense. IQ-TREE and MrBayes estimate branch lengths from the partitioned/concatenated sequence alignment, whereas ASTRAL-IV estimates a species tree from independently inferred gene trees under a coalescent framework.

