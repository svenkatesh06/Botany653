# MAXIMUML LIKELIHOOD PHYLOGENETIC TREE INFERENCE

Maximum likelihood (ML) is a character-based, model-based method for estimating phylogenetic trees from sequence alignments. Instead of first reducing the alignment to pairwise distances, ML evaluates the aligned characters directly and asks which tree topology, branch lengths, and substitution model parameters make the observed sequence data most probable.
However, ML also has limitations: it depends on the adequacy of the chosen substitution model, can be computationally more intensive than distance methods, and usually uses heuristic tree searches rather than exhaustively evaluating every possible tree.

## CHOICE OF ML TOOL: IQTREE2

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :--- | :--- | :--- | :--- | :--- |
| IQ-TREE 2 | ML phylogenetic inference software that estimates the tree topology, branch lengths, and substitution model parameters that maximize the probability of the observed sequence alignment. It also includes `ModelFinder` for model selection and fast branch-support methods such as ultrafast bootstrap and `SH-aLRT`. | Uses aligned characters directly rather than reducing the data to distances; allows explicit models of molecular evolution; can automatically select best-fitting substitution models; supports partitioned multi-gene analyses; provides fast branch-support estimation; efficient and practical for phylogenomic datasets. | Results depend on the efectiveness of the selected substitution model; heuristic tree searches may not guarantee finding the absolute global optimum; more computationally intensive than distance methods; concatenated analyses may hide gene-tree discordance; branch support can be misleading if alignments, orthology, or model choice are poor. | The alignment is homologous and correctly aligned; the selected substitution model is an adequate approximation of sequence evolution; sites are usually treated as evolving independently; for concatenated analyses, partitions are assumed to share the same underlying species tree; branch lengths and substitution processes adequately describe evolutionary change. | Choice of input alignment type, such as nucleotide or amino acid; partitioning scheme, model-selection option, such as `-m MFP` or `-m MFP+MERGE`; branch-support method, such as ultrafast bootstrap `-B 1000` and SH-aLRT `-alrt 1000`; number of threads using `-T`; rooting strategy after inference. |

### RUNNING IQ-TREE 2

Useful Links:
1. https://iqtree.github.io/doc/Home#documentation
2. https://iqtree.github.io/doc/Quickstart#minimal-command-line-examples
3. https://iqtree.github.io/doc/

```bash
## pwd = scripts/ml_iqtree
conda install -c bioconda::iqtree    # install IQ-TREE2 via conda
mkdir -p ../../results/ml_iqtree     # create output directory for IQ-TREE results
iqtree --help 2>&1 | tee iqtree_help.txt  	 # check IQ-TREE help and save to file
iqtree --version
"IQ-TREE version 3.1.1 for Linux x86 64-bit built Apr  8 2026
Developed by Bui Quang Minh, Thomas Wong, Nhan Ly-Trong, Huaiyan Ren
Contributed by Lam-Tung Nguyen, Dominik Schrempf, Chris Bielow,
Olga Chernomor, Michael Woodhams, Diep Thi Hoang, Heiko Schmidt"
```
The installed IQ-TREE command supports directory input using `-p` FILE|DIR, automatic model selection using `-m MFP`, partition merging using `--merge` or `MFP+MERGE`, ultrafast bootstrap using `-B`, SH-aLRT support using `--alrt`, outgroup rooting using `-o`, and reproducible runs using `--seed`.

In IQ-TREE, `-p` can take either a partition file or a directory of alignment files. When a directory is supplied, IQ-TREE treats each alignment file as a separate initial partition and concatenates the alignments internally. This is useful for multi-gene datasets like this one because it avoids manually calculating gene start and end coordinates in a supermatrix.

The `-m MFP+MERGE` option performs automatic model selection and partition-scheme optimization. `MFP` means ModelFinder Plus, where IQ-TREE tests substitution models and chooses the best-fitting model before tree inference. `+MERGE` tells IQ-TREE to start with the initial partitions — in this case, one partition per gene file — and then merge partitions with similar evolutionary patterns if doing so improves the model fit.
> This is useful because different genes may evolve under different rates or substitution patterns, but keeping every gene as a fully separate partition can over-parameterize the analysis, especially with only 5 taxa.

The `-B` 1000 option runs 1000 ultrafast bootstrap replicates. These values provide branch-support estimates for the inferred maximum-likelihood tree. The `-alrt` 1000 option additionally calculates SH-aLRT support values using 1000 replicates.
>This will give us a clearer picture of which branches are well-supported and which are not, especially since we have a small number of taxa and a large number of genes.

The `-o Ip` option writes the final tree with Ip as the outgroup. This is appropriate because ML tree searches produce unrooted trees by default.

#### WHY PARTITIONING?
The dataset contains multiple genes, and each gene may have a different evolutionary rate or substitution pattern. If all genes were forced into one unpartitioned model, the analysis would assume that every site across all genes evolved under the same model. That is probably too simplistic for a multi-gene coding-sequence dataset.



```bash
chmod +x run*_iqtree_partitioned.sh
## run iqtree on NT sequences with partitioning
./run_nt_iqtree_partitioned.sh

## run iqtree on AA sequences with partitioning
./run_aa_iqtree_partitioned.sh
```

I ran both because nucleotide and amino acid sequences can capture slightly different evolutionary signals. Nucleotide sequences retain synonymous substitutions and contain more raw variation, which can be useful for closely related taxa. However, nucleotide alignments can also be more affected by multiple substitutions, codon-position effects etc., especially across more divergent species. *AA sequences remove variations and focus on protein-level changes, which may provide a more conservative signal for deeper evolutionary relationships.*
In this project, this comparison is useful because the dataset is based on protein-coding electric-organ-associated genes from multiple species.

### RESULTS AND INTERPRETATION


