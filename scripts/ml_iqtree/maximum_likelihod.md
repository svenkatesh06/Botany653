# MAXIMUML LIKELIHOOD PHYLOGENETIC TREE INFERENCE

Maximum likelihood (ML) is a character-based, model-based method for estimating phylogenetic trees from sequence alignments. Instead of first reducing the alignment to pairwise distances, ML evaluates the aligned characters directly and asks which tree topology, branch lengths, and substitution model parameters make the observed sequence data most probable.
However, ML also has limitations: it depends on the adequacy of the chosen substitution model, can be computationally more intensive than distance methods, and usually uses heuristic tree searches rather than exhaustively evaluating every possible tree.

## CHOICE OF ML TOOL: IQTREE2

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :--- | :--- | :--- | :--- | :--- |
| IQ-TREE 2 | ML phylogenetic inference software that estimates the tree topology, branch lengths, and substitution model parameters that maximize the probability of the observed sequence alignment. It also includes `ModelFinder` for model selection and fast branch-support methods such as ultrafast bootstrap and `SH-aLRT`. | Uses aligned characters directly rather than reducing the data to distances; allows explicit models of molecular evolution; can automatically select best-fitting substitution models; supports partitioned multi-gene analyses; provides fast branch-support estimation; efficient and practical for phylogenomic datasets. | Results depend on the efectiveness of the selected substitution model; heuristic tree searches may not guarantee finding the absolute global optimum; more computationally intensive than distance methods; concatenated analyses may hide gene-tree discordance; branch support can be misleading if alignments, orthology, or model choice are poor. | The alignment is homologous and correctly aligned; the selected substitution model is an adequate approximation of sequence evolution; sites are usually treated as evolving independently; for concatenated analyses, partitions are assumed to share the same underlying species tree; branch lengths and substitution processes adequately describe evolutionary change. | Choice of input alignment type, such as nucleotide or amino acid; whether to analyze individual genes or a concatenated supermatrix; partitioning scheme, such as one partition per gene; model-selection option, such as `-m MFP` or `-m MFP+MERGE`; branch-support method, such as ultrafast bootstrap `-B 1000` and SH-aLRT `-alrt 1000`; number of threads using `-nt AUTO`; rooting strategy after inference. |

### RUNNING IQ-TREE 2




