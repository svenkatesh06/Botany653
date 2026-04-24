# DISTANCE AND PARSIMONY METHODS

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :--- | :--- | :--- | :--- | :--- |
| `ape` [distance method is BioNJ] | Distance-based phylogenetic method. Aligned DNA sequences are converted into a pairwise genetic distance matrix. Then BioNJ constructs an unrooted tree using a neighbor-joining framework, but incorporates information about the variance of distance estimates when choosing which taxa/clusters to join. | Fast and simple; works well as an initial tree estimate; does not require a strict molecular clock; useful for comparing overall genetic divergence among the five taxa. | Reduces the full alignment to pairwise distances, so site-pattern information is partly lost; tree quality depends strongly on the chosen distance model and alignment quality; does not explicitly optimize likelihood or parsimony; branch lengths/topology can be affected by missing data and gaps. | Homologous sites are correctly aligned; the selected substitution model reasonably describes sequence divergence; distances are sufficiently tree-like/additive; each gene tree reflects the evolutionary signal in that gene; the consensus tree summarizes clades repeatedly recovered across genes. | Choice of input alignment - MACSE NT alignments; DNA distance model; treatment of gaps/missing data; root |
| `phangorn` [parsimony method is maximum parsimony using `optim.parsimony()`] | Parsimony-based phylogenetic method. Aligned DNA sequences are treated as discrete character states, and the method searches for the tree that requires the fewest total nucleotide changes. Following the class lecture, a raw-distance NJ tree is first used as the starting tree, and `optim.parsimony()` improves this tree under the parsimony criterion. | Simple and intuitive; uses site-pattern information directly rather than reducing the alignment only to pairwise distances; useful as a comparison to the BioNJ distance-based tree; does not require specifying an explicit nucleotide substitution model. | Can be sensitive to long-branch attraction; does not model unequal substitution rates across sites or lineages; results can depend on the starting tree and tree-search strategy; alignment errors and poorly informative sites can affect the inferred topology. | Homologous sites are correctly aligned; the best tree is the one requiring the fewest evolutionary changes; nucleotide characters are informative about shared ancestry; repeated clades across per-gene trees reflect phylogenetic signal rather than convergence, rate variation, or alignment artifacts. | Choice of input alignment - MACSE NT alignments; starting tree - raw-distance NJ tree; parsimony optimization method - `optim.parsimony()`; treatment of gaps/missing data; root |

## DISTANCE-BASED TREE

Useful links:

1. https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf
2. https://cran.r-project.org/web/packages/ape/ape.pdf

I chose BioNJ rather than WPGMA because WPGMA assumes a clock-like or ultrametric tree, where all taxa are equally distant from the root. I don't believe this assumption will hold for the electric fish data, especially because the taxa come from different lineages and may have unequal rates of sequence evolution. BioNJ is still a distance-based method and therefore has limitations, but like Neighbor-Joining, it does not require equal evolutionary rates across lineages.

### BioNJ based tree using `ape` and `phangorn`

For the distance-based analysis, I used the MACSE nucleotide alignments because MACSE is designed for coding sequences and accounts for frameshifts and stop codons during alignment. Instead of concatenating the genes into a single supermatrix, I estimated one pairwise DNA distance matrix per gene using `dist.dna()` and inferred one per-gene BioNJ tree using `bionj()`. I then summarized the collection of per-gene trees using a majority-rule consensus tree.

> > Because this analysis estimates 1 tree per MACSE-aligned gene, the individual outputs are gene trees rather than one concatenated species tree. I used a majority-rule consensus tree to summarize the clades that were repeatedly recovered across genes. This consensus tree can be interpreted as a species-level summary of the selected marker-gene signal, but not as a definitive species tree. Individual genes can have different evolutionary histories due to incomplete lineage sorting, gene duplication/loss, selection, or missing-data/alignment effects.

I tested 2 NT substitution models for estimating the distance matrix: **F84 and TN93**. Both models are more flexible than the simplest distance models because they allow unequal base frequencies. This is useful for real nucleotide sequence data, where A, C, G, and T are often not present at equal frequencies.
- F84 also allows transitions and transversions to occur at different rates.
- TN93 is even more flexible because it allows unequal base frequencies and also allows the two classes of transitions, A↔G and C↔T, to have different rates.

Finally, I used TN93 as the main distance model and F84 as a comparison to check whether the per-gene BioNJ trees and resulting consensus topology were sensitive to the choice of distance model.

In `dist.dna()`, I used `pairwise.deletion = TRUE`. Although the input is a multiple sequence alignment, the MACSE alignments can still contain gaps, ambiguous characters, or missing positions for some taxa. Pairwise deletion means that, for each pair of taxa, the distance is calculated using only the sites that are available and comparable for that particular pair. This allows more of the alignment to contribute to the distance estimates instead of removing every column that contains a gap or missing value in any taxon. This is useful for per-gene alignments because a small number of problematic sites in one taxon would otherwise cause those sites to be excluded for all taxa within that gene alignment.

BioNJ estimates unrooted trees from the distance matrices. I first inferred the per-gene BioNJ trees as unrooted trees, then created a majority-rule consensus tree from them. To make the final consensus tree biologically interpretable, I rooted the consensus tree using *Ictalurus punctatus* as the outgroup because it is the non-electric fish in this dataset. Rooting with the outgroup places the root on the branch separating *I. punctatus* from the electric fish taxa and provides directionality for interpreting relationships among the electric lineages. After rooting, I ladderized the consensus tree to make the plotted topology easier to read. Ladderizing only changes the visual ordering of branches; it does not change the inferred relationships or branch lengths.

```bash
## pwd = scripts/distance_parsimony
mkdir -p ../../results/distance_parsimony/BioNJ
Rscript distance_based_tree.R 2>&1 | tee ../../results/distance_parsimony/BioNJ/_BioNJ_run.log
```

**The F84 and TN93 analyses produced the same overall consensus topology.** For both distance models, per-gene BioNJ trees were inferred and then summarized using a majority-rule consensus tree. In both consensus trees, Ee grouped with Bg, and Bb grouped with Pk, while Ip was separated as the outgroup. This indicates that the distance-based consensus result is not strongly sensitive to whether F84 or TN93 is used.

**The recovered consensus topology is biologically reasonable.** Ee and Bg are both gymnotiform electric fishes, while Bb and Pk are both mormyrid electric fishes. Ip, the channel catfish, was used as the non-electric outgroup to root the consensus tree. The BioNJ consensus tree grouped the two gymnotiforms together and the two mormyrids together, which is consistent with their known evolutionary relationships.

The distance-based per-gene consensus approach was able to capture the repeated genetic divergence patterns among the taxa, even though BioNJ reduces each gene alignment to pairwise distances and does not explicitly optimize likelihood or parsimony.

## PARSIMONY-BASED TREE

Useful links:

1. https://cran.r-project.org/web/packages/phangorn/refman/phangorn.html
2. https://klausvigo.github.io/phangorn/
3. https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf

For the parsimony-based analysis, I used the same MACSE nucleotide alignments used for the distance-based analysis. Maximum parsimony does not first reduce the alignment to a corrected genetic distance matrix. Instead, it treats aligned nucleotide sites as discrete characters and searches for the tree requiring the fewest evolutionary changes.

I first generated a raw-distance NJ tree as the starting tree. I then converted the DNA alignment to a `phyDat` object and used `parsimony()` to calculate the starting parsimony score. Finally, I used `optim.parsimony()` from `phangorn` to improve the starting tree under the maximum-parsimony criterion. I repeated this for each MACSE-aligned gene and summarized the per-gene maximum-parsimony trees using a majority-rule consensus tree.

Maximum parsimony is simple and directly uses aligned site-pattern information, but it does not use an explicit nucleotide substitution model. It assumes that the tree requiring the fewest character changes is preferred. This can be problematic when there is homoplasy, unequal rates of evolution, or long-branch attraction. Therefore, I use the parsimony consensus tree as a comparison to the BioNJ distance-based result rather than as a definitive species tree.


```bash
## pwd = scripts/distance_parsimony
mkdir -p ../../results/distance_parsimony/MP
Rscript parsimony_based_tree.R 2>&1 | tee ../../results/distance_parsimony/parsimony/_MP_run.log
```

I used `optim.parsimony()` with Fitch parsimony and SPR rearrangements. Fitch parsimony was appropriate because I did not define a custom nucleotide-change cost matrix, so all nucleotide changes were treated with equal cost. I used SPR rearrangements because they search tree space more broadly than NNI by pruning and regrafting subtrees, which is computationally manageable for this five-taxon dataset.

I also checked the result using `bab()`, an exact branch-and-bound maximum-parsimony search suitable for small numbers of taxa. The exact search produced the same parsimony scores as the optimized trees, indicating that the inferred trees had already reached the best parsimony score for these alignments.

As an additional summary of the per-gene maximum-parsimony trees, I generated a consensus network using `consensusNet()`. This was done before rooting, because consensus networks summarize unrooted splits across trees rather than a single rooted history. In my result, the consensus network is fully tree-like and does not show boxes or reticulations, indicating little or no strongly supported conflict among the per-gene maximum-parsimony trees at the chosen split-frequency threshold. The network consistently groups **Bg with Ee** and **Bb with Pk**, with **Ip** connecting between these two pairs.
>This suggests that the marker genes used here recover a largely consistent phylogenetic signal under the maximum-parsimony criterion, and that the consensus network agrees well with the corresponding consensus tree rather than revealing major disagreement among genes.