# DISTANCE AND PARSIMONY METHODS

Useful links:

1. https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf
2. https://cran.r-project.org/web/packages/ape/ape.pdf

## DISTANCE-BASED TREE

I chose BioNJ rather than WPGMA because WPGMA assumes a clock-like or ultrametric tree, where all taxa are equally distant from the root. I don't believe this assumption will hold for the electric fish data, especially because the taxa come from different lineages and may have unequal rates of sequence evolution. BioNJ is still a distance-based method and therefore has limitations, but like Neighbor-Joining, it does not require equal evolutionary rates across lineages.

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :--- | :--- | :--- | :--- | :--- |
| `ape` [distance method is BioNJ] | Distance-based phylogenetic method. Aligned DNA sequences are converted into a pairwise genetic distance matrix. Then BioNJ constructs an unrooted tree using a neighbor-joining framework, but incorporates information about the variance of distance estimates when choosing which taxa/clusters to join. | Fast and simple; works well as an initial tree estimate; does not require a strict molecular clock; useful for comparing overall genetic divergence among the five taxa. | Reduces the full alignment to pairwise distances, so site-pattern information is partly lost; tree quality depends strongly on the chosen distance model and alignment quality; does not explicitly optimize likelihood or parsimony; branch lengths/topology can be affected by missing data and gaps. | Homologous sites are correctly aligned; the selected substitution model reasonably describes sequence divergence; distances are sufficiently tree-like/additive; concatenated genes are treated as representing a shared evolutionary history. | Choice of input alignment - MACSE NT alignments; DNA distance model; treatment of gaps/missing data; root |


### BioNJ based tree using `ape` and `phangorn`

For the distance-based analysis, I used the MACSE nucleotide alignments because MACSE is designed for coding sequences and accounts for frameshifts and stop codons during alignment. I concatenated the per-gene alignments into a single supermatrix, estimated pairwise DNA distances using `dist.dna()`, and inferred a BioNJ tree using `bionj()`.

> Because this analysis concatenates multiple MACSE-aligned coding sequences into one supermatrix, the resulting BioNJ tree is treated as a species-level phylogenetic estimate for the 5 sampled taxa rather than as an individual gene tree. A gene tree would be inferred separately from one gene alignment, while this analysis combines information across many selected genes. However, this interpretation assumes that the selected genes broadly reflect the same underlying species history. This may not always be true because individual genes can have different evolutionary histories due to incomplete lineage sorting, gene duplication/loss, selection, or missing-data/alignment effects. Therefore, I interpret this result as a distance-based species-level tree estimate from the selected marker genes, not as a definitive species tree.

I tested 2 NT substitution models for estimating the distance matrix: **F84 and TN93**. Both models are more flexible than the simplest distance models because they allow unequal base frequencies. This is useful for real nucleotide sequence data, where A, C, G, and T are often not present at equal frequencies.
- F84 also allows transitions and transversions to occur at different rates.
- TN93 is even more flexible because it allows unequal base frequencies and also allows the two classes of transitions, A↔G and C↔T, to have different rates.

Finally, I used TN93 as the main distance model and F84 as a comparison to check whether the inferred BioNJ topology was sensitive to the choice of distance model.

In `dist.dna()`, I used `pairwise.deletion = TRUE`. Although the input is a multiple sequence alignment, the MACSE alignments can still contain gaps, ambiguous characters, or missing positions for some taxa. Pairwise deletion means that, for each pair of taxa, the distance is calculated using only the sites that are available and comparable for that particular pair. This allows more of the alignment to contribute to the distance estimates instead of removing every column that contains a gap or missing value in any taxon. This is useful for a concatenated multi-gene alignment, where a small number of problematic sites in one taxon would otherwise cause those sites to be excluded for all taxa.

BioNJ estimates an unrooted tree from the distance matrix. To make the tree biologically interpretable, I rooted the BioNJ tree using *Ictalurus punctatus* as the outgroup because it is the non-electric fish in this dataset. Rooting with the outgroup places the root on the branch separating *I.punctatus* from the electric fish taxa and provides directionality for interpreting relationships among the electric lineages. After rooting, I ladderized the tree to make the plotted topology easier to read. Ladderizing only changes the visual ordering of branches; it does not change the inferred relationships or branch lengths.

```bash
## pwd = scripts/distance_parsimony
Rscript distance_based_tree.R 2>&1 | tee ../../results/distance_parsimony/NJ/NJ_run.log
```

**The F84 and TN93 analyses produced the same overall tree topology.** In both trees, Ee grouped with Bg, and Bb grouped with Pk, while Ip was separated as the outgroup. The branch lengths and pairwise distances were also very similar between models. For example, the Bb–Pk distance was 0.06238 under F84 and 0.06242 under TN93, while the Bg–Ee distance was 0.14149 under F84 and 0.14157 under TN93. This indicates that the distance-based result is not strongly sensitive to whether F84 or TN93 is used.

**The recovered topology is biologically reasonable.** Ee and Bg are both gymnotiform electric fishes, while Bb and Pk are both mormyrid electric fishes. Ip, the channel catfish, was used as the non-electric outgroup to root the tree. The BioNJ tree correctly grouped the two gymnotiforms together and the two mormyrids together, which is consistent with their known evolutionary relationships.

The distance-based method was able to capture the overall genetic divergence patterns among the taxa, even though it does not explicitly model sequence evolution or optimize likelihood/parsimony.

## PARSIMONY-BASED TREE

For the parsimony-based analysis, I also used the MACSE nucleotide alignments. I concatenated the per-gene alignments into a single supermatrix and inferred a parsimony tree using `pratchet()` from the `phangorn` package. The parsimony method seeks the tree topology that minimizes the total number of character state changes (substitutions) required to explain the observed sequence data.

