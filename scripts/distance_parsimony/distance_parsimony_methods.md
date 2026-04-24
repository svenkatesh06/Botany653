# DISTANCE AND PARSIMONY METHODS

Useful links:

1. https://adegenet.r-forge.r-project.org/files/MSc-intro-phylo.1.1.pdf
2. https://cran.r-project.org/web/packages/ape/ape.pdf

## DISTANCE-BASED TREE

I chose BioNJ rather than WPGMA because WPGMA assumes a clock-like or ultrametric tree, where all taxa are equally distant from the root. I don't believe this assumption will hold for the electric fish data, especially because the taxa come from different lineages and may have unequal rates of sequence evolution. BioNJ is still a distance-based method and therefore has limitations, but like Neighbor-Joining, it does not require equal evolutionary rates across lineages.

| Software | Description | Strengths | Weaknesses | Assumptions | User choices |
| :---: | :--- | :--- | :--- | :--- | :--- |
| `ape` [distance method is BioNJ] | Distance-based phylogenetic method. Aligned DNA sequences are converted into a pairwise genetic distance matrix. Then BioNJ constructs an unrooted tree using a neighbor-joining framework, but incorporates information about the variance of distance estimates when choosing which taxa/clusters to join. | Fast and simple; works well as an initial tree estimate; does not require a strict molecular clock; useful for comparing overall genetic divergence among the five taxa. | Reduces the full alignment to pairwise distances, so site-pattern information is partly lost; tree quality depends strongly on the chosen distance model and alignment quality; does not explicitly optimize likelihood or parsimony; branch lengths/topology can be affected by missing data and gaps. | Homologous sites are correctly aligned; the selected substitution model reasonably describes sequence divergence; distances are sufficiently tree-like/additive; each gene tree reflects the evolutionary signal in that gene; the consensus tree summarizes clades repeatedly recovered across genes. | Choice of input alignment - MACSE NT alignments; DNA distance model; treatment of gaps/missing data; root |


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

For the parsimony-based analysis, I also used the MACSE nucleotide alignments.


