library(ape)
library(adegenet)
library(phangorn)

indir <- "../../results/msa_macse/cleaned_nt"
files <- list.files(indir, pattern = "\\.macse\\.cleaned\\.fna$", full.names = TRUE)

outdir <- "../../results/distance_parsimony/parsimony"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species <- c("Bg", "Bb", "Pk", "Ip", "Ee")

cat("Calculating per-gene maximum parsimony trees\n")
gene_trees <- list()
scores <- data.frame()

## estimate one maximum parsimony tree per gene, ensuring that the order of species is consistent across all files
for (f in files) {

	cat("Processing file:", basename(f), "\n")

	x <- read.dna(f, format = "fasta")
	rownames(x) <- sapply(strsplit(rownames(x), "\\|"), function(z) z[1])

	if (!all(species %in% rownames(x))) {
	stop("Missing species in file: ", basename(f))
	}
	x <- x[species, ]

	# convert the DNA alignment to a phangorn phyDat object for parsimony analysis
	x_phy <- as.phyDat(x)

	# calculate a raw-distance NJ tree as the starting tree for maximum parsimony search
	tre_ini <- nj(dist.dna(x, model = "raw", pairwise.deletion = TRUE))

	# calculate the parsimony score of the starting tree
	start_score <- parsimony(tre_ini, x_phy,method = "fitch")
	# search for a maximum parsimony tree by optimizing the starting tree
	tre <- optim.parsimony(tre_ini, x_phy,method = "fitch", rearrangements = "SPR", trace = 1)
	# calculate the parsimony score of the optimized tree
	final_score <- parsimony(tre, x_phy,method = "fitch")
	# reorganizes the internal structure of the tree to get the ladderized effect when plotted
	tre <- ladderize(tre)

	gene <- sub("\\.macse\\.cleaned\\.fna$", "", basename(f))
	gene_trees[[gene]] <- tre

	scores <- rbind(scores, data.frame( gene = gene, start_score = start_score, final_score = final_score ))
}

class(gene_trees) <- "multiPhylo"

# write all per-gene trees to one file
all_tree_file <- file.path(outdir, "All_genes_macse_MP_trees.tree")

con <- file(all_tree_file, open = "w")
for (gene in names(gene_trees)) {
  writeLines(paste0("# ", gene), con)
  writeLines(write.tree(gene_trees[[gene]]), con)
}
close(con)


# create a consensus network from the per-gene maximum parsimony trees
cnet <- consensusNet(gene_trees, prob = 0.1)
pdf(file.path(outdir, "MP_consensus_network.pdf"), width = 7, height = 6)
plot(cnet, main = "Maximum parsimony consensus network")
dev.off()


# write the parsimony scores for the starting and optimized trees
write.table(scores, file.path(outdir, "All_genes_macse_MP_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# create a majority-rule consensus tree from the per-gene maximum parsimony trees
consensus_tree <- consensus(gene_trees, p = 0.8)

# root the consensus tree using Ip as the outgroup
consensus_tree <- root(consensus_tree, outgroup = "Ip", resolve.root = TRUE)
consensus_tree <- ladderize(consensus_tree)
write.tree(consensus_tree, file.path(outdir, "Consensus_macse_MP.tree"))

# plot the consensus tree and save it as a PDF
pdf(file.path(outdir, "Consensus_macse_MP.pdf"), width = 7, height = 6)
plot(consensus_tree, main = "MACSE maximum parsimony consensus tree")
dev.off()