library(ape)
library(adegenet)
library(phangorn)

indir <- "../../results/msa_macse/reordered_cleaned_nt"
files <- list.files(indir, pattern = "\\.fna$", full.names = TRUE)

outdir <- file.path("../../results/distance_parsimony/BioNJ")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


for (model in c("F84", "TN93")) {
	cat("Calculating per-gene distance matrices and BioNJ trees for model:", model, "\n")
	gene_trees <- list()

	## estimate one BioNJ tree per gene
	for (f in files) {

		cat("Processing file:", basename(f), "\n")
		x <- read.dna(f, format = "fasta")

		# calculate the distance matrix for this gene
		d <- dist.dna(x, model = model, pairwise.deletion = TRUE)
		# calculate the BioNJ tree
		tre <- bionj(d)
		# reorganizes the internal structure of the tree to get the ladderized effect when plotted
		tre <- ladderize(tre)

		gene <- sub("\\.macse\\.cleaned\\.fna$", "", basename(f))
		gene_trees[[gene]] <- tre

	}

	class(gene_trees) <- "multiPhylo"
	# write all per-gene trees to one file
	all_tree_file <- file.path(outdir, paste0("All_genes_", model, "_macse_BioNJ_trees.tree"))
	con <- file(all_tree_file, open = "w")
	for (gene in names(gene_trees)) {
		writeLines(paste0("# ", gene), con)
		writeLines(write.tree(gene_trees[[gene]]), con)
	}
	close(con)

	# create a majority-rule consensus tree from the per-gene BioNJ trees
	consensus_tree <- consensus(gene_trees, p = 0.8)
	# root the consensus tree using Ip as the outgroup
	consensus_tree <- root(consensus_tree, outgroup = "Ip", resolve.root = TRUE)

	consensus_tree <- ladderize(consensus_tree)
	write.tree(consensus_tree, file.path(outdir, paste0("Consensus_", model, "_macse_BioNJ.tree")))

	# plot the consensus tree and save it as a PDF
	pdf(file.path(outdir, paste0("Consensus_", model, "_macse_BioNJ.pdf")), width = 7, height = 6)
	plot(consensus_tree, main = paste("MACSE", model, "BioNJ consensus tree"))
	dev.off()
}