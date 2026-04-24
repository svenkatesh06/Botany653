library(ape)
library(adegenet)
library(phangorn)

indir <- "../../results/msa_macse/cleaned_nt"
files <- list.files(indir, pattern = "\\.macse\\.cleaned\\.fna$", full.names = TRUE)

outdir <- "../../results/distance_parsimony/NJ"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species <- c("Bg", "Bb", "Pk", "Ip", "Ee")
alns <- list()
## concatenate the alignments for the specified species across all genes, ensuring that the order of species is consistent across all files
for (f in files) {

	cat("Processing file:", basename(f), "\n")
	x <- read.dna(f, format = "fasta")
	rownames(x) <- sapply(strsplit(rownames(x), "\\|"),'[', 1)

	if (!all(species %in% rownames(x))) {
		stop("Missing species in file: ", basename(f))
	}
	x <- x[species, ]
	alns[[basename(f)]] <- x
}

dna <- do.call(cbind, alns)
write.dna(dna, file = file.path(outdir, "concatenated_macse_alignment.fasta"), format = "fasta")

for (model in c("F84", "TN93")) {
	cat("Calculating distance matrix for model:", model, "\n")
	d <- dist.dna(dna, model = model, pairwise.deletion = TRUE)

	# calculate the neighbor-joining tree
	tre <- bionj(d)

	# root the tree using the outgroup (Ip)
	tre <- root(tre, outgroup = "Ip", resolve.root = TRUE)

	# reorganizes the internal structure of the tree to get the ladderized effect when plotted
	tre <- ladderize(tre)
	# write the distance matrix and tree to files
	write.csv(as.matrix(d), file.path(outdir, paste0(model, "_macse_distance_matrix.csv")))
	write.tree(tre, file.path(outdir, paste0(model, "_macse_NJ.tree")))

	# plot the tree and save it as a PDF
	pdf(file.path(outdir, paste0("macse_", model, "_NJ.pdf")), width = 7, height = 6)
	plot(tre, main = paste("MACSE", model, "Neighbor-Joining tree"))
	dev.off()
}
