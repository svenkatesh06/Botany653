library(ape)

indir <- "../../results/msa_macse/"
# reorder sequences so that all genes have the same order of species
species <- c("Bg", "Bb", "Pk", "Ip", "Ee")

#-------------------------------------
# Format nucleotide alignments
#-------------------------------------
nt_dir <- file.path(indir, "cleaned_nt")
files <- list.files(nt_dir, pattern = "\\.macse\\.cleaned\\.fna$", full.names = TRUE)
outdir <- file.path(indir, "reordered_cleaned_macse_nt")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (f in files) {
	cat("Processing file:", basename(f), "\n")
	x <- read.FASTA(f, type = "DNA")

	# read.FASTA gives a list, so names(x) are the sequence names,
	# and x[[i]] is the sequence for the i-th entry
	names(x) <- sapply(strsplit(names(x), "\\|"), function(z) z[1])

	if (!all(species %in% names(x))) {
		stop("Missing species in file: ", basename(f))
	}

	x <- x[species]

	seq_lens <- sapply(x, length)
	if (length(unique(seq_lens)) != 1) {
		stop("Sequences have different lengths in file: ", basename(f))
	}

	out_file <- file.path(outdir, paste0(sub("\\.macse\\.cleaned\\.fna$", "", basename(f)), ".fna"))
	write.FASTA(x, file = out_file)
}

#-------------------------------------
# Format AA alignments
#-------------------------------------
aa_dir <- file.path(indir, "cleaned_aa")
files <- list.files(aa_dir, pattern = "\\.macse\\.cleaned\\.faa$", full.names = TRUE)
outdir <- file.path(indir, "reordered_cleaned_macse_aa")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (f in files) {
	cat("Processing file:", basename(f), "\n")
	x <- read.FASTA(f, type = "AA")

	# read.FASTA gives a list, so names(x) are the sequence names,
	# and x[[i]] is the sequence for the i-th entry
	names(x) <- sapply(strsplit(names(x), "\\|"), function(z) z[1])

	if (!all(species %in% names(x))) {
		stop("Missing species in file: ", basename(f))
	}

	x <- x[species]

	seq_lens <- sapply(x, length)
	if (length(unique(seq_lens)) != 1) {
		stop("Sequences have different lengths in file: ", basename(f))
	}

	out_file <- file.path(outdir, paste0(sub("\\.macse\\.cleaned\\.faa$", "", basename(f)), ".faa"))
	write.FASTA(x, file = out_file)
}