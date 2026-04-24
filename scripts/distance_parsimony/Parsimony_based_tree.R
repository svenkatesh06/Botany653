library(ape)
library(adegenet)
library(phangorn)

indir <- "../../results/msa_macse/cleaned_nt"
files <- list.files(indir, pattern = "\\.macse\\.cleaned\\.fna$", full.names = TRUE)

outdir <- "../../results/distance_parsimony/max_parsimony"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species <- c("Bg", "Bb", "Pk", "Ip", "Ee")
alns <- list()