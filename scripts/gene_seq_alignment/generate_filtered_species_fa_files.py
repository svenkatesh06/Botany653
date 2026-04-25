import pandas as pd
import os
import re
#----------------------------------------
# define functions
#----------------------------------------
# function to get the set of genes present in a fasta file, filtering by an allowed set of genes
def get_genes_in_fasta(fasta_file, allowed_genes=None):

	found = set()
	for rec in SeqIO.parse(fasta_file, "fasta"):
		 m = gene_pattern.search(rec.description)
		if m:
			gene = m.group(1).strip()
			gene_cmp = gene.lower()

			if (allowed_genes is None) or (gene_cmp in allowed_genes):
				found.add(gene_cmp)

	return found

# function to write out the filtered fasta file for a species
def write_filtered_fasta(fasta_file, out_file, keep_genes):
	kept_records = 0
	with open(out_file, "w") as fout:
		for rec in SeqIO.parse(fasta_file, "fasta"):
			m = gene_pattern.search(rec.description
			if m:
				gene = m.group(1).strip().lower(
				if gene in keep_genes:
					SeqIO.write(rec, fout, "fasta")
					kept_records +=
	return kept_records

# define a class to hold species information
class SpeciesInfo:
	def __init__(self, species_tag: str, file_path: str):
		self.species_tag = species_tag
		self.file_path = file_path

electric_species_info = [
	SpeciesInfo(species_tag = "Bb", file_path = "Bb_assembly_BBRACH_0.4/ncbi_dataset/data/GCF_023856365.1/cds_from_genomic.fna"),
	SpeciesInfo(species_tag = "Bg", file_path = "Bg_assembly_BGAUD_0.2/ncbi_dataset/data/GCF_052324685.1/cds_from_genomic.fna"),
	SpeciesInfo(species_tag = "Ee", file_path = "Ee_genome_ASM4190279v1_2024/ncbi_dataset/data/GCF_041902795.1/cds_from_genomic.fna"),
	SpeciesInfo(species_tag = "Pk", file_path = "Pk_assembly_PKINGS_0.4/ncbi_dataset/data/GCF_048594095.1/cds_from_genomic.fna"),
	SpeciesInfo(species_tag = "Ip", file_path = "Ip_assembly_Coco_2.0/ncbi_dataset/data/GCF_001660625.3/cds_from_genomic.fna")
]

datadir= "../../data/ncbi_genome_downloads/"
outdir = "../../data/filtered_per_species_CDS_fa/"
os.makedirs(outdir, exist_ok=True)

## load the reference markers obtained from electric eel
geneset_up = pd.read_csv("../../data/Ee_filtered_markers_gallant_etal2014/Ee_cluster9_top100.tsv", sep="\t",header=0)
geneset_down = pd.read_csv("../../data/Ee_filtered_markers_gallant_etal2014/Ee_cluster1_top100.tsv", sep="\t",header=0)

geneset_up = geneset_up[["match_name"]].copy()
geneset_up["direction"] = "up"

geneset_down = geneset_down[["match_name"]].copy()
geneset_down["direction"] = "down"

ref_markers_df = pd.concat([geneset_up, geneset_down], ignore_index=True)
ref_markers_df = ref_markers_df.drop_duplicates(subset=["match_name"], keep="first").copy()
ref_markers_df["match_name_cmp"] = ref_markers_df["match_name"].str.lower()

ref_markers = set(ref_markers_df["match_name_cmp"].tolist())

## pattern to extract the gene name from fasta header
gene_pattern = re.compile(r"\[gene=([^\]]+)\]")

#----------------------------------------
# find intersection markers across species
#----------------------------------------
species_gene_sets = {}

for sp in electric_species_info:
	fasta_path = os.path.join(datadir, sp.file_path)
	genes_present = get_genes_in_fasta(fasta_path, allowed_genes = ref_markers)
	species_gene_sets[sp.species_tag] = genes_present
	print(f"{sp.species_tag}: {len(genes_present)} EO markers found")


intersect_genes = set.intersection(*species_gene_sets.values())
print(f"\nGenes present in ALL species: {len(intersect_genes)}")

intersect_df = ref_markers_df[ref_markers_df["match_name_cmp"].isin(intersect_genes)].copy()
intersect_df["match_name"] = intersect_df["match_name"].str.lower()
intersect_df = intersect_df.drop(columns=["match_name_cmp"])
intersect_df.to_csv(os.path.join(outdir, "intersect_genes_all_species.txt"), sep="\t", index=False )


# write out filtered fasta files for each species
for sp in electric_species_info:
	fasta_path = os.path.join(datadir, sp.file_path)
	out_fa = os.path.join(outdir, f"{sp.species_tag}_intersect_genes_cds.fna")
	nrec = write_filtered_fasta(fasta_path, out_fa, intersect_genes)
	print(f"{sp.species_tag}: wrote {nrec} CDS records -> {out_fa}")