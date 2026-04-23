import pandas as pd
import os
import re

#----------------------------------------
# define functions
#----------------------------------------
# function to get the longest cds sequence per gene from a fasta file
def get_longest_cds_per_gene(fasta_file, species_tag):

	gene_to_best = {}
	current_header = None
	current_seq_lines = []

	with open(fasta_file, "r") as f:
		for line in f:
			line = line.rstrip("\n")

			if line.startswith(">"):
				if current_header is not None:
					seq = "".join(current_seq_lines)

					gene_match = gene_pattern.search(current_header)
					protein_match = protein_pattern.search(current_header)
					protID_match = protID_pattern.search(current_header)

					if gene_match is not None:
						gene = gene_match.group(1).strip().lower()
						protein = protein_match.group(1).strip() if protein_match else "NA"
						protID = protID_match.group(1).strip() if protID_match else "NA"
						lcl = current_header.split()[0].replace(">", "")

						rec = {
							"species": species_tag,
							"gene": gene,
							"protein": protein,
							"protID": protID,
							"lcl": lcl,
							"seq": seq,
							"length": len(seq)
						}

						if (gene not in gene_to_best) or (rec["length"] > gene_to_best[gene]["length"]):
							gene_to_best[gene] = rec

				current_header = line
				current_seq_lines = []

			else:
				current_seq_lines.append(line)

	if current_header is not None:
		seq = "".join(current_seq_lines)

		gene_match = gene_pattern.search(current_header)
		protein_match = protein_pattern.search(current_header)
		protID_match = protID_pattern.search(current_header)

		if gene_match is not None:
			gene = gene_match.group(1).strip().lower()
			protein = protein_match.group(1).strip() if protein_match else "NA"
			protID = protID_match.group(1).strip() if protID_match else "NA"
			lcl = current_header.split()[0].replace(">", "")

			rec = {
				"species": species_tag,
				"gene": gene,
				"protein": protein,
				"protID": protID,
				"lcl": lcl,
				"seq": seq,
				"length": len(seq)
			}

			if (gene not in gene_to_best) or (rec["length"] > gene_to_best[gene]["length"]):
				gene_to_best[gene] = rec

	return gene_to_best


# function to write a fasta record
def write_fasta_record(fout, rec):

	protein_clean = rec["protein"].replace(" ", "_")
	header = f">{rec['species']}|gene={rec['gene']}|protein={protein_clean}|protID={rec['protID']}|lcl={rec['lcl']}"
	fout.write(header + "\n")

	seq = rec["seq"]
	for i in range(0, len(seq), 80):
		fout.write(seq[i:i+80] + "\n")


# define a class to hold species information
class SpeciesInfo:
	def __init__(self, species_tag: str, file_path: str):
		self.species_tag = species_tag
		self.file_path = file_path


electric_species_info = [
	SpeciesInfo(species_tag = "Bb", file_path = "../../data/filtered_per_species_CDS_fa/Bb_intersect_genes_cds.fna"),
	SpeciesInfo(species_tag = "Bg", file_path = "../../data/filtered_per_species_CDS_fa/Bg_intersect_genes_cds.fna"),
	SpeciesInfo(species_tag = "Ee", file_path = "../../data/filtered_per_species_CDS_fa/Ee_intersect_genes_cds.fna"),
	SpeciesInfo(species_tag = "Pk", file_path = "../../data/filtered_per_species_CDS_fa/Pk_intersect_genes_cds.fna"),
	SpeciesInfo(species_tag = "Ip", file_path = "../../data/filtered_per_species_CDS_fa/Ip_intersect_genes_cds.fna")
]

#----------------------------------------
# input / output
#----------------------------------------
marker_file = "../../data/filtered_per_species_CDS_fa/intersect_genes_all_species.txt"
outdir = "../../data/input_to_MSA_perGene_fa/"
os.makedirs(outdir, exist_ok=True)

#----------------------------------------
# regex patterns
#----------------------------------------
gene_pattern = re.compile(r"\[gene=([^\]]+)\]")
protein_pattern = re.compile(r"\[protein=([^\]]+)\]")
protID_pattern = re.compile(r"\[protein_id=([^\]]+)\]")

#----------------------------------------
# get top 25 up and top 25 down markers
#----------------------------------------
markers = pd.read_csv(marker_file, sep="\t", header=0)

top_up = markers[markers["direction"] == "up"].head(25).copy()
top_down = markers[markers["direction"] == "down"].head(25).copy()

final_markers = pd.concat([top_up, top_down], ignore_index=True)
final_genes = final_markers["match_name"].tolist()

#----------------------------------------
# get longest sequence per gene in each species
#----------------------------------------
species_gene_records = {}

for sp in electric_species_info:
	gene_to_best = get_longest_cds_per_gene(sp.file_path, sp.species_tag)
	species_gene_records[sp.species_tag] = gene_to_best
	print(f"{sp.species_tag}: {len(gene_to_best)} genes with longest CDS selected")

#----------------------------------------
# write one fasta file per gene
#----------------------------------------
for gene in final_genes:
	outfile = os.path.join(outdir, f"{gene}.fna")

	with open(outfile, "w") as fout:
		for sp in electric_species_info:
			if gene in species_gene_records[sp.species_tag]:
				rec = species_gene_records[sp.species_tag][gene]
				write_fasta_record(fout, rec)

	print(f"wrote {outfile}")