import pandas as pd
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#----------------------------------------
# define functions
#----------------------------------------

def get_longest_cds_per_gene(fasta_file, species_tag):
    gene_to_best = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description

        gene_match = gene_pattern.search(header)
        protein_match = protein_pattern.search(header)
        protID_match = protID_pattern.search(header)

        if gene_match is not None:
            gene = gene_match.group(1).strip().lower()
            protein = protein_match.group(1).strip() if protein_match else "NA"
            protID = protID_match.group(1).strip() if protID_match else "NA"
            lcl = record.id

            rec = {
                "species": species_tag,
                "gene": gene,
                "protein": protein,
                "protID": protID,
                "lcl": lcl,
                "seq": record.seq,
                "length": len(record.seq)
            }

            if (gene not in gene_to_best) or (rec["length"] > gene_to_best[gene]["length"]):
                gene_to_best[gene] = rec

    return gene_to_best


def make_seqrecord(rec):
    protein_clean = rec["protein"].replace(" ", "_")
    header = f"{rec['species']}|gene={rec['gene']}|protein={protein_clean}|protID={rec['protID']}|lcl={rec['lcl']}"

    return SeqRecord(
        rec["seq"],
        id=header,
        description=""
    )


class SpeciesInfo:
    def __init__(self, species_tag: str, file_path: str):
        self.species_tag = species_tag
        self.file_path = file_path


electric_species_info = [
    SpeciesInfo(species_tag="Bb", file_path="../../data/filtered_per_species_CDS_fa/Bb_intersect_genes_cds.fna"),
    SpeciesInfo(species_tag="Bg", file_path="../../data/filtered_per_species_CDS_fa/Bg_intersect_genes_cds.fna"),
    SpeciesInfo(species_tag="Ee", file_path="../../data/filtered_per_species_CDS_fa/Ee_intersect_genes_cds.fna"),
    SpeciesInfo(species_tag="Pk", file_path="../../data/filtered_per_species_CDS_fa/Pk_intersect_genes_cds.fna"),
    SpeciesInfo(species_tag="Ip", file_path="../../data/filtered_per_species_CDS_fa/Ip_intersect_genes_cds.fna")
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
final_genes = final_markers["match_name"].str.lower().tolist()

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
    records_to_write = []

    for sp in electric_species_info:
        if gene in species_gene_records[sp.species_tag]:
            rec = species_gene_records[sp.species_tag][gene]
            records_to_write.append(make_seqrecord(rec))

    SeqIO.write(records_to_write, outfile, "fasta")
    print(f"wrote {outfile}")