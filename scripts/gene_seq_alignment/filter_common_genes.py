import pandas as pd
import os, sys


## load Ee file and mapping.
topdir="../../data/metadata_process_filter"

ee_meta = pd.read_csv(f"{topdir}/Ee_Expr.tsv", sep="\t", header=0, usecols=["match_name", "cluster_ID", "main EO", "Sachs' EO", "Hunter's EO"])
ee_meta.rename(columns={"match_name": "gene"}, inplace=True)
ee_meta["gene"] = ee_meta["gene"].astype(str).str.strip()

ee_map = pd.read_csv("../../data/fasta_seq/e_eel_ncbi_data/GCF_003665695.1/gene_transcript_map.txt", sep="\t", header=None, names=["gene", "transcript_ID"])
ee_map["gene"] = ee_map["gene"].astype(str).str.strip()

# keep only non-empty gene names
ee_meta = ee_meta[ee_meta["cluster_ID"].isin([1, 9])].copy()  # keep only up and down clusters
ee_meta = ee_meta[ee_meta["gene"].notna() & (ee_meta["gene"] != "")].copy()

## merge transcript mapping first, then keep only rows with transcript_ID available
ee_meta = ee_meta.merge(ee_map, on="gene", how="left")
ee_meta = ee_meta[ee_meta["transcript_ID"].notna() & (ee_meta["transcript_ID"].astype(str).str.strip() != "")].copy()

ee_meta_up = ee_meta[ee_meta["cluster_ID"] == 9].copy()
ee_meta_down = ee_meta[ee_meta["cluster_ID"] == 1].copy()

ee_up_genes = set(ee_meta_up["gene"].str.lower())
ee_down_genes = set(ee_meta_down["gene"].str.lower())

species_info = {
    "Me": {"file": "Me_Expr.tsv", "usecols": ["group_ID", "annotation", "EO", "SM"]},
    "Sm": {"file": "Sm_Expr.tsv", "usecols": ["subcomponent_ID", "ZFIN_match", "EO_expression", "muscle_expression"]},
    "Ev": {"file": "Ev_Expr.tsv", "usecols": ["seq_ID", "annotation", "EO", "SM"]},
}

species_dfs = {}

for sp, info in species_info.items():
	df = pd.read_csv(f"{topdir}/{info['file']}", sep="\t", header=0, usecols=info["usecols"])
	df.columns = ["seq_ID", "gene", "EO", "SM"]  # standardize column names for gene, SM, EO
	df["gene"] = df["gene"].astype(str).str.strip()

	df = df[df["gene"].notna() & (df["gene"] != "")].copy()

	species_dfs[sp] = df


# find common genes intersection across all species
common_genes_norm = set(ee_meta["gene"].str.lower())

for sp, df in species_dfs.items():
    common_genes_norm &= set(df["gene"].str.lower())

print(f"Total common genes across Ee, Me, Sm, Ev: {len(common_genes_norm)}")

# write the common genes list to a file
common_genes = ee_meta.loc[ee_meta["gene"].str.lower().isin(common_genes_norm), ["gene"]].drop_duplicates().sort_values("gene")
common_genes.to_csv(f"{topdir}/common_genes.txt", sep="\t", index=False, header=False)

# save the Ee table with only the common genes and the transcript mapping
ee_out_mapped = ee_meta.loc[ee_meta["gene"].str.lower().isin(common_genes_norm), ["transcript_ID", "gene", "cluster_ID", "main EO", "Sachs' EO", "Hunter's EO"]]
ee_out_mapped = ee_out_mapped.drop_duplicates().sort_values(["cluster_ID", "gene", "transcript_ID"])
ee_out_mapped.to_csv(f"{topdir}/Ee_common_mapped.txt", sep="\t", index=False)

## save the Me, Sm, Ev tables as in species dfs but only with the common genes
for sp, df in species_dfs.items():
    out = df.loc[df["gene"].str.lower().isin(common_genes_norm), ["seq_ID", "gene", "EO", "SM"]].copy()
    out["regulation"] = out["gene"].str.lower().map(lambda g: "up" if g in ee_up_genes else ("down" if g in ee_down_genes else pd.NA))
    out = out[["seq_ID", "gene", "EO", "SM", "regulation"]].drop_duplicates().sort_values("gene")
    out.to_csv(os.path.join(topdir, f"{sp}_common.txt"), sep="\t", index=False)

print(f"Saved outputs in: {topdir}")