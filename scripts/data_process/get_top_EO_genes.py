import pandas as pd
import numpy as np
import os

# input
topdir = "../../data/gallant_etal2014_metadata/split_metadata"
infile = f"{topdir}/Ee_Expr.tsv"
os.makedirs(topdir, exist_ok=True)
# output
out_up = f"{topdir}/Ee_cluster9_top25.tsv"
out_down = f"{topdir}/Ee_cluster1_top25.tsv"
out_combined = f"{topdir}/Ee_top50_cluster1_9.tsv"

# read
df = pd.read_csv(infile, sep="\t")

# columns
gene_col = "match_name"
cluster_col = "cluster_ID"
muscle_col = "sk. muscle"
eo_cols = ["main EO", "Sachs' EO", "Hunter's EO"]

# keep needed columns and non-empty gene names
keep_cols = ["gene_ID", "match_ID", gene_col, "match_class", cluster_col, muscle_col] + eo_cols
df = df[keep_cols].copy()
df = df[df[gene_col].notna()]
df = df[df[gene_col].astype(str).str.strip() != ""]

# numeric cleanup
for c in [cluster_col, muscle_col] + eo_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# mean EO across the 3 EO tissues
df["EO_mean"] = df[eo_cols].mean(axis=1)

# score: positive means EO > muscle, negative means muscle > EO
df["log2_EO_vs_muscle"] = np.log2((df["EO_mean"] + 1) / (df[muscle_col] + 1))

#keep only clusters 1 and 9
df = df[df[cluster_col].isin([1, 9, 1.0, 9.0])].copy()

# cluster 9: strongest EO-up genes -----
c9 = df[df[cluster_col].isin([9, 9.0])].copy()

# for duplicate gene names, keep the row with highest EO-vs-muscle score
c9 = c9.sort_values("log2_EO_vs_muscle", ascending=False)
c9 = c9.drop_duplicates(subset=[gene_col], keep="first")

top25_up = c9.head(25).copy()

# cluster 1: strongest EO-down / muscle-up genes -----
c1 = df[df[cluster_col].isin([1, 1.0])].copy()

# for duplicate gene names, keep the row with lowest EO-vs-muscle score
c1 = c1.sort_values("log2_EO_vs_muscle", ascending=True)
c1 = c1.drop_duplicates(subset=[gene_col], keep="first")

top25_down = c1.head(25).copy()

# save
top25_up.to_csv(out_up, sep="\t", index=False)
top25_down.to_csv(out_down, sep="\t", index=False)

top50 = pd.concat([top25_up, top25_down], ignore_index=True)
top50.to_csv(out_combined, sep="\t", index=False)

print("cluster 9 selected:", top25_up.shape[0])
print("cluster 1 selected:", top25_down.shape[0])
print("combined selected:", top50.shape[0])