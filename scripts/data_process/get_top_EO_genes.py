import pandas as pd
import numpy as np
import os


# input / output
topdir = "../../data/gallant_etal2014_metadata/split_metadata"
infile = f"{topdir}/Ee_Expr.tsv"

outdir = "../../data/Ee_filtered_markers_gallant_etal2014"
os.makedirs(outdir, exist_ok=True)

# choose how many genes to keep from each cluster
topN = 100

out_up = f"{outdir}/Ee_cluster9_top{topN}.tsv"
out_down = f"{outdir}/Ee_cluster1_top{topN}.tsv"
out_combined = f"{outdir}/Ee_top{2*topN}_cluster1_9.tsv"

# read input
df = pd.read_csv(infile, sep="\t")
# column names
gene_col = "match_name"
cluster_col = "cluster_ID"
muscle_col = "sk. muscle"
eo_cols = ["main EO", "Sachs' EO", "Hunter's EO"]

# ---------------------------
# data cleaup and filtering
# ---------------------------
# keep only needed columns
keep_cols = [gene_col, cluster_col, muscle_col] + eo_cols
df = df[keep_cols].copy()

# keep only non-empty gene names
df = df[df[gene_col].notna()].copy()
df = df[df[gene_col].astype(str).str.strip() != ""].copy()

# drop rows with missing required values
df = df.dropna(subset=[cluster_col, muscle_col] + eo_cols).copy()

# remove gene names containing :, -, (, ), .
bad_pattern = r"[:\-\(\)\.]"
df = df[~df[gene_col].astype(str).str.contains(bad_pattern, regex=True, na=False)].copy()

# numeric cleanup
for c in [cluster_col, muscle_col] + eo_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# keep only clusters 1 and 9
df = df[df[cluster_col].isin([1, 9, 1.0, 9.0])].copy()

# ---------------------------
# compute EO mean and log-ratio
# ---------------------------
df["EO_mean"] = df[eo_cols].mean(axis=1)
df["log2_EO_vs_muscle"] = np.log2((df["EO_mean"] + 1) / (df[muscle_col] + 1))

# ---------------------------
# cluster 9: EO-up
# ---------------------------
c9 = df[df[cluster_col].isin([9, 9.0])].copy()
c9 = c9.sort_values("log2_EO_vs_muscle", ascending=False)
c9 = c9.drop_duplicates(subset=[gene_col], keep="first")
top_up = c9.head(topN).copy()

# ---------------------------
# cluster 1: EO-down / muscle-up
# ---------------------------
c1 = df[df[cluster_col].isin([1, 1.0])].copy()
c1 = c1.sort_values("log2_EO_vs_muscle", ascending=True)
c1 = c1.drop_duplicates(subset=[gene_col], keep="first")
top_down = c1.head(topN).copy()


out_cols = [gene_col, cluster_col, "EO_mean", "log2_EO_vs_muscle"]
top_up = top_up[out_cols].copy()
top_down = top_down[out_cols].copy()

top_combined = pd.concat([top_up, top_down], ignore_index=True)

# sort combined by cluster then score
top_combined = top_combined.sort_values(
    by=[cluster_col, "log2_EO_vs_muscle"],
    ascending=[True, False]
).copy()


top_up.to_csv(out_up, sep="\t", index=False)
top_down.to_csv(out_down, sep="\t", index=False)
top_combined.to_csv(out_combined, sep="\t", index=False)

print("cluster 9 selected:", top_up.shape[0])
print("cluster 1 selected:", top_down.shape[0])
print("combined selected:", top_combined.shape[0])
print("wrote:")
print(out_up)
print(out_down)
print(out_combined)