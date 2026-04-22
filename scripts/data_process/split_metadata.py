import pandas as pd
import re
import os,sys

topdir="../../data/gallant_etal2014_metadata/"
xlsx = f"{topdir}/NIHMS884961-supplement-Table_S1.xlsx"
outdir = f"{topdir}/split_metadata"
os.makedirs(outdir, exist_ok=True)

xls = pd.ExcelFile(xlsx)

for sheet in xls.sheet_names:
    safe = re.sub(r'[\\/*?:"<>|]', "_", sheet)
    df = pd.read_excel(xlsx, sheet_name=sheet)
    df.to_csv(f"{outdir}/{safe}.tsv", sep="\t", index=False)

print(f"Done: {len(xls.sheet_names)} sheets written to {outdir}")