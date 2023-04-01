#!/usr/bin/env python3
import epiclust as ec
import os
import pandas as pd

adata = anndata.read("TileMatrix.h5ad")
V = ec.gene_distance.df_to_pyranges(ec.gene_distance.peak_names_to_var(adata.var.loc[adata.var["selected"],:].index.values), chrom="seqname", left="start", right="end", name="Tile")
Dname="/home/benjames/data-group/Xushen_Sharing/ToBen/ATAC.matrix/matrix.by.major_celltype.20220816/"
tbl = {}
for x in os.listdir(Dname):
    if os.path.exists(os.path.join(Dname, x, "peak.annotation.txt")):
        df = pd.read_csv(os.path.join(Dname, x, "peak.annotation.txt"), sep="\t")
        df.index = ["%s:%d-%d" % (chrom, begin, end) for chrom, begin, end in zip(df["seqnames"], df["start"], df["end"])]
        tbl[x] = df
        pr = ec.gene_distance.df_to_pyranges(df, chrom="seqnames", left="start", right="end", state="peakType", name="Peak") ### todo: kwargs for item assignment
        pf = V.join(pr, report_overlap=True).df
        tbl[x] = pf.sort_values("Overlap", ascending=False).drop_duplicates("Tile")
        tbl[x]["percent_overlap"] = tbl[x]["Overlap"] / (tbl[x]["End"] - tbl[x]["Start"])
        tbl[x].index = tbl[x]["Tile"].values


for k, df in tbl.items():
    adata.var[k] = ""
    df = df.loc[df["percent_overlap"] > 0.80, :]
    adata.var.loc[df.index.values, k] = df["state"]
        #FL = {x: pd.read_csv(os.path.join(Dname, x, "peak.annotation.txt") for x in os.listdir(Dname) if os.path.exists(os.path.join(Dname, x, "peak.annotation.txt"))]
