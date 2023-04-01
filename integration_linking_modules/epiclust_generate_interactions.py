#!/usr/bin/env python3
import argparse
from epiclust.gene_distance import distance_weight_all, peak_names_to_var, shuffle_peaks
from epiclust.gtf import load_gtf
import pandas as pd
import anndata

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-n", "--num-shuffle", default=5, type=int)
    args = vars(ap.parse_args())
    var = anndata.read(args["input"], backed="r").var
    dw = distance_weight_all(peak_names_to_var(var.index.values), load_gtf(args["gtf"]))
    dw["interaction_index"] = ["%s|%s" % (peak, gene) for peak, gene in zip(dw["peak"], dw["gene"])]
    dw = dw.loc[dw["gene"].isin(var.index.values), :] ### filter genes to those used
    nw = shuffle_peaks(dw, size=args["num_shuffle"])
    df = pd.concat({1: dw, -1: nw}).reset_index(0).rename({"level_0": "label"}, axis=1)
    df.to_csv(args["output"], sep="\t", index=False)
