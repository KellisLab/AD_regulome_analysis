#!/usr/bin/env python3

import argparse
import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
import epiclust as ec

def run_linking(adata, interaction_df, power, covariates=[], squared_correlation=False):
    ec.fit(adata, power=power, batch="feature_types",
           squared_correlation=squared_correlation,
           covariates=list(np.intersect1d(covariates, adata.obs.columns)))
    print("Finding correlations")
    return ec.linking(adata, interaction_df["gene"].values, interaction_df["peak"].values)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("--interactions", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-p", "--power", type=float, default=0.)
    ap.add_argument("--covariates", nargs="+", default=["pmi", "n_genes_by_counts", "log1p_total_counts", "TSSEnrichment", "Library", "projid", "nFrags", "FRIP"])
    ap.add_argument("--squared-correlation", dest="squared_correlation", action="store_true")
    ap.add_argument("--nonsquared-correlation", dest="squared_correlation", action="store_false")
    ap.set_defaults(squared_correlation=False)
    args = vars(ap.parse_args())
    kwargs = args.copy()
    del kwargs["input"], kwargs["output"], kwargs["interactions"]
    out = run_linking(anndata.read(args["input"]),
                      interaction_df=pd.read_csv(args["interactions"], sep="\t"),
                      **kwargs)
    out.to_csv(args["output"], sep="\t", index=False)
