#!/usr/bin/env python3

import scanpy as sc
import numpy as np
import anndata
import pandas as pd
import argparse

def combine(pdata, gdata="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/TSS6/GeneMatrix_estimated_0.h5ad", output="combined.h5ad"):
    if type(gdata) == str:
        gdata = sc.read(gdata)
    if type(pdata) == str:
        pdata = sc.read(pdata)
    cdata = anndata.concat({"Gene Accessibility": gdata, "Peaks": pdata}, axis=1, label="feature_types", merge="same")
    del pdata, gdata
    sc.pp.log1p(cdata)
    sc.pp.calculate_qc_metrics(cdata, inplace=True, percent_top=[])
    cdata.X = cdata.X.astype(np.float32)
    sc.pp.pca(cdata, n_comps=100)
    cdata.write_h5ad(output, compression="gzip")
    return output

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--peaks", required=True)
    ap.add_argument("-a", "--accessibility", default="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/TSS6/GeneMatrix_estimated_0.h5ad")
    ap.add_argument("-o", "--output", required=True)
    args = vars(ap.parse_args())
    combine(args["peaks"], gdata=args["accessibility"], output=args["output"])
