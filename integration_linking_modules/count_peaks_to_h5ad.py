#!/usr/bin/env python3
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
import argparse
from sklearn.feature_extraction.text import TfidfTransformer

def run(fname, metadata, peakdata, max_val=4):
    import numpy as np
    import scipy.sparse
    import anndata
    df = pd.read_csv(fname, sep="\t", header=None)
    df.columns = ["barcode", "peak", "count"]
    df = df.loc[metadata.index.get_indexer(df["barcode"].values) >= 0, :]
    df = df.loc[peakdata.index.get_indexer(df["peak"].values) >= 0, :]
    ubc, bc_inv = np.unique(df["barcode"].values, return_inverse=True)
    pk_inv = peakdata.index.get_indexer(df["peak"].values)
    S = scipy.sparse.csr_matrix((df["count"].values,
                                 (bc_inv, pk_inv)), dtype=int,
                                shape=(len(ubc), peakdata.shape[0]))
    S.data[S.data > max_val] = max_val
    return anndata.AnnData(S, obs=metadata.loc[ubc, :], var=peakdata, dtype=np.uint8)

def combine(k, dirname):
    pdata = sc.read(os.path.join(dirname, k))
    cdata = anndata.concat({"Gene Accessibility": adata, "Peaks": pdata}, axis=1, label="feature_types", merge="same")
    sc.pp.log1p(cdata)
    sc.pp.calculate_qc_metrics(cdata, inplace=True, percent_top=[])
    cdata.X = cdata.X.astype(np.float32)
    sc.pp.pca(cdata, n_comps=100)
    oname = os.path.join("pseudobulk", k)
    cdata.write_h5ad(oname, compression="gzip")
    return oname

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, nargs="+")
    ap.add_argument("-m", "--metadata", required=True)
    ap.add_argument("-p", "--peaks", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = vars(ap.parse_args())
    metadata=pd.read_csv(args["metadata"], sep="\t")
    peaks=pd.read_csv(args["peaks"], sep="\t")
    peaks.columns = ["chrom", "begin", "end", "peak_name"]
    peaks.index = peaks["peak_name"].values
    if "obsname" in metadata.columns:
        metadata.index = metadata["obsname"].values
    if "projid" in metadata.columns:
        metadata["projid"] = pd.Categorical(metadata["projid"])
    load = partial(run, metadata=metadata, peakdata=peaks)
    with multiprocessing.Pool() as p:
        out = p.map(load, args["input"])
    adata = anndata.concat(out, merge="same")
    adata.raw = adata.copy()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = TfidfTransformer(smooth_idf=False).fit_transform(adata.raw.X)
    sc.pp.normalize_total(adata, target_sum=10000)
    adata.X = adata.X.astype(np.float32)
    adata.write_h5ad(args["output"], compression="gzip")
