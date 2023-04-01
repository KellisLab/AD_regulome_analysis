#!/usr/bin/env python3
import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import argparse
import scanpy as sc
def pseudobulk(integ_df, adata, col=["leiden"], obsm="NMF_W", varm="NMF_H"):
    I = np.intersect1d(integ_df.index.values, adata.obs.index.values)
    adata = adata[I, :]
    integ_df = integ_df.loc[I, :]
    cls_idx = integ_df.groupby(col).ngroup()
    ucls, cls_idx, cls_inv, cls_cnt = np.unique(cls_idx, return_index=True, return_inverse=True, return_counts=True)
    S = scipy.sparse.csr_matrix((1/cls_cnt[cls_inv],
                                 (cls_inv, np.arange(len(cls_inv)))),
                                shape=(len(ucls), len(I)))
    if obsm in adata.obsm.keys() and varm in adata.varm.keys():
        X = S.dot(adata.obsm[obsm]).dot(adata.varm[varm].T)
    elif adata.isbacked:
        X = S.dot(adata.to_memory().X)
    else:
        X = S.dot(adata.X)
    obs = pd.DataFrame(index=ucls)
    for x in np.setdiff1d(integ_df.columns, col):
        for i, cls in enumerate(ucls):
            allval = integ_df[x].values[i == cls_inv]
            if not np.all(allval == allval[0]):
                break
        else:
            obs[x] = integ_df[x].values[cls_idx]
    return anndata.AnnData(X, dtype=np.float32,
                           obs=obs,
                           var=adata.var)

def extract_pca(adata, n_pcs=None):
        """we know that the S component is always positive
        so it can be recontructed from s^2/(DoF)"""
        Us = adata.obsm["X_pca"]
        s = adata.uns["pca"]["variance"]
        s = np.sqrt(s * (adata.shape[0] - 1)).clip(1e-20, np.inf)
        U = Us @ np.diag(1/s)
        VT = adata.varm["PCs"].T
        if n_pcs is not None and n_pcs <= len(s):
                U = U[:, range(n_pcs)]
                s = s[range(n_pcs)]
                VT = VT[range(n_pcs), :]
        return U, s, VT

def read_integration(tsv_list, col):
    L = {f: pd.read_csv(f, sep="\t", index_col=0) for f in tsv_list}
    df = pd.concat(L).reset_index(level=0)
    for x in col:
        df[x] = df["level_0"].astype(str) + "_" + df[x].astype(str)
    return df

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--atac", required=True)
    ap.add_argument("--rna", required=True)
    ap.add_argument("-i", "--integration", required=True, nargs="+")
    ap.add_argument("-c", "--col", required=True, nargs="+")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--rank", dest="rank", action="store_true")
    ap.add_argument("--no-rank", dest="rank", action="store_false")
    ap.add_argument("--filter-atac", default="selected")
    ap.add_argument("-p", "--power", default=1., type=float)
    ap.add_argument("--obsm", default="NMF_W")
    ap.add_argument("--varm", default="NMF_H")
    ap.add_argument("--n-comps", default=50, type=int)
    ap.set_defaults(rank=False)
    args = vars(ap.parse_args())
    df = read_integration(args["integration"], col=args["col"])
    atac = pseudobulk(df,
                      adata=anndata.read(args["atac"], backed="r"),
                      col=args["col"], obsm=args["obsm"], varm=args["varm"])
    if args["filter_atac"] in atac.var.columns:
        print("Filtering ATAC to \"%s\"" % args["filter_atac"])
        atac = atac[:, atac.var[args["filter_atac"]]].copy()
    rna = pseudobulk(df,
                     adata=anndata.read(args["rna"], backed="r"),
                     col=args["col"], obsm=args["obsm"], varm=args["varm"])
    adata = anndata.concat({"Peaks": atac, "Gene Expression": rna},
                           label="feature_types",
                           merge="same",
                           axis=1)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    if scipy.sparse.issparse(adata.X):
        adata.X = np.asarray(adata.X.todense())
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    if np.abs(args["power"] - 1.0) >= 1e-10:
        print("Running PCA")
        sc.pp.pca(adata, n_comps=min(args["n_comps"], np.min(adata.shape)-1))
        U, s, VT = extract_pca(adata)
        adata.X = U @ np.diag(s**args["power"]) @ VT
        del adata.obsm["X_pca"]
        del adata.uns["pca"]
        del adata.varm["PCs"]
    if args["rank"]:
        print("Ranking data")
        adata.X = scipy.stats.rankdata(adata.X, axis=0)
    adata.X = adata.X - np.ravel(adata.X.mean(0))[None, :] ### subtract mean
    norm = np.linalg.norm(adata.X, ord=2, axis=0) ### normalize
    adata.X = np.divide(adata.X, norm[None, :],
                        where=norm[None, :] != 0, out=np.zeros_like(adata.X))
    ### then divide by norm so that outer product is Spearman correlation
    adata.write_h5ad(args["output"], compression="gzip")
