#!/usr/bin/env python3

import anndata
import scipy.sparse
import numpy as np
import pandas as pd

def find_pearson_cor(genes, peaks, adata, batch_size=500, eps=1e-16):
    from tqdm.auto import tqdm
    assert len(genes) == len(peaks)
    G = adata.var.index.get_indexer(genes)
    P = adata.var.index.get_indexer(peaks)
    out = np.zeros(len(genes))
    for begin in tqdm(np.arange(0, len(genes), batch_size)):
        end = min(begin + batch_size, len(genes))
        xg = adata.X[:, G[begin:end]]
        xp = adata.X[:, P[begin:end]]
        cor = np.multiply(xg, xp).sum(0)
        out[begin:end] = np.clip(cor, a_min=-1+eps, a_max=1-eps)
    return np.arctanh(out)

def find_mean_std_cor(feat_on, feat_over, adata, eps=1e-16, batch_size=500, frac=0.01, min_std=0.01):
    from tqdm.auto import tqdm
    G = adata.var.index.get_indexer(feat_on)
    G = np.unique(G)
    P = adata.var.index.get_indexer(feat_over)
    P = np.unique(P)
    m_out = np.zeros(len(G))
    s_out = np.zeros(len(G))
    for begin in tqdm(np.arange(0, len(G), batch_size)):
        end = min(begin + batch_size, len(G))
        idx_p = np.random.choice(P, size=int(frac*len(P)), replace=False)
        allcor = adata.X[:, G[begin:end]].T @ adata.X[:, idx_p]
        allcor = np.arctanh(np.clip(allcor, a_min=-1+eps, a_max=1-eps))
        m_out[begin:end] = np.mean(allcor, axis=1)
        s_out[begin:end] = np.std(allcor, axis=1).clip(min=min_std, max=np.inf)
    return pd.DataFrame({"mean": m_out, "std": s_out}, index=adata.var.index.values[G])

def compute_cor(pgdf, pdata_list, frac=0.01, min_std=0.01):
    """todo: 1. shuffle pgdf
     2. for each pdata
         a. find mean,std
         b. find cor for pos & neg and adjust
     3. NGBoost/logistic for distance bins
     4.
    then find_pearson_cor for each"""
    cor_p = np.zeros((pgdf.shape[0], len(pdata_list)))
    cor_g = np.zeros((pgdf.shape[0], len(pdata_list)))
    # pgdf = pgdf.loc[pgdf["gene"].isin(pdata_list.var.index.values), :]
    # pgdf = pgdf.loc[pgdf["peak"].isin(pdata_list.var.index.values), :]
    for i, pdata in enumerate(pdata_list):
        print("[%d] Finding mean, std for peak->gene" % i)
        ctrl_p = find_mean_std_cor(pd.unique(pgdf["peak"]), pd.unique(pgdf["gene"]), pdata, frac=frac)
        print("[%d] Finding mean, std for gene->peak" % i)
        ctrl_g = find_mean_std_cor(pd.unique(pgdf["gene"]), pd.unique(pgdf["peak"]), pdata, frac=frac)
        print("[%d] Computing Pearson cor" % i)
        X = find_pearson_cor(genes=pgdf["gene"].values, peaks=pgdf["peak"].values, adata=pdata)
        diff = (X - ctrl_p.loc[pgdf["peak"].values, "mean"])
        cor_p[:, i] = diff / ctrl_p.loc[pgdf["peak"].values, "std"]
        diff = (X - ctrl_g.loc[pgdf["gene"].values, "mean"])
        cor_g[:, i] = diff / ctrl_g.loc[pgdf["gene"].values, "std"]
    return np.hstack((cor_p, cor_g))

def gen_train_test(pgdf, pdata_list, frac=0.01, min_std=0.01):
    ### 1. select to valid genes only
    all_genes = pgdf["gene"].values
    all_peaks = pgdf["peak"].values
    for pdata in pdata_list:
        all_genes = np.intersect1d(all_genes, pdata.var.index.values)
        all_peaks = np.intersect1d(all_peaks, pdata.var.index.values)
    pgdf = pgdf.loc[pgdf["gene"].isin(all_genes) & pgdf["peak"].isin(all_peaks), :]
    ug, ginv = np.unique(pgdf["gene"].values, return_inverse=True)
    up, pinv = np.unique(pgdf["peak"].values, return_inverse=True)
    rg = ug.copy()
    rp = up.copy()
    print("Shuffling peaks and genes")
    np.random.shuffle(rg)
    np.random.shuffle(rp)
    X_pos = compute_cor(pgdf, frac=frac, min_std=min_std, pdata_list=pdata_list)
    X_neg_rp = compute_cor(pd.DataFrame({"gene": ug[ginv], "peak": rp[pinv]}),
                           frac=frac, min_std=min_std, pdata_list=pdata_list)
    X_neg_rg = compute_cor(pd.DataFrame({"gene": rg[ginv], "peak": up[pinv]}),
                           frac=frac, min_std=min_std, pdata_list=pdata_list)
    print("Combining data")
    X = np.vstack((X_pos, X_neg_rp, X_neg_rg))
    y = np.hstack((np.ones(X_pos.shape[0]),
                   np.zeros(X_neg_rp.shape[0] + X_neg_rg.shape[0])))
    return pd.concat([pgdf, pgdf, pgdf]), X, y.astype(int)

from sklearn.linear_model import LogisticRegressionCV
def predict_links(pgdf, X, y, bin_size=1000, train_bins=5):
    dbin = np.round(pgdf["distance"].values / bin_size).astype(int)
    udbin, inv = np.unique(dbin, return_inverse=True)
    out = np.zeros(X.shape[0])
    for i, dbin in enumerate(udbin):
        print("distance:", dbin*bin_size)
        model = LogisticRegressionCV(n_jobs=-1, class_weight="balanced")
        train_flag = np.abs(i - inv) < train_bins
        model.fit(X[train_flag, :], y[train_flag])
        out[i == inv] = model.predict_proba(X[i == inv, :])[:, 1]
    pgdf["prob"] = out
    return pgdf.loc[y > 0, :].sort_values("prob", ascending=False)

import argparse
from gtf import load_gtf
from gene_distance import distance_weight_all, peak_names_to_var

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--pseudobulk", nargs="+")
    ap.add_argument("-g", "--gtf", required=True)
    ap.add_argument("-d", "--max-distance", type=int, default=1000*1000)
    ap.add_argument("-o", "--output", required=True)
    args = vars(ap.parse_args())
    print("Reading pseudobulk")
    AL = [anndata.read(x) for x in args["pseudobulk"]]
    peaks = peak_names_to_var(AL[0].var.loc[AL[0].var["feature_types"] == "Peaks", :].index.values)
    print("Reading GTF")
    gtf = load_gtf(args["gtf"])
    print("Computing distances")
    dw = distance_weight_all(peaks, gtf, max_distance=args["max_distance"])
    pgdf, X, y = gen_train_test(dw, AL)
    del AL
    links = predict_links(pgdf, X, y)
    del pgdf
    links.to_csv(args["output"], sep="\t", index=False)
