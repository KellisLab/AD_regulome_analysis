#!/usr/bin/env python3
import anndata
import numpy as np
from sklearn.feature_extraction.text import TfidfTransformer
from multiprocessing import Pool
from kneed import KneeLocator
import argparse
from datetime import datetime

def run_nmf(arguments):
    from scopen.MF import NMF
    data, n_components, alpha, max_iter, verbose, random_state, init = arguments
    model = NMF(n_components=n_components,
                random_state=random_state,
                init=init,
                alpha=alpha,
                l1_ratio=0,
                max_iter=max_iter,
                verbose=verbose)
    w_hat = model.fit_transform(X=data)
    h_hat = model.components_
    print(f"ranks: {n_components}, fitting error: {model.reconstruction_err_}")
    return [w_hat, h_hat, model.reconstruction_err_]

def estimate_nmf(data, a_min=2, a_max=30, a_step=1, alpha=1.0, init="nndsvd", random_state=42, verbose=True, max_iter=200, nproc=16):
    arguments_list = list()
    n_components_list = np.arange(a_min, a_max+1, a_step)
    w_hat_dict = {}
    h_hat_dict = {}
    error_list = []
    for n_components in n_components_list:
        arguments = (data, n_components, alpha, max_iter, verbose, random_state, init)
        arguments_list.append(arguments)
    with Pool(processes=nproc) as pool:
        res = pool.map(run_nmf, arguments_list)
    print("res:")
    print(res)
    for i, n_components in enumerate(n_components_list):
            w_hat_dict[n_components] = res[i][0]
            h_hat_dict[n_components] = res[i][1]
            error_list.append(res[i][2])
    if len(error_list) > 1:
        kl = KneeLocator(n_components_list, error_list,
                         S=1.0, curve="convex", direction="decreasing")
        return w_hat_dict[kl.knee], h_hat_dict[kl.knee]
    else:
        return w_hat_dict[n_components_list[0]], h_hat_dict[n_components_list[0]]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--alpha", default=1., type=float)
    ap.add_argument("--min", default=2, type=int)
    ap.add_argument("--max", default=30, type=int)
    ap.add_argument("--step", default=2, type=int)
    ap.add_argument("--random", default=42, type=int)
    ap.add_argument("--max-iter", default=200, type=int)
    ap.add_argument("--nproc", default=16, type=int)
    ap.add_argument("--blacklist", default="blacklist")
    args = vars(ap.parse_args())
    adata = anndata.read(args["input"])
    if args["blacklist"] in adata.var.columns:
        print("Filtering out blacklist", args["blacklist"])
        X = adata[:, ~adata.var[args["blacklist"]]].X
    else:
        X = adata.X
    X = TfidfTransformer(sublinear_tf=True).fit_transform(X)
    W, H = estimate_nmf(X, alpha=args["alpha"],
                        a_min=args["min"], a_max=args["max"],
                        a_step=args["step"], random_state=args["random"],
                        max_iter=args["max_iter"], nproc=args["nproc"])
    adata.obsm["NMF_W"] = W
    if args["filter_peaks"] in adata.var.columns:
        adata.varm["NMF_H"] = np.zeros((adata.shape[1], H.shape[0]), dtype=H.dtype)
        adata.varm["NMF_H"][adata.var[args["filter_peaks"]].values, :] = H.T
    else:
        adata.varm["NMF_H"] = H.T
    adata.write_h5ad(args["output"], compression="gzip")
