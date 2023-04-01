#!/usr/bin/env python3
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import bbknn
import argparse

def integrate_atac_rna_prepare(atac, rna, hvg=10000, rna_plus=None, atac_plus=None):
    adata = anndata.concat({"ATAC": atac, "RNA": rna}, label="assay", merge="same")
    if rna_plus:
        for col in rna_plus:
            if rna.obs[col].dtype.kind == "O":
                adata.obs[col] = ""
            else:
                adata.obs[col] = np.nan
            adata.obs.loc[rna.obs.index.values, col] = rna.obs[col].values
    if atac_plus:
        for col in atac_plus:
            if atac.obs[col].dtype.kind == "O":
                adata.obs[col] = ""
            else:
                adata.obs[col] = np.nan
            adata.obs.loc[atac.obs.index.values, col] = atac.obs[col].values
    del atac, rna
    sc.pp.normalize_total(adata, exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    if hvg > 0:
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg, batch_key="assay")
        adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.combat(adata, "assay")
    return adata

def integrate_atac_rna_full(adata, n_comps=50):
    sc.pp.pca(adata, n_comps=n_comps)
    sc.external.pp.harmony_integrate(adata, "assay", max_iter_harmony=50)
    sc.external.pp.bbknn(adata, "assay", use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    return adata

def read_anndata(filename, subset=[], celltype_from="CellType", celltype_to="major_celltype"):
    adata = anndata.read(filename, backed="r")
    if celltype_from in adata.obs.columns:
        adata.obs[celltype_to] = adata.obs[celltype_from]
    flag = np.repeat(True, adata.shape[0])
    if subset is not None:
        for ss in subset:
            S = ss.split("=")
            k = S[0]
            if k in adata.obs and len(S) == 2:
                v = S[1]
                print("Filtering", k, "to", v)
                flag &= np.isin(adata.obs[k].values, v.split(","))
    return adata[adata.obs.index.values[flag], :]

def transfer_labels(adata, label="cell_type_high_resolution", batch_from="RNA", batch="assay", key_added="pred"):
    import scipy.sparse
    import scipy.stats
    flag = adata.obs[batch] == batch_from
    ulab, lab_inv = np.unique(adata.obs.loc[flag, label], return_inverse=True)
    S = scipy.sparse.csr_matrix((np.ones_like(lab_inv), (np.arange(len(lab_inv)), lab_inv)))
    X = adata.obsp["connectivities"][:, flag].dot(S).todense()
    Z = (X.max(1) - X.mean(1))/(X.std(1).clip(1e-5, np.inf) * np.sqrt(X.shape[1]))
    adata.obs[key_added] = pd.Categorical(ulab[np.ravel(X.argmax(1))])
    adata.obs["%s_prob" % key_added] = scipy.stats.norm.cdf(Z)
    if "%s_colors" % label in adata.uns:
        kv = {k: v for k, v in zip(adata.obs[label].cat.categories, adata.uns["%s_colors" % label])}
        adata.uns["%s_colors" % key_added] = [kv[cat] for cat in adata.obs[key_added].cat.categories]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--atac", required=True)
    ap.add_argument("--rna", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-t", "--tsv", default="output.tsv.gz")
    ap.add_argument("--hvg", type=int, default=10000)
    ap.add_argument("--n-comps", type=int, default=50)
    ap.add_argument("--bbknn", action="store_true", dest="bbknn")
    ap.add_argument("--no-bbknn", action="store_false", dest="bbknn")
    ap.add_argument("--diffmap", action="store_true", dest="diffmap")
    ap.add_argument("--no-diffmap", action="store_false", dest="diffmap")
    ap.add_argument("--atac-celltype", default="Celltype4")
    ap.add_argument("--rna-celltype", default="major.celltype")
    ap.add_argument("--rna-plus", nargs="+", default=["cell_type_high_resolution"])
    ap.add_argument("--atac-plus", nargs="+")
    ap.add_argument("--subset-atac", nargs="+")
    ap.add_argument("--subset-rna", nargs="+")
    ap.add_argument("--leiden", type=int, nargs="+")
    ap.add_argument("--use-celltype", dest="use_celltype", action="store_true")
    ap.set_defaults(bbknn=True, diffmap=True, use_celltype=False)
    args = vars(ap.parse_args())
    atac = read_anndata(args["atac"], subset=args["subset_atac"], celltype_from=args["atac_celltype"], celltype_to="CellType")
    print("ATAC:", atac.shape)
    rna = read_anndata(args["rna"], subset=args["subset_rna"], celltype_from=args["rna_celltype"], celltype_to="CellType")
    print("RNA:", rna.shape)
    adata = integrate_atac_rna_prepare(atac=atac.to_memory(),
                                       rna=rna.to_memory(),
                                       hvg=args["hvg"], rna_plus=args["rna_plus"], atac_plus=args["atac_plus"])
    adata = integrate_atac_rna_full(adata, n_comps=args["n_comps"])
    if args["diffmap"]:
        print("DiffMap")
        sc.tl.diffmap(adata, n_comps=args["n_comps"])
        sc.external.pp.harmony_integrate(adata, "assay", basis="X_diffmap", max_iter_harmony=50)
        sc.external.pp.bbknn(adata, "assay", use_rep="X_pca_harmony")
        sc.tl.umap(adata)
    if args["leiden"] is None:
        print("No clusters prescribed. Adding")
        args["leiden"] = [1.]
    for cls in args["leiden"]:
        sc.tl.leiden(adata, resolution=cls, key_added="leiden_%d" % cls)
    if args["bbknn"]:
        print("BBKNN regression")
        ckey = ["leiden_%d" % k for k in args["leiden"]]
        if args["use_celltype"] and len(pd.unique(adata.obs["CellType"])) > 1:
            ckey.append("CellType")
        bbknn.ridge_regression(adata, batch_key=["assay"], confounder_key=ckey)
        adata = integrate_atac_rna_full(adata, n_comps=args["n_comps"])
        for cls in args["leiden"]:
            sc.tl.leiden(adata, resolution=cls, key_added="leiden_%d" % cls)
    sc.tl.embedding_density(adata, groupby="assay")
    for label in args["rna_plus"]:
        transfer_labels(adata, label=label, batch_from="RNA", batch="assay", key_added="pred_%s" % label)
    adata.obs.to_csv(args["tsv"], sep="\t")
    adata.write_h5ad(args["output"], compression="gzip")
