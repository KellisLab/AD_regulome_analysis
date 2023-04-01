
def agg2df(adata, ind="Sample"):
    import numpy as np
    import scipy.sparse
    import pandas as pd
    us, sinv, scnt = np.unique(adata.obs[ind].values, return_inverse=True, return_counts=True)
    S = scipy.sparse.csr_matrix((np.ones(len(sinv)), (sinv, np.arange(len(sinv)))), dtype=np.uint64)
    X = S.dot(scipy.sparse.csr_matrix(adata.layers["raw"])).tocoo()
    md = adata.obs.drop_duplicates(ind).loc[:, [ind, "Pathology", "msex", "pmi", "projid", "braaksc", "cogdx", "niareagansc", "age_death", "apoe_genotype"]]
    df = pd.DataFrame({ind: us[X.row],
                       "module": adata.var.index.values[X.col],
                       "count": X.data})
    return pd.merge(df, md)


L = {x.split("_")[1].split(".")[0]: anndata.read(x) for x in os.listdir(".") if x.startswith("agg")}
df = pd.concat({k: agg2df(v) for k, v in L.items()}).reset_index().rename({"level_0": "CellType"}, axis=1)
del df["level_1"]

mf = pd.concat({x.replace(".bed.gz", ""): pd.read_csv(x, sep="\t", header=None) for x in os.listdir(".") if x.endswith(".bed.gz")})
mf.columns = ["chrom", "begin", "end", "module"]
mf = mf.reset_index().rename({"level_0": "CellType"}, axis=1)
cf = mf.groupby(["CellType", "module"])["chrom"].count().reset_index().rename({"chrom": "n_peaks"}, axis=1)

tf = pd.merge(cf, df)
