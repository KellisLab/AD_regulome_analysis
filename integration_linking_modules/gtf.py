#!/usr/bin/env python3

def load_gtf(gtf, max_score=5, min_score=1, peaks=None):
    import pandas as pd
    import gtfparse
    import anndata
    import numpy as np
    if type(gtf) != pd.DataFrame:
        gtf = gtfparse.read_gtf(gtf)
    gtf = gtf.loc[gtf["feature"] == "gene", :].drop_duplicates(["gene_id"])
    gtf.index = gtf["gene_name"].values
    gtf.index = anndata.utils.make_index_unique(gtf.index)
    gtf["tss"] = gtf["start"].values
    gtf["tts"] = gtf["end"].values
    gtf.loc[gtf["strand"] == "-", "tss"] = gtf["end"]
    gtf.loc[gtf["strand"] == "-", "tts"] = gtf["start"]
    gtf["feature_types"] = "Gene Expression"
    diff = 1/np.abs(gtf["end"].values - gtf["start"].values).clip(1, np.inf)
    raw_score = (diff - np.min(diff))/(np.max(diff) - np.min(diff))
    gtf["gene_length_score"] = (max_score - min_score) * raw_score + min_score
    return gtf

