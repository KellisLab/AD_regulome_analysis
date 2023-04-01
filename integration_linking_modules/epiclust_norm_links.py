#!/usr/bin/env python3

import scipy.stats
import numpy as np
import pandas as pd

def activity_zscore(df, col="mean_counts"):
    import numpy as np
    import scipy.stats
    log_counts = np.log(df.drop_duplicates("peak")[col].values)
    mlogp = -scipy.stats.norm.logcdf(np.mean(log_counts) - np.log(df[col].values),
                                     scale=np.std(log_counts))
    return mlogp

def abc_distance(distance, hic_power=0.87):
    import numpy as np
    scale = -4.80 + 11.63 * hic_power
    offset = np.clip(np.abs(distance), 5000, np.Inf)
    return np.exp(scale + -1 * hic_power * np.log(offset + 1))

def predict(df, prob_factor=1, geometric=True):
    import numpy as np
    df = df.loc[df["label"] == 1, :].copy()
    ### exp(prob) / (1-exp(prob)) =
    df["adj_prob"] = np.exp(df["prob"].values).clip(-1+1e-16, 1-1e-16)
    df["adj_prob"] = np.arctanh(df["adj_prob"].values) ### 0-1 --> 0-18.7
    df["adj_act"] = activity_zscore(df)
    df["abc_dist"] = abc_distance(df["distance"].values)
    if geometric:
        df["base"] = (df["adj_prob"] * df["adj_act"])**0.5
    else:
        df["base"] = 2/(1/df["adj_prob"] + 1/df["adj_act"])
    df["base"] *= df["abc_dist"]
    gf = df.groupby("gene").agg(total_base=("base", "sum"))
    df["total_base"] = gf.loc[df["gene"].values, "total_base"].values
    df["score"] = df["base"] / df["total_base"]
    df.index = df["interaction_index"].values
    return df.loc[:, ["peak", "gene", "score", "prob", "adj_prob", "adj_act", "abc_dist", "distance"]].sort_values("score", ascending=False)

def generate_predictors(df, vf, act_col="log1p_total_counts"):
    df = df.loc[df["label"] == 1,:].copy()
    df.index = df["interaction_index"].values
    df["one"] = 1
    df["eprob"] = np.exp(df["prob"].values)
    df["aprob"] = df["eprob"] / (1 - df["eprob"].clip(0, 1-1e-16))
    df["geom"] = (df[act_col] * df["eprob"])
    df["harm"] = 2/(1/df[act_col] + 1/df["eprob"])
    df["abc_dist"] = abc_distance(df["distance"])
    out = {}
    out["distance_raw"] = auprc(df, vf, col="abc_dist")
    out["distance_abc"] = auprc(abc(df, col="one"), vf)
    out["abc"] = auprc(abc(df, col=act_col), vf)
    out["eprob_abc"] = auprc(abc(df, col="eprob"), vf)
    out["eprob_raw"] = auprc(df, vf, col="eprob")
    out["aprob_abc"] = auprc(abc(df, col="aprob"), vf)
    out["aprob_raw"] = auprc(df, vf, col="aprob")
    out["activity_raw"] = auprc(df, vf, col=act_col)
    out["geom"] = auprc(abc(df, col="geom"), vf)
    out["harm"] = auprc(abc(df, col="harm"), vf)
    return out

def abc(df, col="log1p_total_counts"):
    import numpy as np
    if "label" in df.columns:
        df = df.loc[df["label"] == 1, :].copy()
    df["base"] = df[col] * abc_distance(df["distance"].values)
    gf = df.groupby("gene").agg(total_base=("base", "sum"))
    df["total_base"] = gf.loc[df["gene"].values, "total_base"].values
    df["score"] = df["base"] / df["total_base"]
    if "interaction_index" in df.columns:
        df.index = df["interaction_index"].values
    return df.loc[:, ["peak", "gene", "score", "distance"]].sort_values("score", ascending=False)

def auprc(df, tf, col="score"):
    I = np.intersect1d(df.index.values, tf.index.values)
    scores = df.loc[I, col]
    labels = tf.loc[I, "ClusterSummit"].values
    return average_precision_score(labels, scores)

def ttest(df, col="cor", nbin=5):
    df["dbin"] = np.round(df["distance"] / 1000).astype(int)
    udbin, dbin_inv = np.unique(df["dbin"].values, return_inverse=True)
    means, stds, counts = {}, {}, {}
    for i in tqdm(np.arange(len(udbin))):
        dbin = udbin[i]
        flag = np.abs(i - dbin_inv) <= nbin
        means[dbin] = df[col].values[flag].mean()
        stds[dbin] = df[col].values[flag].std()
        counts[dbin] = np.sum(flag)
    return pd.DataFrame({"counts": counts, "means": means, "stds": stds})


def load_plac(plac_dir="/net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/evaluation/plac"):
    tbl = {}
    for fname in os.listdir(plac_dir):
        if fname.endswith(".tsv.gz"):
            tbl[fname.replace(".tsv.gz", "")] = pd.read_csv(os.path.join(plac_dir, fname), sep="\t", index_col=0)
    return tbl

def auprc_results(linking_dir=".", score_cutoff=0.001, plac_dir="/net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/evaluation/plac"):
    auprc = {}
    plac = load_plac(plac_dir)
    for fname in os.listdir(linking_dir):
        if fname.startswith("predicted_") and fname.endswith(".tsv.gz"):
            ct = fname.replace("predicted_", "").replace(".tsv.gz", "")
            df = pd.read_csv(os.path.join(linking_dir, fname), sep="\t", index_col=0)
            if ct in plac:
                print(ct)
                auprc[ct] = {}
                I = np.intersect1d(df.index.values, plac[ct].index.values)
                labels = plac[ct].loc[I, "ClusterSummit"].values
                auprc[ct]["distance"] = average_precision_score(labels, np.abs(df.loc[I, "distance"].values).clip(1, np.inf)**(-1))
                I = np.intersect1d(df.loc[df["score"] >= score_cutoff, :].index.values, plac[ct].index.values)
                labels = plac[ct].loc[I, "ClusterSummit"].values
                auprc[ct]["linking"] = average_precision_score(labels, df.loc[I, "score"].values)
    return auprc
