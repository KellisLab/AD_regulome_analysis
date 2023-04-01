#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from functools import reduce
from tqdm.auto import tqdm
import scipy.stats

def abc_distance(distance, hic_power=0.87):
    scale = -4.80 + 11.63 * hic_power
    offset = np.clip(np.abs(distance), 5000, np.Inf)
    return np.exp(scale + -1 * hic_power * np.log(offset + 1))

def read_interactions(ifname, distance_bin=1000.):
    dfi = pd.read_csv(ifname, sep="\t")
    dfi["dbin"] = np.round(dfi["distance"] / distance_bin).astype(int)
    return dfi

def read_activity(dfi, fname, celltype, columns_to_copy):
    dfp = pd.read_csv(fname, sep="\t")
    dfp = dfp.loc[dfp["CellType"] == celltype, :]
    dfp.index = dfp["peak"].values
    for col in columns_to_copy:
        dfi[col] = dfp.loc[dfi["peak"].values, col].values
    return dfi

def add_cor(dfi, fname, activity="log1p_total_counts", abc=False):
    df = pd.read_csv(fname, sep="\t")
    df = df.loc[df["original_index"].isin(dfi.index.values), :]
    dfi[fname] = np.nan
    dfi.loc[df["original_index"].values, fname] = df["cor"].values
    dfi = dfi.loc[~dfi[fname].isna(), :].copy()
    if activity in dfi.columns:
        dfi[fname] *= dfi[activity].values
    if abc:
        dfi[fname] *= abc_distance(dfi["distance"].values)
    return dfi

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


def train(df, columns, C=1e-4):
    X = df.loc[:, columns].values
    y = df["label"].values
    model = LogisticRegression(class_weight="balanced", n_jobs=-1, C=C)
    model.fit(X, y)
    py = model.predict_proba(X)[:, list(model.classes_).index(1)]
    df["py"] = py * df["log1p_total_counts"]
    return abc(df.loc[df["label"] == 1, :], col="py")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", nargs="+", required=True)
    ap.add_argument("--activity", required=True)
    ap.add_argument("--interactions", required=True)
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--activity-columns", default=["mean_counts", "log1p_total_counts"], nargs="+")
    ap.add_argument("-o", "--output", required=True)
    args = vars(ap.parse_args())
    ### read interactions
    dfi = read_interactions(args["interactions"])
    read_activity(dfi, fname=args["activity"],
                  celltype=args["celltype"],
                  columns_to_copy=args["activity_columns"])
    for fname in args["input"]:
        dfi = add_cor(dfi, fname)
    out = train(dfi, columns=args["input"])
    out.to_csv(args["output"], sep="\t")
