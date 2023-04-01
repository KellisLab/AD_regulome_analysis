#!/usr/bin/env python3

import argparse
import os
import re
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score, precision_recall_curve
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.interpolate import interp1d

from matplotlib import rcParams

# figure size in inches
rcParams['figure.figsize'] = 11.7,8.27

def run(link_df, plac, cutoff=0.001, linking_method="Linking"):
    pf = pd.merge(link_df, plac)
    tbl = {"Closest Gene": pf.loc[pf["is_closest_gene"], :],
           "Distal": pf.loc[~pf["is_closest_gene"], :],
           "All": pf}
    pr_curve = []
    auprc = []
    for vname, val_f in tbl.items():
        methods = {}
        methods["Distance"] = {
            "labels": val_f["ClusterSummit"].values,
            "scores": np.abs(val_f["distance"].values).clip(1, np.inf) ** -0.87
        }
        link_f = val_f.loc[val_f["score"] >= cutoff, :]
        methods[linking_method] = {
            "labels": link_f["ClusterSummit"].values,
            "scores": link_f["score"].values
        }
        for mname, mvals in methods.items():
            average_pr = average_precision_score(mvals["labels"], mvals["scores"])
            precision, recall, _ = precision_recall_curve(mvals["labels"], mvals["scores"])
            pr_curve.append(pd.DataFrame({
                "precision": precision,
                "recall": recall,
                "method": mname,
                "validation": vname
            }))
            auprc.append(pd.DataFrame({"AUPRC": average_pr,
                                       "method": mname,
                                       "validation": vname}, index=["%s %s" % (mname, vname)]))
    return pd.concat(pr_curve).reset_index(drop=True), pd.concat(auprc).reset_index(drop=True)



# def forceAspect(ax,aspect=1):
#     extent = ax.get_images().get_extent()
#     ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--plac", required=True)
    ap.add_argument("--link", nargs="+")
    ap.add_argument("--cutoff", default=0.001, type=float)
    ap.add_argument("--interpolation", default="nearest")
    ap.add_argument("--palette", default="deep")
    ap.add_argument("--validation-prcurve", default="All")
    ap.add_argument("--figdir", default=".")
    args = vars(ap.parse_args())
    plac = pd.read_csv(args["plac"])
    if "distance" in plac.columns:
        del plac["distance"]
    prcurve, auprc = {}, {}
    for link_tsv in args["link"]:
        celltype = os.path.basename(link_tsv)
        celltype = re.sub(r".tsv.gz$", "", celltype)
        celltype = re.sub(r"^predicted_", "", celltype)
        if celltype in pd.unique(plac["CellType"]):
            print(celltype)
            link_df = pd.read_csv(link_tsv, sep="\t", index_col=0)
            prcurve[celltype], auprc[celltype] = run(link_df, plac=plac.loc[plac["CellType"] == celltype, :], cutoff=args["cutoff"])
    prcurve = pd.concat(prcurve).reset_index(0).rename({"level_0": "CellType"}, axis=1)
    auprc = pd.concat(auprc).reset_index(0).rename({"level_0": "CellType"}, axis=1)
    ### Have separate subplots per cell type
    ### Then do Distal/Nearest/All for x, Method for hue
    method_colors = sns.color_palette(args["palette"])
    method_colors = {k: method_colors[i] for i, k in enumerate(pd.unique(auprc["method"]))}
    celltypes = pd.unique(auprc["CellType"])
    fig_nrow = 1
    fig_ncol = 4
    ### Plot AUPRC barplots
    fig, axes = plt.subplots(fig_nrow, fig_ncol, sharey=True, sharex=False) ### want explicit labels
    for i, ct in enumerate(celltypes):
        B = sns.barplot(data=auprc.loc[auprc["CellType"] == ct,:], x="validation", y="AUPRC", hue="method", palette=method_colors, ax=axes[i])
        B.set(title=ct, xlabel=None)

    fig.tight_layout()
    fig.savefig(os.path.join(args["figdir"], "auprc.pdf"), dpi=600, bbox_inches="tight")
    ### Plot PR curves
    fig, axes = plt.subplots(fig_nrow, fig_ncol, figsize=(8, 32), sharey=True)
    legend = {}
    for i, ct in enumerate(celltypes):
        df = prcurve.loc[prcurve["CellType"] == ct, :]
        df = df.loc[df["validation"] == args["validation_prcurve"], :]
        interp = {k: interp1d(df.loc[df["method"] == k, "recall"].values,
                              df.loc[df["method"] == k, "precision"].values,
                              kind=args["interpolation"])
                  for k in pd.unique(df["method"])}
        rec = np.sort(pd.unique(df["recall"]))
        prec_d = interp["Distance"](rec)
        prec_l = interp["Linking"](rec)
        axes[i].plot(rec, prec_d, color=method_colors["Distance"], label="Distance")
        axes[i].plot(rec, prec_l, color=method_colors["Linking"], label="Linking")
        legend["Distance"] = axes[i].fill_between(rec, prec_d, prec_l, where=prec_l < prec_d)
        legend["Linking"] = axes[i].fill_between(rec, prec_d, prec_l, where=prec_l > prec_d, color=method_colors["Linking"])
        axes[i].set(title=ct)
        axes[i].set_aspect(1)
        #axes[i].legend()
    #fig.legend([legend[k] for k in legend.keys()], [k for k in legend.keys()])
    fig.tight_layout()
    fig.savefig(os.path.join(args["figdir"], "pr_curve.pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(args["figdir"], "pr_curve.png"), dpi=800, bbox_inches="tight")

        #B = sns.barplot(data=prcurve.loc[prcurve["CellType"] == ct,:], x="validation", y="AUPRC", hue="method", ax=axes[row_idx, col_idx])
