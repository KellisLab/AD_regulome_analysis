#!/usr/bin/env python3

import argparse
import os
import sys
import re
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.interpolate import interp1d

from matplotlib import rcParams

# figure size in inches
rcParams['figure.figsize'] = 11.7,8.27

def run(link_df, plac, cutoff=0, linking_method="Linking"):
    pf = pd.merge(link_df, plac)
    tbl = {"Closest Gene": pf.loc[pf["is_closest_gene"], :],
           "Distal": pf.loc[~pf["is_closest_gene"], :],
           "All": pf}
    pr_curve = []
    auprc = []
    roc = []
    auroc = []
    for vname, val_f in tbl.items():
        print(val_f)
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
            auprc.append(pd.DataFrame({"AUPRC": average_pr,
                                       "method": mname,
                                       "validation": vname}, index=["%s %s" % (mname, vname)]))
            pr_curve.append(pd.DataFrame({
                "precision": precision,
                "recall": recall,
                "method": mname,
                "validation": vname
            }))
            fpr, tpr, _ = roc_curve(mvals["labels"], mvals["scores"])
            roc.append(pd.DataFrame({
                "FPR": fpr,
                "TPR": tpr,
                "method": mname,
                "validation": vname
            }))
            auroc.append(pd.DataFrame({
                "AUROC": roc_auc_score(mvals["labels"], mvals["scores"]),
                "method": mname,
                "validation": vname
            }, index=["%s %s" % (mname, vname)]))
    return pd.concat(pr_curve).reset_index(drop=True), pd.concat(auprc).reset_index(drop=True), pd.concat(roc).reset_index(drop=True), pd.concat(auroc).reset_index(drop=True)

# def forceAspect(ax,aspect=1):
#     extent = ax.get_images().get_extent()
#     ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--plac", required=True)
    ap.add_argument("--link", nargs="+")
    ap.add_argument("--cutoff", default=0, type=float)
    ap.add_argument("--interpolation", default="nearest")
    ap.add_argument("--palette", default="deep")
    ap.add_argument("--validation-prcurve", default="All")
    ap.add_argument("--figdir", default=".")
    args = vars(ap.parse_args())
    plac = pd.read_csv(args["plac"])
    prcurve, auprc, auroc, roc = [], [], [], []
    for link_tsv in args["link"]:
        celltype = os.path.basename(link_tsv)
        celltype = re.sub(r".bed.gz$", "", celltype).split("_")
        method = celltype[0]
        celltype = celltype[1]
        if celltype in pd.unique(plac["CellType"]):
            print(celltype)
            link_df = pd.read_csv(link_tsv, sep="\t", header=None)
            link_df = link_df.iloc[:, [0, 1, 2, 3, 4]]
            link_df.columns = ["chrom", "start", "end", "gene", "score"]
            link_df["peak"] = ["%s:%d-%d" % (chrom, start, end) for chrom, start, end in zip(link_df["chrom"], link_df["start"], link_df["end"])]
            local_prcurve, local_auprc, local_roc, local_auroc = run(link_df, plac=plac.loc[plac["CellType"] == celltype, :], linking_method=method)
            local_prcurve["CellType"] = celltype
            local_auprc["CellType"] = celltype
            local_roc["CellType"] = celltype
            local_auroc["CellType"] = celltype
            local_roc["ParentMethod"] = method
            local_prcurve["ParentMethod"] = method
            prcurve.append(local_prcurve)
            auprc.append(local_auprc)
            auroc.append(local_auroc)
            roc.append(local_roc)
    prcurve = pd.concat(prcurve).reset_index(drop=True)
    auprc = pd.concat(auprc).reset_index(drop=True)
    roc = pd.concat(roc).reset_index(drop=True)
    auroc = pd.concat(auroc).reset_index(drop=True)
    print(auprc)
    ### Have separate subplots per cell type
    ### Then do Distal/Nearest/All for x, Method for hue
    method_colors = sns.color_palette(args["palette"])
    method_colors = {k: method_colors[i] for i, k in enumerate(pd.unique(auprc["method"]))}
    celltypes = pd.unique(auprc["CellType"])
    fig_nrow = 1
    fig_ncol = 4
    ### Plot AUPRC barplots
    fig, axes = plt.subplots(fig_nrow, fig_ncol, sharey=True, sharex=False, figsize=(15, 6)) ### want explicit labels
    for i, ct in enumerate(celltypes):
        B = sns.barplot(data=auprc.loc[auprc["CellType"] == ct,:], x="validation", y="AUPRC", hue="method", palette=method_colors, ax=axes[i])
        B.set(title=ct, xlabel=None)

    fig.tight_layout()
    print("Saving figure to", args["figdir"])
    fig.savefig(os.path.join(args["figdir"], "auprc.pdf"), dpi=600, bbox_inches="tight")
    fig_nrow = 1
    fig_ncol = 4
    ### Plot AUPRC barplots
    fig, axes = plt.subplots(fig_nrow, fig_ncol, sharey=True, sharex=False, figsize=(15, 6)) ### want explicit labels
    for i, ct in enumerate(celltypes):
        B = sns.barplot(data=auroc.loc[auprc["CellType"] == ct,:], x="validation", y="AUROC", hue="method", palette=method_colors, ax=axes[i])
        axes[i].set_ylim(0.5, 0.65)
        B.set(title=ct, xlabel=None)
        sns.move_legend(B, "upper right")
    fig.tight_layout()
    print("Saving figure to", args["figdir"])
    fig.savefig(os.path.join(args["figdir"], "auroc.pdf"), dpi=600, bbox_inches="tight")
    ### Plot PR curves
    fig, axes = plt.subplots(fig_nrow, fig_ncol, figsize=(8, 32), sharey=True)
    legend = {}
    for i, ct in enumerate(celltypes):
        df = prcurve.loc[prcurve["CellType"] == ct, :]
        df = df.loc[df["validation"] == args["validation_prcurve"], :]
        #GOOD_METHODS = df.sort_values("method", ascending=False).drop_duplicates("ParentMethod")
        #df.groupby(["ParentMethod", "method"]).
        interp = {k: interp1d(df.loc[df["method"] == k, "recall"].values,
                              df.loc[df["method"] == k, "precision"].values,
                              kind=args["interpolation"])
                  for k in pd.unique(df["method"])}
        rec = np.sort(pd.unique(df["recall"]))
        print(interp.keys())
        for method in interp.keys():
            prec = interp[method](rec)
            axes[i].plot(rec, prec, color=method_colors[method], label=method)
            #legend[method] = axes[i].fill_between(rec, prec_d, prec_l, where=prec_l < prec_d)
            #legend["Linking"] = axes[i].fill_between(rec, prec_d, prec_l, where=prec_l > prec_d, color=method_colors["Linking"])
        axes[i].set(title=ct)
        axes[i].set_aspect(1)
        #axes[i].legend()
    #fig.legend([legend[k] for k in legend.keys()], [k for k in legend.keys()])
    fig.tight_layout()
    fig.savefig(os.path.join(args["figdir"], "pr_curve.pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(args["figdir"], "pr_curve.png"), dpi=800, bbox_inches="tight")
    ### Plot ROC curves
    fig, axes = plt.subplots(fig_nrow, fig_ncol, figsize=(8, 32), sharey=True)
    legend = {}
    for i, ct in enumerate(celltypes):
        df = roc.loc[roc["CellType"] == ct, :]
        df = df.loc[df["validation"] == args["validation_prcurve"], :]
        interp = {k: interp1d(df.loc[df["method"] == k, "FPR"].values,
                              df.loc[df["method"] == k, "TPR"].values,
                              kind=args["interpolation"])
                  for k in pd.unique(df["method"])}
        FPR = np.sort(pd.unique(df["FPR"]))
        print(interp.keys())
        for method in interp.keys():
            TPR = interp[method](FPR)
            axes[i].plot(FPR, TPR, color=method_colors[method], label=method)
            #legend[method] = axes[i].fill_between(rec, prec_d, prec_l, where=prec_l < prec_d)
            #legend["Linking"] = axes[i].fill_between(rec, prec_d, prec_l, where=prec_l > prec_d, color=method_colors["Linking"])
        axes[i].set(title=ct)
        axes[i].set_aspect(1)
        #axes[i].legend()
    #fig.legend([legend[k] for k in legend.keys()], [k for k in legend.keys()])
    fig.tight_layout()
    fig.savefig(os.path.join(args["figdir"], "roc_curve.pdf"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(args["figdir"], "roc_curve.png"), dpi=800, bbox_inches="tight")

        #B = sns.barplot(data=prcurve.loc[prcurve["CellType"] == ct,:], x="validation", y="AUPRC", hue="method", ax=axes[row_idx, col_idx])
