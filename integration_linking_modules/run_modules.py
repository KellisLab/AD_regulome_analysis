#!/usr/bin/env python3
import argparse
import anndata
import epiclust as ec
from epiclust.gene_distance import df_to_pyranges, peak_names_to_var
import os
def run_module(adata, power_list, batch="feature_types", n_neighbors=10,
               filter_z=2, filter_pct=0.0, gtf="",
               resolution=1.0, max_comm_size=10000):
    import json
    glist = []
    if os.path.exists(gtf):
        print("Annotating peaks")
        pr = df_to_pyranges(peak_names_to_var(adata.var.index.values), chrom="seqname", left="start", right="end")
        ft = adata.var[batch].values.astype(str)
        del adata.var[batch]
        adata.var[batch] = ft
        adata.var.loc[pr.df["name"].values, batch] = ec.gtf.annotate_ranges(pr, gtf)
    for power in power_list:
        print("Power:", power, "Fitting")
        ec.fit(adata, power=power, batch=batch)
        print("Power:", power, "Finding neighbors")
        ec.neighbors(adata, key_added="pow_%.2f" % power, n_neighbors=n_neighbors)
        glist.append("pow_%.2f" % power)
    print("Filtering nodes")
    V = ec.filter_var(adata, glist, z=filter_z, pct=filter_pct)
    adata = adata[:, V].copy()
    largs = {"resolution": resolution}
    if max_comm_size > 0:
        largs["max_comm_size"] = max_comm_size
    print("Clustering with args %s" % json.dumps(largs))
    ec.leiden(adata, ["pow_%.2f" % power for power in power_list], **largs)
    return adata

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--gtf", default="")
    ap.add_argument("--batch", default="feature_types")
    ap.add_argument("-p", "--power", type=float, nargs="+",
                    default=[0, 0.25, 0.5, 0.75, 1])
    ap.add_argument("-n", "--n-neighbors", type=int, default=10)
    ap.add_argument("-z", "--filter-z", type=float, default=2.)
    ap.add_argument("-f", "--filter-pct", type=float, default=0.)
    ap.add_argument("-r", "--resolution", type=float, default=1.)
    ap.add_argument("--max-comm-size", type=int, default=0)
    args = vars(ap.parse_args())
    adata = run_module(anndata.read(args["input"]), power_list=args["power"],
                       gtf=args["gtf"],
                       n_neighbors=args["n_neighbors"],
                       filter_z=args["filter_z"],
                       filter_pct=args["filter_pct"],
                       resolution=args["resolution"],
                       max_comm_size=args["max_comm_size"])
    print("Writing to H5AD")
    adata.write_h5ad(args["output"], compression="gzip")
