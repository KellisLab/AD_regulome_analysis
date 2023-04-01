#!/usr/bin/env python3

import snapatac2 as snap

import os
import numpy as np
import pandas as pd
import argparse
def process_fragments(fragment_file, bin_size=500,
                      gene_annotation=None,
                      valid_barcodes=None, out_file=None):
    if out_file is None:
        out_file = os.path.basename(fragment_file).replace(".tsv.gz", ".h5ad")
    adata = snap.pp.import_data(
        fragment_file,
        gene_annotation,
        snap.genome.hg38,
        file=out_file,
        sorted_by_barcode=False)
    if valid_barcodes is not None:
        print("Subsetting to barcodes...")
        adata.subset(np.isin(adata.obs_names, [x.split("#")[1] for x in valid_barcodes]))
    else:
        snap.pp.filter_cells(adata, min_counts=5000, min_tsse=10, max_counts=50000)
    snap.pp.make_tile_matrix(adata, bin_size=bin_size)
    adata.close()
    return out_file

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-s", "--sample", required=True)
    ap.add_argument("-a", "--annotation", required=True)
    ap.add_argument("--sample-column", default="Sample")
    ap.add_argument("--barcode-column", default="obsname")
    ap.add_argument("-g", "--gene-annotation", default="/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
    ap.add_argument("-b", "--bin-size", type=int, default=500)
    args = vars(ap.parse_args())
    obs = pd.read_csv(args["annotation"], sep="\t")
    bc = obs.loc[obs[args["sample_column"]] == args["sample"], args["barcode_column"]].values
    process_fragments(fragment_file=args["input"],
                      bin_size=args["bin_size"],
                      gene_annotation=args["gene_annotation"],
                      valid_barcodes=bc,
                      out_file=args["output"])
