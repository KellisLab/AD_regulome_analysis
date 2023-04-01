#!/usr/bin/env Rscript
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

counts = Read10X_h5("./filtered_feature_bc_matrix.h5")
fragpath = "./atac_fragments.tsv.gz"


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# create a Seurat object containing the RNA adata
s.obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
s.obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks[,colnames(s.obj)],
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(s.obj) = "RNA"
s.obj = SCTransform(s.obj)
s.obj = RunPCA(s.obj)
DefaultAssay(s.obj) = "ATAC"
s.obj = FindTopFeatures(s.obj, min.cutoff = 5)
s.obj = RunTFIDF(s.obj)
s.obj = RunSVD(s.obj)

s.obj = FindMultiModalNeighbors(
  object = s.obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)


s.obj <- RunUMAP(
  object = s.obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)
