#!/usr/bin/env Rscript
library(Seurat)
library(Signac)

loadPeaks <- function(peak.dir) {
    PM = readRDS(paste0(peak.dir, "/Peak.matrix.re.rds"))
    rownames(PM) = sapply(strsplit(rowRanges(PM)$peakid, "_"), function(x) {
        paste0(x[1], ":", x[2], "-", x[3])
    })
    return(PM)
}
celltypes = list.files("~/data-group/Xushen_Sharing/ToBen/ATAC.matrix/matrix.by.major_celltype.20220816/", full.names=T, pattern="^[A-Z][A-Za-z]+$")
PML = lapply(setNames(celltypes, basename(celltypes)), loadPeaks)

cds = lapply(L, function(s.obj) { make_cicero_cds(as.cell_data_set(s.obj), reduced_coordinates=s.obj@reductions$umap@cell.embeddings) })
genome.df = as.data.frame(L[[1]]@assays$peaks@seqinfo)
genome.df = genome.df[-grep("^H[GS]", rownames(genome.df)),]
genome.df = data.frame(chr=paste0("chr", rownames(genome.df)), length=genome.df$seqlengths)
#rng = with(, data.frame(chr=paste0("chr", rownames(
