.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16)
addArchRGenome("hg38")   ## use hg38

load("ArrowProject.step2.RData")

proj <- proj[proj@cellColData$TSSEnrichment>=1]

ncell.pre = 100000
ncell.final = 2*ncell.pre
proj <- addIterativeLSI(ArchRProj = proj,dimsToUse = 1:50, useMatrix = "TileMatrix", iterations = 3, clusterParams = list(resolution = c(1), sampleCells = ncell.pre, n.start = 10), name = "IterativeLSI",varFeatures = 80000,sampleCellsPre=ncell.pre,projectCellsPre=TRUE,sampleCellsFinal=ncell.final,force = TRUE)

proj <- addClusters(input = proj, dimsToUse = 1:50, reducedDims = "IterativeLSI",force = TRUE ,method = "Seurat",maxClusters=50)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force = TRUE)
#proj <- addTSNE(ArchRProj = proj, reducedDims = "IterativeLSI",perplexity = 30, name = "TSNE",force = TRUE)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
pdf("ArchR.TSS1.iter3.sampCellAll.sampFeat80000.Seurat.res1.dims50.sampling.pdf")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
dev.off()

saveRDS(proj,"ArchR.TSS1.iter3.sampCellAll.sampFeat80000.Seurat.res1.dim50.sampling.step3.rds")
