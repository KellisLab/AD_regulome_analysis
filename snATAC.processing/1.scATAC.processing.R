.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
addArchRThreads(threads = 16)
addArchRGenome("hg38")   ## use hg38

## read in data
inputFiles <- list.files("./data", pattern = ".tsv.gz$", full.names = TRUE)
library(stringr)
names(inputFiles) <- str_extract(inputFiles,"D19-[0-9_]+")      ## extract names of the samples
inputFiles

## create Arrow files that store large data
## Do QC with TSS enrichment >= 6 & nFrag > 1000
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 6, 
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

print(ArrowFiles)  # make sure you get the right ArrowFiles

## calculate doubletscores for each sample
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

## create an ArchR object
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Save-Proj1",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)  # make sure that you get the tile matrix

# remove doublets
proj <- filterDoublets(ArchRProj = proj)

# do LSI dimension reduction
ncell = length(proj$cellNames) 
iter = 3  # iteration of LSI; basically 3 or 4 should be good enough
# several parameters that can test
res = 0.2 # resolution
nVar = 50000  # try 10k, 20k, 50k, 80, etc
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", iterations = iter, clusterParams = list(resolution = c(res), sampleCells = (0.2 * ncell), n.start = 10), name = "IterativeLSI",varFeatures = nVar,force = TRUE)

# clustering and UMAP embedding 
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force = TRUE ,method = "Seurat")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force = TRUE)

# visualize 
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

pdf("test.pdf")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
dev.off()

# you can export the gene score matrix at the cluster level to help define cell types
geneMat <- getGroupSE(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  divideN = TRUE,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

mat.tmp<-(assay(geneMat2,"GeneScoreMatrix"))
genelist<-getFeatures(proj)
rownames(mat.tmp)<-genelist
write.table(data.frame("gene"=rownames(mat.tmp),mat.tmp),file = "GeneMat_by_clusters.Feat50k.txt",quote = FALSE, sep = "\t",row.names=FALSE)

# get the aggregated PsychENCODE marker matrix
system('perl get_markerMat.PEC_markers_2.byCluster.re.pl')

## save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj2", load = FALSE)




