.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 16)
addArchRGenome("hg38")   ## use hg38

## read in data
inputFiles <- list.files("./0.data", pattern = ".tsv.gz$", full.names = TRUE)
library(stringr)
names(inputFiles) <- str_extract(inputFiles,"SRR[0-9]+")      ## extract name of the sample
inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 1 ,  ## set to 1
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Save-Proj1",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

save.image(file = "ArrowProject.step2.RData")



