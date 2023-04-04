# do TF footprinting
.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38
pathToMacs2 <- findMacs2()

proj <- loadArchRProject(path = "../B.Save-Proj6")

# get the positions of each TF
motifPositions <- getPositions(proj)

# generate pseudobulk
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Celltype4")

# read in TF regulator candidates from motif enrichment 
atac_chromvar_rna = read.table("../C.16.Positive.TF/ATAC_RNA_ChromVar.cor.0919.xls",head=T)
tf.list = as.character(atac_chromvar_rna[atac_chromvar_rna$Is.Cand.Reglator=='Yes',]$GeneScoreMatrix_matchName)
seZ.mat2 = read.table("../C.16.Positive.TF/TF.chromVar.matrix.0919.xls",head=T)
rownames(seZ.mat2) = seZ.mat2$TF
seZ.mat2 = seZ.mat2[tf.list,]
seZ.mat = seZ.mat2[,2:8]

reord<-order(apply(seZ.mat,1,which.max),decreasing = F)
seZ.mat.ord = seZ.mat2[reord,]

ct.list = unique(seZ.mat.ord$ct.max)

for(ct in ct.list){
print(ct)
tf.list.ct = as.character(seZ.mat.ord[seZ.mat.ord$ct.max == ct,]$TF)

tf.list.tmp = paste0(tf.list.ct,"_") # to avoid mis-match
markerMotifs <- unlist(lapply(tf.list.tmp, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "BATF3_120"]
print(data.frame(a=tf.list.tmp,b=markerMotifs))

# compute footprints
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Celltype4"
)

# plot with Tn5 bias corrected (Subtract)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  height = 5,
  plotName = paste0("Footprints-Subtract-Bias.",ct),
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "Divide",
  plotName = paste0("Footprints-Divide-Bias.",ct),
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj,
  normMethod = "None",
  plotName = paste0("Footprints-NoNorm-Bias.",ct),
  addDOC = FALSE,
  smoothWindow = 5
)
}









