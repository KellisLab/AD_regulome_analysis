.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38

proj <- loadArchRProject(path = "Save-Proj3")

## export gene matrix
mtx = getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
GeneMtx<-(assay(mtx,"GeneScoreMatrix"))
genelist = getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  ignoreCase = TRUE
) 
rownames(GeneMtx) = genelist
saveRDS(GeneMtx,"scATAC.gene.mtx.rds")

## peak calling
proj <- addGroupCoverages(ArchRProj = proj,groupBy = "Celltype")
pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    peakMethod = "Macs2",
    groupBy = "Celltype",
    force = TRUE,
    maxPeaks=300000,
    cutOff = 0.01,
    pathToMacs2 = pathToMacs2
)
getPeakSet(proj)

proj <- addPeakMatrix(proj)  # add peak matrix

peakmtx = getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(peakMtx,"scATAC.peak.mtx.rds")

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj4", load = FALSE)

# export peaks of each cell type
celltypes<-unique(proj$Celltype)
for (i in 1:length(celltypes)){
        OutName=paste0("./Exported_peaks/Export.peak.",celltypes[[i]],".xls")
        InName=paste0("/net/bmc-lab5/data/kellis/users/xiongxs/Proj_snATAC_AD/20210205_PFC.reprocess/2.UseArchR.Re2/Save-Proj5/PeakCalls/",celltypes[[i]],"-reproduciblePeaks.gr.rds");
        peak<-readRDS(InName);
        peak.gr<-as.data.frame(peak);
        write.table(peak.gr,file = OutName,append=FALSE,quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".",row.names = FALSE, col.names = TRUE)
}

# cell type specific peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")  # define cell type specific peaks
saveRDS(markerList,'Celltype.specific.peaks.rds')

## use rGREAT to annotate the cell-type specific peaks, use all the peaks as background
library(ggplot2)
library(rGREAT)
library("bedr")

df = read.table("B.Exported_peaks/Export.peak.All.xls",head=T)
df = df[,c(1,2,3)]
allpeaks = df
allpeaks.loc = paste0(allpeaks[,1],":",allpeaks[,2],"-",allpeaks[,3])
allpeaks.loc = bedr.sort.region(allpeaks.loc)

for(ct in celltypes){
 df = read.table(paste0("B.Exported_peaks/Celltype.specific.peaks.",ct,".xls"),head=T)

# do an overlap
 df.loc = paste0(df[,1],":",df[,2],"-",df[,3])
 df.loc = bedr.sort.region(df.loc)
 df.intersect = bedr.join.region(df.loc,allpeaks.loc)

 df = unique(df.intersect[,2:4])
 colnames(df) = colnames(allpeaks)
 df = df[df$seqnames != ".",]
 df$start = as.numeric(df$start)
 df$end = as.numeric(df$end)

 job = submitGreatJob(df,bg=allpeaks,species="hg38")
 tb = getEnrichmentTables(job)
 saveRDS(tb,paste0("B.Exported_peaks/Celltype.specific.peaks.",ct,".specific.rGREAT.out.rds"))
}





