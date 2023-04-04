.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38

TSS=1

proj <- loadArchRProject(path = "../Save-Proj2-TSS1")
meta = proj@cellColData

anno = read.table("/net/bmc-lab5/data/kellis/users/xiongxs/Reference/human/hg38/RoadMap/ChromHMM.anno.txt",head=T,sep="\t")
PFC = read.table("/net/bmc-lab5/data/kellis/users/xiongxs/Reference/human/hg38/RoadMap/E073_18_core_K27ac_hg38lift_segments.bed.gz",head=F,sep="\t")
PFC = PFC[PFC$V1 != "chrM" & PFC$V1 != "chrY",]
PFC$state = gsub("E","",PFC$V4)
PFC = merge(PFC,anno,by.x="state",by.y="STATE.NO.",all.x=T)

state.list = levels(PFC$MNEMONIC)

for(state in state.list){
print(state)

peak = PFC[PFC$MNEMONIC == state,]

peak.gr = GRanges(
  peak$V1,
  IRanges(peak$V2,peak$V3)
)

proj = addPeakSet(
  ArchRProj = proj,
  peakSet = peak.gr,
  force = T
)

proj = addPeakMatrix(
  ArchRProj = proj,
  force = TRUE
)

meta0 = proj@cellColData
meta[,state] = meta0$FRIP
}

saveRDS(meta,"ChromHMM_Morabito.TSS1.meta.clean.rds")
write.table(meta,paste0("ChromHMM_Morabito.TSS",TSS,".meta.clean.txt"),sep="\t",quote=F,row.names=F)







