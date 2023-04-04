.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")

library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38

proj <- loadArchRProject(path = "B.Save-Proj3")

proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = F
)
cA.df = as.data.frame(cA)

gr = metadata(cA)[[1]]
gr.df = data.frame(gr)

write.table(cA.df,"CoAccessibility.0904.txt",quote=F,row.names=F,sep="\t")
write.table(data.frame(idx=rownames(gr.df),gr.df),"CoAccessibility.peak.idx.0904.txt",quote=F,row.names=F,sep="\t")




