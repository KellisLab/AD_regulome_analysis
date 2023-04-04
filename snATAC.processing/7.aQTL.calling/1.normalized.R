.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(SummarizedExperiment)
library(edgeR)
library(ggplot2)
library(qvalue)

ct.list = c("Ast","Ex","In","Microglia","Oligo","OPC","PerEndo")
anno = read.table("/net/bmc-lab5/data/kellis/users/xiongxs/Proj_snATAC_AD/Samp_info/PFC.sampInfo.Info.txt",head=T)
anno$sampid = gsub("-",".",anno$SampID)
anno$ProjId2 = paste0(anno$ProjId,"_",anno$ProjId)

anno = anno[!duplicated(anno$ProjId),]
rownames(anno) = anno$sampid

# GT-available samples
gt.samp = read.table("/net/bmc-lab5/data/kellis/users/xiongxs/Proj_snATAC_AD/Genotypes/GT.projecid.added.txt",head=F)
gt.samp$samp = paste0(gt.samp$V1,"_",gt.samp$V1)

for(ct in ct.list){
 mtx.obj = readRDS(paste0(ct,".peakMat.bysample.pseudobulk.rds"))
 mtx = assays(mtx.obj)$PeakMatrix
 peaks = as.data.frame(rowData(mtx.obj))
 peaks$peakid = paste(peaks$seqnames,peaks$start,peaks$end,sep="_")
 peaks = peaks[,c(-2)]
 write.table(peaks,paste0(ct,".peaks.hg38.bed"),row.names=F,sep="\t",quote=F)
 rownames(mtx) = peaks$peakid
 filter = apply(mtx,2,function(x) sum(x)>=10000)
 mtx = mtx[,filter]
 filter <- apply(mtx, 1, function(x) length(x[x>=2])>=ncol(mtx)*0.5)
 mtx <- mtx[filter,]
 peaks <- rownames(mtx)
 y <- DGEList(counts=mtx)
 ## perform M-value normalization
 y <- calcNormFactors(y)
 y = voom(y)
 y.mtx = y$E
 
 anno.samp = anno[anno$SampID %in% colnames(y.mtx),]
 y.mtx = y.mtx[,as.character(anno.samp$SampID)] 
 all(colnames(y.mtx) == anno.samp$SampID)   ## check names are the same
 colnames(y.mtx) = anno.samp$ProjId2
 y.mtx = y.mtx[,colnames(y.mtx) %in% gt.samp$samp]

 write.table(data.frame(peakid=rownames(y.mtx),y.mtx),paste0(ct,".TMM.voom.normed.mtx"),row.names=F,sep="\t",quote=F)
}

