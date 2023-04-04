.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ggplot2)
library(BuenColors)
library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38

proj <- readRDS("3.test/ArchR.TSS1.iter3.sampCellAll.sampFeat80000.Seurat.res1.dim50.sampling.step3.rds")

mytheme <- theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio=1)

mytheme5<-theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio = 1,axis.text = element_text(face="bold", size=rel(1.5)),axis.title = element_text(face="bold", size=rel(1.5)),legend.text=element_text(face="bold",size=15),plot.title=element_text(face="bold",size=rel(1.5)),legend.title =element_text(face="bold", size=rel(1.5)))

# assign
## assign cell type
proj$Celltype[proj$Clusters %in% c("C12","C13","C14","C11")]<-"Oligo"
proj$Celltype[proj$Clusters %in% c("C2","C4","C6")] <- "Microglia"
proj$Celltype[proj$Clusters %in% c("C19","C20")] <- "Ex"
proj$Celltype[proj$Clusters %in% c("C3")] <- "PerEndo"
proj$Celltype[proj$Clusters %in% c("C15","C16","C17","C18")] <- "In"
proj$Celltype[proj$Clusters %in% c("C10")] <- "OPC"
proj$Celltype[proj$Clusters %in% c("C7","C8","C9")] <- "Ast"
proj$Celltype[is.na(proj$Celltype)] = "Deidentified"

meta = proj@cellColData
meta$cellID = rownames(meta)
## attach info
info = read.table("0.data/SraRunTable.txt",head=T,sep=",")
info = info[,c(1:2,15,25,27,28)]
metaInfo = merge(meta,info,by.x="Sample",by.y="Run",all.x=T)

## check before adding
all(metaInfo$cellID == rownames(meta)) # false
rownames(metaInfo) = metaInfo$cellID
metaInfo = metaInfo[rownames(meta),]
all(metaInfo$cellID == rownames(meta)) # TRUE

proj<-addCellColData(
  ArchRProj = proj,
  data = as.character(metaInfo$AGE),
  name = "AGE",
  cells = proj$cellNames,
  force = FALSE
)
proj<-addCellColData(
  ArchRProj = proj,
  data = as.character(metaInfo$diagnosis),
  name = "diagnosis",
  cells = proj$cellNames,
  force = FALSE
)
proj<-addCellColData(
  ArchRProj = proj,
  data = as.character(metaInfo$pmi),
  name = "pmi",
  cells = proj$cellNames,
  force = FALSE
)
proj<-addCellColData(
  ArchRProj = proj,
  data = as.character(metaInfo$sex),
  name = "sex",
  cells = proj$cellNames,
  force = FALSE
)

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj2-TSS1", load = FALSE)







# plot against TSS
samp.list = as.character(unique(tmp$Sample.Name))
coord = data.frame()
for (i in 1:length(samp.list)){
        ind=which(tmp$Sample.Name == samp.list[i])
	ncell = nrow(tmp[ind,])
        mid.tss = median(tmp[ind,]$TSSEnrichment)
        mid.nFrags = median(tmp[ind,]$nFrags)
        coord = rbind(coord,data.frame(Sample=samp.list[i],TSSEnrichment.mid=mid.tss,nFrags.mid=mid.nFrags,nCell=ncell))
}

coord.info = merge(coord,info,by.x="Sample",by.y="Sample.Name",all.x=T)
coord.info = coord.info[!duplicated(coord.info$Sample),]

write.table(coord.info,"Sample.TSS.by.nFrag.mid.txt",sep="\t",row.names=F,quote=F)

sum = coord.info

pdf("NG2021.reprocess.0919.nUMI.vs.TSS.pdf")
ggplot(sum,aes(x=log10(nFrags.mid),y=TSSEnrichment.mid,color=diagnosis,size=nCell)) + geom_point() + mytheme5+
  scale_size_continuous(range = c(4,8))
ggplot(sum,aes(x=log10(nFrags.mid),y=TSSEnrichment.mid,color=sex,size=nCell)) + geom_point() + mytheme5+
  scale_size_continuous(range = c(4,8))
dev.off()

## cell proportion vs. sample
Proportion = as.data.frame(table(tmp$Sample.Name,tmp$Clusters))

prop = Proportion
colnames(prop) = c("Samp","Cls","Freq")
library(plyr)
prop<-ddply(prop,"Cls",transform,perc=Freq/sum(Freq)*100)
prop<-ddply(prop,"Samp",transform,perc2=Freq/sum(Freq)*100)
pdf("NG2021.reprocess.0919.proportion.pdf")
ggplot(prop, aes(x=Cls, y=perc, fill=Samp)) +
  geom_bar(stat="identity", colour="black") +
  theme(aspect.ratio = 0.5)

Proportion = as.data.frame(table(tmp$diagnosis,tmp$Clusters))

prop = Proportion
colnames(prop) = c("Samp","Cls","Freq")
library(plyr)
prop<-ddply(prop,"Cls",transform,perc=Freq/sum(Freq)*100)
prop<-ddply(prop,"Samp",transform,perc2=Freq/sum(Freq)*100)
ggplot(prop, aes(x=Cls, y=perc, fill=Samp)) +
  geom_bar(stat="identity", colour="black") +
  theme(aspect.ratio = 0.5)

dev.off()




