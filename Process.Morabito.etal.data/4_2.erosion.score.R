.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38
library(ggplot2)
library(BuenColors)

mytheme <-theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio = 1,axis.text = element_text(face="bold", size=rel(1.5)),axis.title = element_text(face="bold", size=rel(1.5)),legend.text=element_text(face="bold",size=15),plot.title=element_text(face="bold",size=rel(1.5)),legend.title =element_text(face="bold", size=rel(1.5)))
mytheme2 <-theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio = 1,axis.text = element_text(face="bold", size=rel(1)),axis.title = element_text(face="bold", size=rel(1)),legend.text=element_text(face="bold",size=10),plot.title=element_text(face="bold",size=rel(1)),legend.title =element_text(face="bold", size=rel(1)))

########### Full ChromHMM annotations
TSS = 1
print(TSS)

tmp = readRDS("ChromHMM_Morabito.TSS1.meta.clean.rds")

chrom.state = colnames(tmp)[22:39]

# clr; all cell types together
tmp2 = t(as.data.frame(tmp[,22:39]))
library(compositions)
tmp2.clr = clr(tmp2)
tmp2.clr = data.frame(tmp2.clr)
tmp2.clr.t = t(tmp2.clr)
# check if same order
all(gsub("X","",rownames(tmp2.clr.t)) == gsub("-|#",".",rownames(tmp)))
colnames(tmp2.clr.t) = paste0(colnames(tmp2.clr.t),".clr")
tmp = cbind(tmp,tmp2.clr.t)

## adding up the "deviation to average" of the 18 status, as the "erosion score"
## because the data is centered (CLR), making the mean to 0, therefore does not need to subtract to the mean again 

# try "undirectional" - adding up all the abs(deviation)
tmp$abs_deviation_sum = rowSums(abs(tmp2.clr.t))

# directed; assigning repressive regions to +, and the opened regions to -; therefore, a higher score means more opened in the represive regions
tmp2.clr.t.direct = tmp2.clr.t
opened = grep("Enh|Tss|Tx",colnames(tmp2.clr.t))
for(k in opened){    ## the opened regions (Enh, TSS)
 print(colnames(tmp2.clr.t)[k])
 tmp2.clr.t.direct[,k] = -tmp2.clr.t.direct[,k]
}
tmp$directional_deviation_sum = rowSums(tmp2.clr.t.direct)

# save rds
saveRDS(tmp,"Morabito.TSS1.erosion.score.metadata.cleaned.rds")

# color uMAP by erosion
proj <- loadArchRProject(path = "../Save-Proj2-TSS1") # read in proj for umap coordindate
meta = data.frame(proj@cellColData)

proj$IterativeLSI.UMAP_Dimension_1 = proj@embeddings$UMAP$df[,1]
proj$IterativeLSI.UMAP_Dimension_2 = proj@embeddings$UMAP$df[,2]
meta = data.frame(proj@cellColData)
all(rownames(meta) == rownames(tmp))
tmp$IterativeLSI.UMAP_Dimension_1 = meta$IterativeLSI.UMAP_Dimension_1
tmp$IterativeLSI.UMAP_Dimension_2 = meta$IterativeLSI.UMAP_Dimension_2
saveRDS(tmp,"Morabito.TSS1.erosion.score.metadata.cleaned.rds")

##### redraw umap
df <- as.data.frame(shuf(tmp))

png(paste0("Morabito.TSS",TSS,".by_erosion.clr_aross_all_ct.UMAP.cleaned.png"),width=30, height=12,units='in', res=200)
p.list = list()

# by cell type
df$Celltype <- as.factor(df$Celltype)
coord = c()
types<-levels(df$Celltype)   #factor(unique(tmp$Celltype1))
for (i in 1:length(types)){
        ind=which(df$Celltype == types[i])
        loc.x = median(df[ind,]$IterativeLSI.UMAP_Dimension_1)
        loc.y = median(df[ind,]$IterativeLSI.UMAP_Dimension_2)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

cls.cols1 = c("#E31A1C", "grey30","#33A02C", "#B2DF8A", "#CAB2D6", "#FDBF6F", "#B15928", "grey70")
tsp.col = function(x, alpha=0.5){
  rr = col2rgb(x)
  rgb(t(rr)/255, alpha=alpha)
}
cls.cols =  tsp.col(cls.cols1)

p.list[[1]] = ggplot(df,aes(x=IterativeLSI.UMAP_Dimension_1, y=IterativeLSI.UMAP_Dimension_2,colour=Celltype)) + geom_point(shape=16,size=0.5) + mytheme2 + ggtitle(paste0("TSS",TSS," by cell type")) + theme(legend.position = "none") + scale_colour_manual(values=cls.cols)

# by erosion score
p.list[[2]] = ggplot(df,aes(x=IterativeLSI.UMAP_Dimension_1, y=IterativeLSI.UMAP_Dimension_2,colour=directional_deviation_sum)) + geom_point(shape=16,size=0.5) + mytheme2 + ggtitle(paste0("TSS",TSS," by erosion score")) + theme(legend.position = "none") + scale_colour_gradientn(colours = c("cyan3","yellow", "darkred"))

# by erosion score
p.list[[3]] = ggplot(df,aes(x=IterativeLSI.UMAP_Dimension_1, y=IterativeLSI.UMAP_Dimension_2,colour=directional_deviation_sum)) + geom_point(shape=16,size=0.5) + mytheme2 + ggtitle(paste0("TSS",TSS," by erosion score")) + scale_colour_gradientn(colours = c("cyan3","yellow3", "darkred"))

library(gridExtra)
do.call("grid.arrange", c(p.list, nrow=1))
dev.off()

# erosion vs. AD status
summ =  data.frame(table(df$Sample,df$Celltype))
meta2 = unique(df[,c('Sample','diagnosis')])
summ2 = merge(summ,meta2,by.x='Var1',by.y='Sample')
library(plyr)
summ2 = ddply(summ2, "Var1", transform,
        Percent = Freq / sum(Freq) * 100)

summ2.y = summ2[summ2$Var2 == 'Deidentified',]
summ2.y = summ2.y[order(summ2.y$Percent),]

wilcox.out = wilcox.test(summ2.y[summ2.y$diagnosis=="AD",]$Percent,summ2.y[summ2.y$diagnosis=="Control",]$Percent,alternative="greater") # testing whether AD is higher
pdf('Morabito.deid.fraction.pdf',useDingbats=F)
summ2.y$diagnosis = factor(summ2.y$diagnosis,levels=c('Control','AD'))
ggplot(summ2.y,aes(x=diagnosis,y=Percent,color=diagnosis)) + geom_boxplot() +geom_jitter(position=position_jitter(0.2))  + 
 theme_classic() + scale_colour_manual(values=c('cyan3','pink4')) + ggtitle(paste0('Percentage of de-identified cells, p=',wilcox.out$p.value)) + theme(aspect.ratio=2.5)

mu.mean <- ddply(df, "Sample", summarise, grp.mean=mean(directional_deviation_sum))
mu.median <- ddply(df, "Sample", summarise, grp.median=median(directional_deviation_sum))
mu.mean = merge(mu.mean,meta2,by='Sample')
mu.median = merge(mu.median,meta2,by="Sample")
wilcox.out1 = wilcox.test(mu.mean[mu.mean$diagnosis=="AD",]$grp.mean,mu.mean[mu.mean$diagnosis=="Control",]$grp.mean,alternative="greater") # testing whether AD is higher
wilcox.out2 = wilcox.test(mu.median[mu.median$diagnosis=="AD",]$grp.median,mu.median[mu.median$diagnosis=="Control",]$grp.median,alternative="greater") # testing whether AD is higher

mu.mean$diagnosis = factor(mu.mean$diagnosis,levels=c('Control','AD'))
mu.median$diagnosis = factor(mu.median$diagnosis,levels=c('Control','AD'))

ggplot(mu.mean,aes(x=diagnosis,y=grp.mean,color=diagnosis)) + geom_boxplot() +geom_jitter(position=position_jitter(0.2)) +
 theme_classic() + scale_colour_manual(values=c('cyan3','pink4')) + ggtitle(paste0('Erosion score, p=',wilcox.out1$p.value)) + theme(aspect.ratio=2.5)
ggplot(mu.median,aes(x=diagnosis,y=grp.median,color=diagnosis)) + geom_boxplot() +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + scale_colour_manual(values=c('cyan3','pink4')) + ggtitle(paste0('Erosion score, p=',wilcox.out2$p.value)) + theme(aspect.ratio=2.5)


dev.off()

# use group defined by Carles; 0327
tmp = readRDS("Morabito.TSS1.erosion.score.metadata.cleaned.rds") # read in data
meta = read.table("../Carles.metadata.txt",head=T,sep="\t")
# no id provided in table s1, use pmi + age + sex to locate samples
meta$sex2 = "male"
meta$sex2[meta$Sex=="F"] = "female"
meta$idx = paste0(meta$PMI,"-",meta$Age,"-",meta$sex2)
meta = meta[,c("idx","group1","group2")]

tmp$idx = paste0(tmp$pmi,"-",tmp$AGE,"-",tmp$sex)

tmp = merge(tmp,meta)

# erosion vs. AD status
df = data.frame(tmp)
summ =  data.frame(table(df$Sample,df$Celltype))
meta2 = data.frame(unique(df[,c('Sample','group1')]))
summ2 = merge(summ,meta2,by.x='Var1',by.y='Sample')
library(plyr)
summ2 = ddply(summ2, "Var1", transform,
        Percent = Freq / sum(Freq) * 100)

summ2.y = summ2[summ2$Var2 == 'Deidentified',]
summ2.y = summ2.y[order(summ2.y$Percent),]

wilcox.out = wilcox.test(summ2.y[summ2.y$group1=="AD",]$Percent,summ2.y[summ2.y$group1=="Normal",]$Percent,alternative="greater") # testing whether AD is higher
pdf('Morabito.deid.fraction.Carles_group1.pdf',useDingbats=F)
ggplot(summ2.y,aes(x=group1,y=Percent,color=group1)) + geom_boxplot() +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Percentage of de-identified cells, p=',wilcox.out$p.value)) + theme(aspect.ratio=1)

mu.mean <- ddply(df, "Sample", summarise, grp.mean=mean(directional_deviation_sum))
mu.median <- ddply(df, "Sample", summarise, grp.median=median(directional_deviation_sum))
mu.mean = merge(mu.mean,meta2,by='Sample')
mu.median = merge(mu.median,meta2,by="Sample")
wilcox.out1 = wilcox.test(mu.mean[mu.mean$group1=="AD",]$grp.mean,mu.mean[mu.mean$group1=="Normal",]$grp.mean,alternative="greater") # testing whether AD is higher
wilcox.out2 = wilcox.test(mu.median[mu.median$group1=="AD",]$grp.median,mu.median[mu.median$group1=="Normal",]$grp.median,alternative="greater") # testing whether AD is higher

ggplot(mu.mean,aes(x=group1,y=grp.mean,color=group1)) + geom_boxplot() +geom_jitter(position=position_jitter(0.2)) +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out1$p.value)) + theme(aspect.ratio=1)
ggplot(mu.median,aes(x=group1,y=grp.median,color=group1)) + geom_boxplot() +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out2$p.value)) + theme(aspect.ratio=1)


dev.off()

df = data.frame(tmp)
summ =  data.frame(table(df$Sample,df$Celltype))
meta2 = data.frame(unique(df[,c('Sample','group2')]))
summ2 = merge(summ,meta2,by.x='Var1',by.y='Sample')
library(plyr)
summ2 = ddply(summ2, "Var1", transform,
        Percent = Freq / sum(Freq) * 100)

summ2.y = summ2[summ2$Var2 == 'Deidentified',]
summ2.y = summ2.y[order(summ2.y$Percent),]

wilcox.out = wilcox.test(summ2.y[summ2.y$group2=="AD",]$Percent,summ2.y[summ2.y$group2=="Normal",]$Percent,alternative="greater") # testing whether AD is higher
pdf('Morabito.deid.fraction.Carles_group2.pdf',useDingbats=F)
ggplot(summ2.y,aes(x=group2,y=Percent,color=group2)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Percentage of de-identified cells, p=',wilcox.out$p.value)) + theme(aspect.ratio=1)

mu.mean <- ddply(df, "Sample", summarise, grp.mean=mean(directional_deviation_sum))
mu.median <- ddply(df, "Sample", summarise, grp.median=median(directional_deviation_sum))
mu.mean = merge(mu.mean,meta2,by='Sample')
mu.median = merge(mu.median,meta2,by="Sample")
wilcox.out1 = wilcox.test(mu.mean[mu.mean$group2=="AD",]$grp.mean,mu.mean[mu.mean$group2=="Normal",]$grp.mean,alternative="greater") # testing whether AD is higher
wilcox.out2 = wilcox.test(mu.median[mu.median$group2=="AD",]$grp.median,mu.median[mu.median$group2=="Normal",]$grp.median,alternative="greater") # testing whether AD is higher

ggplot(mu.mean,aes(x=group2,y=grp.mean,color=group2)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2)) +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out1$p.value)) + theme(aspect.ratio=1)
ggplot(mu.median,aes(x=group2,y=grp.median,color=group2)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out2$p.value)) + theme(aspect.ratio=1)

dev.off()

### exclude the unclear samples
df = df[df$group1 %in% c("AD","Normal","Normal_stage2"),]
summ =  data.frame(table(df$Sample,df$Celltype))
meta2 = data.frame(unique(df[,c('Sample','group1')]))
summ2 = merge(summ,meta2,by.x='Var1',by.y='Sample')
library(plyr)
summ2 = ddply(summ2, "Var1", transform,
        Percent = Freq / sum(Freq) * 100)

summ2.y = summ2[summ2$Var2 == 'Deidentified',]
summ2.y = summ2.y[order(summ2.y$Percent),]

wilcox.out = wilcox.test(summ2.y[summ2.y$group1=="AD",]$Percent,summ2.y[summ2.y$group1=="Normal",]$Percent,alternative="greater") # testing whether AD is higher
pdf('Morabito.deid.fraction.Carles_group1_clean.pdf',useDingbats=F)
summ2.y$group1 = factor(summ2.y$group1,levels=c("Normal","AD"))
ggplot(summ2.y,aes(x=group1,y=Percent,color=group1)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Percentage of de-identified cells, p=',wilcox.out$p.value)) + theme(aspect.ratio=2.5) + scale_color_manual(values=c("cyan3","pink4"))

mu.mean <- ddply(df, "Sample", summarise, grp.mean=mean(directional_deviation_sum))
mu.median <- ddply(df, "Sample", summarise, grp.median=median(directional_deviation_sum))
mu.mean = merge(mu.mean,meta2,by='Sample')
mu.median = merge(mu.median,meta2,by="Sample")
wilcox.out1 = wilcox.test(mu.mean[mu.mean$group1=="AD",]$grp.mean,mu.mean[mu.mean$group1=="Normal",]$grp.mean,alternative="greater") # testing whether AD is higher
wilcox.out2 = wilcox.test(mu.median[mu.median$group1=="AD",]$grp.median,mu.median[mu.median$group1=="Normal",]$grp.median,alternative="greater") # testing whether AD is higher

mu.mean$group1 = factor(mu.mean$group1,levels=c("Normal","AD"))
mu.median$group1 = factor(mu.median$group1,levels=c("Normal","AD"))
ggplot(mu.mean,aes(x=group1,y=grp.mean,color=group1)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2)) +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out1$p.value)) + theme(aspect.ratio=2.5) + scale_color_manual(values=c("cyan3","pink4"))
ggplot(mu.median,aes(x=group1,y=grp.median,color=group1)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out2$p.value)) + theme(aspect.ratio=2.5)  + scale_color_manual(values=c("cyan3","pink4"))
dev.off()

summ =  data.frame(table(df$Sample,df$Celltype))
meta2 = data.frame(unique(df[,c('Sample','group2')]))
summ2 = merge(summ,meta2,by.x='Var1',by.y='Sample')
library(plyr)
summ2 = ddply(summ2, "Var1", transform,
        Percent = Freq / sum(Freq) * 100)

summ2.y = summ2[summ2$Var2 == 'Deidentified',]
summ2.y = summ2.y[order(summ2.y$Percent),]

wilcox.out = wilcox.test(summ2.y[summ2.y$group2=="AD",]$Percent,summ2.y[summ2.y$group2=="Normal",]$Percent,alternative="greater") # testing whether AD is higher
pdf('Morabito.deid.fraction.Carles_group2_clearn.pdf',useDingbats=F)
summ2.y$group2 = factor(summ2.y$group2,levels=c("Normal","Normal_stage2","AD"))
ad.color = c("cyan3","pink","pink4")
ggplot(summ2.y,aes(x=group2,y=Percent,color=group2)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Percentage of de-identified cells, p=',wilcox.out$p.value)) + theme(aspect.ratio=2.5) + scale_color_manual(values=ad.color)

mu.mean <- ddply(df, "Sample", summarise, grp.mean=mean(directional_deviation_sum))
mu.median <- ddply(df, "Sample", summarise, grp.median=median(directional_deviation_sum))
mu.mean = merge(mu.mean,meta2,by='Sample')
mu.median = merge(mu.median,meta2,by="Sample")
wilcox.out1 = wilcox.test(mu.mean[mu.mean$group2=="AD",]$grp.mean,mu.mean[mu.mean$group2=="Normal",]$grp.mean,alternative="greater") # testing whether AD is higher
wilcox.out2 = wilcox.test(mu.median[mu.median$group2=="AD",]$grp.median,mu.median[mu.median$group2=="Normal",]$grp.median,alternative="greater") # testing whether AD is higher

mu.mean$group2 = factor(mu.mean$group2,levels=c("Normal","Normal_stage2","AD"))
mu.median$group2 = factor(mu.median$group2,levels=c("Normal","Normal_stage2","AD"))
ggplot(mu.mean,aes(x=group2,y=grp.mean,color=group2)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2)) +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out1$p.value)) + theme(aspect.ratio=2.5) + scale_color_manual(values=ad.color)
ggplot(mu.median,aes(x=group2,y=grp.median,color=group2)) + geom_boxplot(outlier.shape=NA) +geom_jitter(position=position_jitter(0.2))  +
 theme_classic() + ggtitle(paste0('Erosion score, p=',wilcox.out2$p.value)) + theme(aspect.ratio=2.5) + scale_color_manual(values=ad.color) 
dev.off()



