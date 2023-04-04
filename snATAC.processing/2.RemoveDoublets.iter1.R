## first round of doublet removal at the cluster level 
# remove the clusters that show marker signatures of two cell types, and then re-do clustering + cell type assignment
.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38

proj <- loadArchRProject(path = "Save-Proj6")

# remove doublets
proj = proj[!(proj$Clusters %in% c("C6","C9","C10","C12","C20"))]
proj$Clusters.original = proj$Clusters
proj$UMAP1.original = proj@embeddings$UMAP$df[,1]
proj$UMAP2.original = proj@embeddings$UMAP$df[,2]

# redo clustering
ncell=length(proj$cellNames)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", iterations = 4, clusterParams = list(resolution = c(0.2), sampleCells = ncell, n.start = 10), name = "IterativeLSI",varFeatures = 50000,force = TRUE)

proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force = TRUE ,method = "Seurat")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force = TRUE)

saveRDS(proj,"ArchR.TSS6.1st_doublets_rm.2nd_clustering.rds")

proj = readRDS("ArchR.TSS6.1st_doublets_rm.2nd_clustering.rds")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
pdf("ArchR.TSS6.1st_doublets_rm.2nd_clustering.pdf")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
dev.off()

proj$IterativeLSI.UMAP_Dimension_1 = proj@embeddings$UMAP$df[,1]
proj$IterativeLSI.UMAP_Dimension_2 = proj@embeddings$UMAP$df[,2]

# plot QC traits
mytheme2 <-theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio = 1,axis.text = element_text(face="bold", size=rel(1)),axis.title = element_text(face="bold", size=rel(1)),legend.text=element_text(face="bold",size=10),plot.title=element_text(face="bold",size=rel(1)),legend.title =element_text(face="bold", size=rel(1)))

df = as.data.frame(proj@cellColData)

pdf("ArchR.TSS6.1st_doublets_rm.2nd_clustering.cluster.QC.pdf")
ggplot(df,aes(y=Clusters,x=TSSEnrichment)) + geom_boxplot() + mytheme2
ggplot(df,aes(y=Clusters,x=nFrags)) + geom_boxplot() + mytheme2
ggplot(df,aes(y=Clusters,x=DoubletScore)) + geom_boxplot() + mytheme2
ggplot(df,aes(y=Clusters,x=log(DoubletScore+2))) + geom_boxplot() + mytheme2
ggplot(df,aes(y=Clusters,x=DoubletEnrichment)) + geom_boxplot() + mytheme2
ggplot(df,aes(y=Clusters,x=log(DoubletEnrichment+2))) + geom_boxplot() + mytheme2
ggplot(df,aes(y=Clusters,x=BlacklistRatio)) + geom_boxplot() + mytheme2
dev.off()

p.list = list()
cls.list = sort(unique(df$Clusters))

for(i in 1:length(cls.list)){
 cls = cls.list[[i]]
 df$cls = "Others"
 df$cls[df$Clusters == cls] = cls
 df = df[order(df$cls,decreasing=T),]
 p.list[[i]] = ggplot(df,aes(x=IterativeLSI.UMAP_Dimension_1, y=IterativeLSI.UMAP_Dimension_2,colour=cls)) + geom_point(shape=16,size=0.5) + mytheme2 + ggtitle(cls) + scale_colour_manual(values=c("darkred","grey"))
}

library(gridExtra)
png("ArchR.TSS6.1st_doublets_rm.2nd_clustering.by.cluster.png",width=28, height=28,units='in', res=300)
do.call("grid.arrange", c(p.list, nrow=5))
dev.off()

# density
library(hexbin)
left = min(df$IterativeLSI.UMAP_Dimension_1) * 1.05
right = max(df$IterativeLSI.UMAP_Dimension_1) * 1.05
top = max(df$IterativeLSI.UMAP_Dimension_2) * 1.05
bottom = min(df$IterativeLSI.UMAP_Dimension_2) * 1.05

for(i in 1:length(cls.list)){
 cls = cls.list[[i]]
 df.cls = df[df$Clusters == cls,]
 p.list[[i]] = ggplot(df.cls,aes(x=IterativeLSI.UMAP_Dimension_1, y=IterativeLSI.UMAP_Dimension_2)) + geom_point(shape=16,size=0.5) + mytheme2 + ggtitle(cls) + stat_binhex() + scale_fill_gradient(low="lightblue", high="red") + xlim(left,right) + ylim(bottom,top)
}

library(gridExtra)
png(paste0("ArchR.TSS6.1st_doublets_rm.2nd_clustering.by.cluster.density.png"),width=28, height=28,units='in', res=300)
do.call("grid.arrange", c(p.list, nrow=5))
dev.off()

## marker
markers<-read.table("/home/xiongxs/data-lab5/Proj_snATAC_AD/CellMarker/FromCarles/known_markers.broad_ct.merge.uniq.tsv",header=T)
markerGenes<-as.character(markers$symbol)

####### visualizing marker gene intensity
proj <- addImputeWeights(proj)

p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
system("mkdir -p B.3.test")
pdf(paste0("B.3.test/ArchR.TSS6.1st_doublets_rm.2nd_clustering.MarkerGene.All.plot.umap.pdf"))
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
dev.off()






