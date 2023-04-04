# look for positive TF-Regulators; instead of TF activity
.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")

library(ArchR)
set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")   ## use hg38
library(stringr)

proj <- loadArchRProject(path = "B.Save-Proj3")

###### correlation between TF expression vs. deviation score
# positive regulators at the cell type level
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Celltype4")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
seZ.mat = assays(seZ)$MotifMatrix
rownames(seZ.mat) = str_extract(rowData(seZ)$name,"[a-zA-Z0-9]+")
ct.names = colnames(seZ.mat)
seZ.mat = data.frame(seZ.mat)
seZ.mat$ct.max = ct.names[apply(seZ.mat,1,which.max)]

# get gene score
genescore <- getGroupSE(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "Celltype4")
GS.mat = assays(genescore)$GeneScoreMatrix
rownames(GS.mat) = rowData(genescore)$name
GS.mat = data.frame(GS.mat)
ct.names = colnames(GS.mat)
GS.mat = data.frame(GS.mat)
GS.mat$ct.max = ct.names[apply(GS.mat,1,which.max)]

# order the two matrices
tf.list = intersect(rownames(GS.mat),rownames(seZ.mat))

seZ.mat = seZ.mat[tf.list,]
GS.mat = GS.mat[tf.list,]

# do correlation
corGSM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM = data.frame(corGSM_MM)
rownames(corGSM_MM) = corGSM_MM$GeneScoreMatrix_name
corGSM_MM = corGSM_MM[tf.list,]
corGSM_MM = corGSM_MM[order(corGSM_MM$padj),]

seZ.mat = seZ.mat[rownames(corGSM_MM),]
GS.mat = GS.mat[rownames(corGSM_MM),]
write.table(data.frame(TF=rownames(seZ.mat),seZ.mat),'C.16.Positive.TF/TF.chromVar.matrix.xls',row.names=F,sep="\t",quote=F)
write.table(data.frame(TF=rownames(GS.mat),GS.mat),'C.16.Positive.TF/TF.inferred_exp.matrix.xls',row.names=F,sep="\t",quote=F)

corGSM_MM = corGSM_MM[,c('GeneScoreMatrix_matchName','MotifMatrix_name','cor','padj','pval')]

corGSM_MM$Strongest.ChromVar = seZ.mat$ct.max
corGSM_MM$Strongest.Expression = GS.mat$ct.max

corGSM_MM$Is.Cand.Reglator = 'No'
corGSM_MM$Is.Cand.Reglator[corGSM_MM$padj<0.01 & corGSM_MM$cor>0.6] = "Yes"

write.table(corGSM_MM,'C.16.Positive.TF/TF.motif_vs_exp.correlation.xls',row.names=F,sep="\t",quote=F)

## do some plots!
# the matrix table first!
corGSM_MM.cand = corGSM_MM[corGSM_MM$Is.Cand.Reglator=='Yes',]

########## visualize top correlated TFs
tf.names.df = data.frame(gene=corGSM_MM.cand$GeneScoreMatrix_matchName,label=corGSM_MM.cand$MotifMatrix_name)
cor.sort.pos = tf.names.df

library(RColorBrewer)
test.color = colorRampPalette(brewer.pal(n = 3, name = "OrRd" ))(30)
n_each = 4  # number of columns & rows in each plot
n_total = n_each * n_each
pdf("C.16.Positive.TF/AllCell.chromVar.UMAP.top.correlated.positive.deviation.Z.20220917.pdf")
for(i in 1:5){
 cor.sort.pos.top = cor.sort.pos[((i-1)*n_total+1):(i*n_total),]
 tf.plot = cor.sort.pos.top$label
 tf.plot = paste0("z:",tf.plot)
 p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "MotifMatrix",
    name = tf.plot,
    embedding = "UMAP",
    pal = test.color,
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
 p3 = do.call(cowplot::plot_grid, c(list(ncol = n_each),p2))
 print(p3)
}
dev.off()

## plot the corresponding TF exp umap
pdf("C.16.Positive.TF/AllCell.chromVar.UMAP.top.correlated.positive.Exp.20220917.pdf")
for(i in 1:5){
 cor.sort.pos.top = cor.sort.pos[((i-1)*n_total+1):(i*n_total),]
 tf.plot = as.character(cor.sort.pos.top$gene)
 p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = tf.plot,
    embedding = "UMAP",
    pal = test.color,
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
 p3 = do.call(cowplot::plot_grid, c(list(ncol = n_each),p2))
 print(p3)
}
dev.off()

## side by side heatmap
seZ.mat = seZ.mat[rownames(corGSM_MM.cand),]
GS.mat = GS.mat[rownames(corGSM_MM.cand),]
cts.list = c("Ex","In","Ast","Oligo","OPC","Microglia","PerEndo")
seZ.mat = seZ.mat[,c(cts.list,'ct.max')]
GS.mat = GS.mat[,c(cts.list,'ct.max')]
saveRDS(seZ.mat,'C.16.Positive.TF/motif.mat.rds')
saveRDS(GS.mat,'C.16.Positive.TF/Genescore.mat.rds')
saveRDS(corGSM_MM.cand,'C.16.Positive.TF/corGSM_MM.cand.rds')

library(pheatmap)
library(RColorBrewer)

pdf('Heatmap.of.Positive.TFs.pdf',width=5,useDingbats=F)
corGSM_MM.cand2 = corGSM_MM.cand[corGSM_MM.cand$cor>0.6,]
tf.list = corGSM_MM.cand2$GeneScoreMatrix_matchName

seZ.mat = seZ.mat[tf.list,c(1:7)]
reord<-order(apply(seZ.mat,1,which.max),decreasing = F)
seZ.mat.ord = seZ.mat[reord,]
pheatmap(seZ.mat.ord,scale='row',color = colorRampPalette(c('white','gray90','darkgreen'))(100),cluster_row = F,cluster_cols = F,fontsize = 12,border_color=NA,main='Motif enrichment')

GS.mat = GS.mat[tf.list,c(1:7)]
reord<-order(apply(GS.mat,1,which.max),decreasing = F)
GS.mat.ord = GS.mat[rownames(seZ.mat.ord),]
pheatmap(GS.mat.ord,scale='row',color = colorRampPalette(c('white','gray90','darkred'))(100),cluster_row = F,cluster_cols = F,fontsize = 12,border_color=NA,main='Inferred expression')

dev.off()

###### further consider RNA-ATAC correlation
.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")
library(Seurat)

# read in RNA-seq data and annotation
library(SingleCellExperiment)
rna = readRDS("/net/bmc-lab5/data/kellis/users/xiongxs/Proj_snATAC_AD/20210205_PFC.reprocess/2.UseArchR.Re2/C.01.rna.marker/RNA.subset.seurat.rds")

# read in the GS-chromVar correlation
cor = read.table('C.16.Positive.TF/TF.motif_vs_exp.correlation.xls',head=T)
tf.list = as.character(cor$GeneScoreMatrix_matchName)

Markers = FetchData(rna,tf.list)   # no 'BHLHA15'; remove it from the other two matrices
meta.rna = readRDS('./scRNA.meta.0816.rds')
meta.rna = meta.rna[rownames(Markers),]
all(rownames(meta.rna) == rownames(Markers))
Markers.mean = aggregate(Markers,list(meta.rna$major.celltype),FUN=mean)
rownames(Markers.mean) = Markers.mean$Group.1
Markers.mean.rna = Markers.mean[,-1]
Markers.mean.rna = Markers.mean.rna[c('Exc','Inh','Ast','Oli','Opc','Mic','Vas'),]

# organize the atac-GS matrix and chromVar matrix, and remove the genes not in the RNA genes
seZ.mat2 = seZ.mat[colnames(Markers.mean.rna),c('Ex','In','Ast','Oligo','OPC','Microglia','PerEndo')]
GS.mat2 = GS.mat[colnames(Markers.mean.rna),c('Ex','In','Ast','Oligo','OPC','Microglia','PerEndo')]
# reorder the cell types
GS.mat2.t = t(GS.mat2)
GS.mat2.t = GS.mat2.t[c('Ex','In','Ast','Oligo','OPC','Microglia','PerEndo'),]
seZ.mat2.t = t(seZ.mat2)
seZ.mat2.t = seZ.mat2.t[c('Ex','In','Ast','Oligo','OPC','Microglia','PerEndo'),]

atac_rna = data.frame()
for(tf in rownames(GS.mat2)){
 print(tf)
 out.cor = cor.test(GS.mat2.t[,tf],Markers.mean.rna[,tf])
 df0 = data.frame(TF=tf,cor=as.numeric(out.cor$estimate),pval=out.cor$p.value)
 atac_rna = rbind(atac_rna,df0)
}
atac_rna$padj = p.adjust(atac_rna$pval,method='fdr')

# attach to the correlation file
atac_chromvar_rna = merge(cor,atac_rna,by.x='GeneScoreMatrix_matchName',by.y='TF',suffix=c('_atac_chromVar','_atac_rna'))

# do filtering based on three axes
atac_chromvar_rna$Is.Cand.Reglator = 'No'
atac_chromvar_rna$Is.Cand.Reglator[atac_chromvar_rna$cor_atac_chromVar>0.5 & atac_chromvar_rna$padj_atac_chromVar<0.01 & atac_chromvar_rna$cor_atac_rna>0.5 & atac_chromvar_rna$padj_atac_rna<0.1] = 'Yes' # the atac vs rna correlation test is across 7 cell types only, therefore use a looser cutoff

## do plotting!
# first sort the matrix and make supp tables
atac_chromvar_rna = atac_chromvar_rna[order(atac_chromvar_rna$cor_atac_chromVar,decreasing=T),]
seZ.mat2 = seZ.mat[as.character(atac_chromvar_rna$GeneScoreMatrix_matchName),]
GS.mat2 = GS.mat[as.character(atac_chromvar_rna$GeneScoreMatrix_matchName),]
rna.mat2 = t(Markers.mean.rna[,as.character(atac_chromvar_rna$GeneScoreMatrix_matchName)])
write.table(data.frame(TF=rownames(seZ.mat2),seZ.mat2),'C.16.Positive.TF/TF.chromVar.matrix.0919.xls',row.names=F,sep="\t",quote=F)
write.table(data.frame(TF=rownames(GS.mat2),GS.mat2),'C.16.Positive.TF/TF.inferred_exp.matrix.0919.xls',row.names=F,sep="\t",quote=F)
write.table(data.frame(TF=rownames(rna.mat2),rna.mat2),'C.16.Positive.TF/TF.rna_gene_exp.matrix.0919.xls',row.names=F,sep="\t",quote=F)
write.table(atac_chromvar_rna,'C.16.Positive.TF/ATAC_RNA_ChromVar.cor.0919.xls',row.names=F,sep="\t",quote=F)

## side by side heatmap
library(pheatmap)
pdf('Heatmap.of.Positive.TFs.0919.pdf',width=3,useDingbats=F)
tf.list = as.character(atac_chromvar_rna[atac_chromvar_rna$Is.Cand.Reglator=='Yes',]$GeneScoreMatrix_matchName)
seZ.mat = seZ.mat2[tf.list,1:7]
reord<-order(apply(seZ.mat,1,which.max),decreasing = F)
seZ.mat.ord = seZ.mat[reord,]
pheatmap(seZ.mat.ord,scale='row',color = colorRampPalette(c('white','gray90','darkgreen'))(100),cluster_row = F,cluster_cols = F,fontsize = 12,border_color=NA,main='Motif enrichment')

GS.mat = GS.mat2[tf.list,c(1:7)]
GS.mat.ord = GS.mat[rownames(seZ.mat.ord),]
pheatmap(GS.mat.ord,scale='row',color = colorRampPalette(c('white','gray90','darkred'))(100),cluster_row = F,cluster_cols = F,fontsize = 12,border_color=NA,main='Inferred expression')

rna.mat = rna.mat2[tf.list,c(1:7)]
rna.mat.ord = rna.mat[rownames(seZ.mat.ord),]
pheatmap(rna.mat.ord,scale='row',color = colorRampPalette(c('white','gray90','cyan4'))(100),cluster_row = F,cluster_cols = F,fontsize = 12,border_color=NA,main='RNA expression')

dev.off()








