#!/usr/bin/R
# --------------------------------------------------------------
# Nebula as in minimal example, but using the precomp. matrices:
# Updated: 05/25/21;
# From Carles; edited by Xushen: 08/22/21;
# --------------------------------------------------------------
.libPaths("/net/bmc-lab5/data/kellis/users/xiongxs/R/3.6")

library(tidyr)
library(Matrix)

# For QC; binarize matrix 
library(biclust)
library(dplyr)

# For nebula: 
library(nebula)
library(DESeq2)
library(RUVSeq)
library(qvalue)

# For plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
print(version)

## parameters for cell types
celltype = 'Ast'
path = "AD.path"
pctcut = 0.1
#NRUV = 10

library(argparser, quietly=TRUE)
p <- arg_parser("Cell type diff peak analysis;")
p <- add_argument(p, "ct", help="cell type for analysis")
p <- add_argument(p, "pctcut", help="percent of cells for each gene")
arg <- parse_args(p)
celltype=arg$ct
pctcut = arg$pctcut

for(NRUV in c(1,3,5,10)){
prefstr = paste0(celltype,'.',path,".pctcut",pctcut,".nRUV",NRUV,".non_vs_early")
# Directories:
#topimgdir = paste0(img, 'multiRegion/')
#plotdir = paste0(img, 'multiRegion/difftl/')
regdir = paste0('dereg/')
imgpref = paste0('difftl/')
#cmd = paste('mkdir -p', topimgdir, plotdir, regdir)
cmd = paste('mkdir -p',regdir,imgpref)
system(cmd)

# Building functions for regression:
asform = function(x){ as.formula(paste0(x, collapse='')) }

# Read in pseudo-gene matrix
mat.df <- readRDS(paste0(celltype,".mtx.count.rds"))
mat <-  assays(mat.df)$PeakMatrix
mat = as.matrix(mat)

pathdf <- mat.df@colData
rm(mat.df)
pathdf = pathdf[pathdf$Pathology != "lateAD",]
mat = mat[,rownames(pathdf)]
pathdf$barcode = rownames(pathdf)
pathdf$AD.path = 0
pathdf$AD.path[pathdf$Pathology=="earlyAD"] = 1

# -----------------------------------------------
# Run differential expression using Nebula + RUV:
# -----------------------------------------------
outtsv = paste0(regdir, 'nebula_ruv.', prefstr,'.tsv.gz')
outrda = paste0(regdir, 'nebula_ruv.', prefstr,'.rda')
nsigfile = paste0(regdir, 'nebula_ruv.', prefstr,'.nsig.tsv')
if (!file.exists(outrda)){
    # --------------
    # Subset matrix:
    # --------------

    # calculate pct of cells for each gene
    # binarize to get the pct of cells that express each gene
    mat.count.bi <- binarize(mat,threshold=0.5)
    pctcells = rowSums(mat.count.bi)/ncol(mat)
    
    keep.genes = names(pctcells)[pctcells > pctcut]
    mat = mat[keep.genes,]
    print(paste("[STATUS] Subsetting matrix to", 
                paste0(dim(mat), collapse = ' x '),'(gene x cell)'))

    # --------------------------------------------
    # Run RUV on the data at the individual level:
    # --------------------------------------------
    # Make the individual-aggregate matrix:
#    NRUV = 10
    print(paste0("[STATUS] Running RUV for N=", NRUV))
    
    # Carles use matrix multiplication for aggregation;  
    if(0){
    pids = as.character(unique(pathdf$projid))
    tform = make.tform(pathdf$projid, u=pids)
    indvar = 'projid'
    data_ind = mat %*% tform 
    }

    # I don't have the make.tform function; use dplyr instead
    mat.t = as.data.frame(t(mat))
    indvar = 'projid'
    all(rownames(mat.t)==rownames(pathdf))  # check cell order
    mat.t$projid = pathdf$projid
    mat.t.agg = mat.t %>% 
	  group_by(projid) %>% 
	  summarise(across(everything(), sum))
    mat.t.agg = as.data.frame(mat.t.agg)
    rownames(mat.t.agg) = mat.t.agg$projid
    data_ind = as.matrix(t(mat.t.agg[,-1]))

    # Make the aggregate design matrix:
    uqcols = c(indvar,path)
    dform = asform(c('~',path))
    uqobs = unique(pathdf[,uqcols])
    rownames(uqobs) = uqobs[[indvar]]
    uqobs = uqobs[colnames(data_ind),]
    design = model.matrix(dform, data=uqobs)

    # DESeq2 object
    d_e = DGEList(data_ind, genes=rownames(data_ind))
    keep = rowSums(cpm(d_e)>1) >= 3
    d_e = d_e[keep, , keep.lib.sizes=FALSE]
    d_e = calcNormFactors(d_e, method="TMM")
    d_e = estimateGLMCommonDisp(d_e, design)
    d_e = estimateGLMTagwiseDisp(d_e, design)
    fit1 = glmFit(d_e, design)
    res1 = residuals(fit1, type="deviance")
    ruv_cov = RUVr(round(d_e$counts), 
                    as.character(rownames(d_e$counts)), 
                    k=NRUV, res1)

    # Merge the learned factors back into the data.frame:
    uqobs = cbind(uqobs, ruv_cov$W)
    pathdf = merge(pathdf, uqobs,by=uqcols, all.x=TRUE)
    # Re-order to original
    rownames(pathdf) = pathdf$barcode
    pathdf = pathdf[colnames(mat),]

    # -----------------------------------------------
    # Run NEBULA with the RUV results and other vars:
    # -----------------------------------------------
    ruvw = paste0("W_", 1:NRUV)
    flist = c('~', path)
    flist = c(flist, ' + age_death + msex + pmi')
    if(0){flist = c(flist, ' + TSSEnrichment')}  # don't use TSS ennrichment
#    flist = c(flist, ' + cpg + nFeature_RNA')  ## what's cpg?
    # flist = c(flist, ' + pctMT')
    flist = c(flist, " + ", paste(ruvw, collapse=" + "))
    nb.form = asform(flist)
    print(nb.form)
    mdx = model.matrix(nb.form, data=pathdf)
    
    if (path %in% c('cogdxad','nrad')){
        pathstr = paste0(path, 'AD')
    } else { 
        pathstr = path 
    }
    leff = paste0('logFC_', pathstr)
    peff = paste0('p_', pathstr)

    chunksize=500
    nchunk = ceiling(nrow(mat) / chunksize)
    fulldf = c()
    offset = log10(pathdf$nFrags)
    for (chunk in 1:nchunk){
        print(chunk)
        ind = (1 + (chunk-1) * chunksize):min(c(chunk * chunksize, nrow(mat))) 
        submat = mat[ind,]
        data_g = group_cell(count=submat,id=as.character(pathdf[[indvar]]),pred=mdx)
	re = nebula(data_g$count,data_g$id,pred=data_g$pred,offset=offset, model='PMM')
        rdf = re$summary
        resdf = rdf[order(rdf[[peff]]),c('gene',leff, peff)]
        names(resdf) = c('gene','logFC','p')
        fulldf = rbind(fulldf, resdf)
    }
    fulldf = fulldf[order(fulldf$p),]
    fulldf$padj = p.adjust(fulldf$p, 'fdr')
#    fulldf$q = qvalue(fulldf$p)$q

    pcut = 0.1
    fulldf$col = 1 * (fulldf$padj < pcut) * (2 - 1 * (fulldf$logFC < 0))   ## label significant up as 2, down as 1, non-signif as 0
    fulldf$pc = pctcells[fulldf$gene]
    labdf = fulldf[fulldf$col != 0,]

    # Write out the regression result dataframes:
    write.table(fulldf, gzfile(outtsv), quote=F, row.names=F, sep="\t")
    save(resdf, fulldf, file=outrda)

    pcols = brewer.pal(12,'Paired')
    gplot = ggplot(fulldf, aes(logFC, -log10(p), color=factor(col))) + 
        geom_point(cex=.25) + 
        geom_text_repel(data=labdf, aes(logFC, -log10(p), label=gene, color=factor(col)), size=2, max.overlaps=20) + 
        scale_color_manual(values=c('grey80',pcols[1],pcols[5])) + 
        scale_y_continuous(expand=c(0,0)) + 
        theme_pubr() + theme(legend.position='none')
    ggsave(paste0(imgpref, 'volcano_',prefstr,'_nebula_frommat.png'),gplot, units='in', dpi=450, width=6, height=6)
    ggsave(paste0(imgpref, 'volcano_',prefstr,'_nebula_frommat.pdf'),gplot, units='in', dpi=450, width=6, height=6)

    # Write out nsig: 
    ctvec = c(celltype=celltype, path=path, pctcut=pctcut, pcut=paste0("fdr_",pcut), nc = nrow(fulldf[fulldf$col==0,]),down=nrow(fulldf[fulldf$col==1,]),up=nrow(fulldf[fulldf$col==2,]))
    nsig = ctvec

    pcut = 0.05
    fulldf$col = 1 * (fulldf$p < pcut) * (2 - 1 * (fulldf$logFC < 0))   ## label significant up as 2, down as 1, non-signif as 0
    fulldf$pc = pctcells[fulldf$gene]

    ctvec = c(celltype=celltype, path=path, pctcut=pctcut, pcut=paste0("pval_",pcut), nc = nrow(fulldf[fulldf$col==0,]),down=nrow(fulldf[fulldf$col==1,]),up=nrow(fulldf[fulldf$col==2,]))
    nsig = rbind(nsig,ctvec)

    write.table(nsig, nsigfile, quote=F, row.names=F, sep="\t")
} else {
    print("[STATUS] Regression output files already exist")
}

}

