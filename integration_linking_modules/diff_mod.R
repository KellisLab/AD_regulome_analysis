#!/usr/bin/env Rscript

library(rGREAT)
library(GenomicRanges)
run.sce <- function(sce, func=deg.deseq2) {
    dle = func(sce, individual="Sample", pathology="Pathology", case="lateAD", control="earlyAD", covariates=c("pmi", "nFrags", "msex", "age_death"));
    dln = func(sce, individual="Sample", pathology="Pathology", case="lateAD", control="nonAD", covariates=c("pmi", "nFrags", "msex", "age_death"));
    den = func(sce, individual="Sample", pathology="Pathology", case="earlyAD", control="nonAD", covariates=c("pmi", "nFrags", "msex", "age_death"));
    DE = dplyr::bind_rows(list("late_vs_early"=dle, "late_vs_non"=dln, "early_vs_non"=den), .id="comparison");
    return(DE)
}

load.all <- function(dir="~/data-group/Benjamin/AD_ATAC/modules") {
    FL = list.files(dir, pattern="^agg_.*h5ad$", full.names=T)
    names(FL) = gsub("^agg_", "", gsub(".h5ad$", "", basename(FL)))
    return(lapply(FL, read.anndata.layer))
}

run.rgreat <- function(dir="~/data-group/Benjamin/AD_ATAC/modules") {
    FL = list.files(dir, pattern="^.*.bed$", full.names=T)
    names(FL) = gsub(".bed$", "", basename(FL))
    varlist = lapply(FL, function(x) {
        with(read.table(x, sep="\t"), GRanges(seqnames=V1, ranges=IRanges(V2, V3), cls=V4))
    })
    G = lapply(varlist, function(bg) {
        print("Celltype")
        permod = lapply(split(bg, bg$cls), function(fg) {
            print(fg)
            job = rGREAT::submitGreatJob(fg, bg, species="hg38")
            L = rGREAT::getEnrichmentTables(job, download_by="tsv")
            if (is.list(L)) {
                return(Reduce(rbind, L))
            } else {
                return(L)
            }
        })
        return(dplyr::bind_rows(permod, .id="leiden"))
    })
    return(dplyr::bind_rows(G, .id="CellType"))
}
LA = load.all()
df = dplyr::bind_rows(lapply(LA, function(x) { run.sce(x)}), .id="CellType")
df = df[order(df$FDR),]
rownames(df) = NULL

plot.de.mod <- function(go, dam) {
    colnames(dam) = gsub("^gene$", "leiden", colnames(dam))
    mat = dam[dam$FDR < 0.05,c("leiden", "CellType")]
    mat = mat[!duplicated(mat),]
    mat = mat[order(mat$leiden),]
    mat = mat[order(mat$CellType),]
    rownames(mat) = with(mat, paste0(CellType, "_", leiden))
    go = go[go$Ontology == "GO Biological Process",]
    gotop = go[!duplicated(go[c("leiden", "CellType")]), c("leiden", "CellType", "Desc", "HyperFdrQ")]
    rownames(gotop) = with(gotop,paste0(CellType, "_", leiden))
    mat$Desc = gotop[rownames(mat),]$Desc
    mat$HyperFdrQ = gotop[rownames(mat),]$HyperFdrQ
    rownames(dam) = paste0(dam$CellType, "_", dam$leiden, "_", dam$comparison)
    for (cond in c("early_vs_non", "late_vs_non", "late_vs_early")) {
        mat[[paste0(cond, "_log2FC")]] = dam[paste0(rownames(mat), "_", cond),"log2FC"]
        mat[[paste0(cond, "_FDR")]] = dam[paste0(rownames(mat), "_", cond),"FDR"]
    }
    FC = as.matrix(mat[grep("_log2FC$", colnames(mat))])
    PV = as.matrix(mat[grep("_FDR$", colnames(mat))])
    ht = Heatmap(FC, row_gap=unit(4, "mm"), width=unit(3, "cm"), heatmap_legend_param=list(title="log2FC"), column_title="Differentially accessible modules in Alzheimer's", cluster_columns=F, show_row_dend=F, row_split=mat$CellType, cell_fun=function(j, i, x, y, w, h, fill) {
        if (PV[i, j] < 0.001) {
            grid.text("***", x, y, gp=gpar(fontsize=8))
        } else if (PV[i, j] < 0.01) {
            grid.text("**", x, y, gp=gpar(fontsize=8))
        } else if (PV[i, j] < 0.05) {
            grid.text("*", x, y, gp=gpar(fontsize=8))
        }
    }) + rowAnnotation(modName=anno_text(paste0(mat$leiden), gp=gpar(fontsize=7)), GREAT_mlog10_HyperFdrQ=-log10(mat$HyperFdrQ), topGo=anno_text(mat$Desc, gp=gpar(fontsize=7)))
}
