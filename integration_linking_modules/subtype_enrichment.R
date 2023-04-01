
FL = list.files("~/data-group/Benjamin/AD_ATAC/modules/", pattern="^agg_.*h5ad$", full.names=T)
names(FL) = benj::msub(c("^agg_", ".h5ad$"), c("", ""), basename(FL))
ct = read.table("~/data-group/Benjamin/AD_ATAC/TSS6/cell_type_predicted.tsv.gz",sep="\t",row.names=1, comment="", header=T)
adl = lapply(FL, function(fname) {
    print(fname)
    adata = benj::read_h5ad(fname, layer="raw")
    SummarizedExperiment::colData(adata)$cthr = ct$pred_cell_type_high_resolution[match(colnames(adata),rownames(ct))]
    return(adata)
})


deg = dplyr::bind_rows(lapply(adl, function(adata) {
    adata = se_normalize_total(adata, assay="raw")
    X = as.matrix(SummarizedExperiment::assays(adata)$normalized)
    cthr = as.factor(SummarizedExperiment::colData(adata)$cthr)
    print(str(cthr))
    X = log1p(X)
    ucthr = levels(cthr)
    out = lapply(ucthr, function(ct) {
        gwt = lapply(1:nrow(X), function(i) {
            wilcox.test(X[i,ct == cthr], X[i,ct != cthr], alternative="greater")
        })
        lwt = lapply(1:nrow(X), function(i) {
            wilcox.test(X[i,ct == cthr], X[i,ct != cthr], alternative="less")
        })
        gp = sapply(gwt, function(x) { x$p.value })
        lp = sapply(lwt, function(x) { x$p.value })
        gstat = sapply(gwt, function(x) { x$statistic })
        lstat = sapply(lwt, function(x) { x$statistic })
        odf =data.frame(p.value=apply(cbind(gp, lp), 1, min),
                        module=rownames(X),
                        subtype=ct)
        odf$statistic = ifelse(odf$p.value == gp, gstat, lstat)
        odf$direction = ifelse(odf$p.value == gp, "greater", "less")
        return(odf)
    })
    return(Reduce(rbind, out))
}), .id="CellType")
