
df = readRDS("~/tmp/modules/go.rds")
df = df[(df$HyperFdrQ < 0.01) & (df$FgRegionsHit >= 50) & (df$NumFgGenesHit >= 5) & (df$NumBgGenesHit <= 500),]
M = lapply(split(df, df$CellType), function(mf) {
    peaks = strsplit(mf$FgRegionNames, ",")
    upeaks = unique(Reduce(c, peaks))
    M = Reduce(rbind, lapply(seq_along(peaks), function(i) {
        data.frame(row=i, col=match(peaks[[i]], upeaks))
    }))
    M$GOID = as.factor(mf$ID[M$row])
    M = with(M, Matrix::sparseMatrix(i=as.integer(GOID), j=col, dimnames=list(levels(GOID), upeaks)))
    return(as.matrix(M))
})


se = lapply(M, function(mat) {
    print("motif overlap")
    gr = benj::parse.range(colnames(mat))
    se = benj::motif_overlap(gr, JASPAR2022::JASPAR2022)
    for (rn in rownames(mat)) {
        SummarizedExperiment::rowData(se)[[rn]] = ifelse(c(mat[rn,rownames(se)]),
                                                         rep(rn, nrow(se)),
                                                         rep("NONE", nrow(se)))
    }
    return(se)
})

go2desc = df[!duplicated(df$ID),c("ID", "Desc")]
rownames(go2desc) = go2desc$ID
L = lapply(se, function(xe) { Reduce(rbind, lapply(colnames(rowData(xe)), function(x) { benj::enrich_se(xe, x) })) })
ef = dplyr::bind_rows(L, .id="CellType")
ef = ef[order(-ef$mlog10p),]
ef$Desc = go2desc[ef$ID,]$Desc
