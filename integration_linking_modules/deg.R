#!/usr/bin/env Rscript
library(rhdf5)
library(Matrix)
library(SingleCellExperiment)
library(nebula)
library(RUVSeq)
library(DESeq2) ### for DGEList
library(edgeR)

read.h5ad.dataframe <- function(h5, sub="obs") {
    lst = rhdf5::h5read(h5, sub)
    df = data.frame(row.names=lst[["_index"]])
    for (cn in setdiff(names(lst), c("_index"))) {
        if (is.list(lst[[cn]])) {
            df[[cn]] = factor(lst[[cn]]$categories[lst[[cn]]$codes + 1],
                              levels=lst[[cn]]$categories)
        } else {
            df[[cn]] = lst[[cn]]
        }
    }
    return(df)
}

read.h5ad.raw <- function(h5) {
    hlist = h5ls(h5)
    obs = read.h5ad.dataframe(h5)
    if (("/raw/X" %in% hlist$group) & ("/raw/var" %in% hlist$group)) {
        M = h5read(h5, "/raw/X/", bit64conversion="bit64")
        if (is.list(M)) {
            M$col = as.integer(rep(-1, length(M$indices)))
            for (i in 1:(length(M$indptr)-1)) {
### Add one to begin as R is 1-indexed
                begin = M$indptr[i] + 1
                end = M$indptr[i+1]
                M$col[begin:end] = i ## index already starts at 1
            }
            M = Matrix::sparseMatrix(i=c(M$indices+1), j=c(M$col), x=c(M$data))
        }
        dimnames(M) = list(as.character(h5read(h5, "/raw/var/_index")),
                           as.character(h5read(h5, "/obs/_index")))
        return(SingleCellExperiment(list(counts=M), colData=obs))
    } else {
        print("Using /X")
        M = h5read(h5, "/X/", bit64conversion="bit64")
        if (is.list(M)) {
            M$col = as.integer(rep(-1, length(M$indices)))
            for (i in 1:(length(M$indptr)-1)) {
### Add one to begin as R is 1-indexed
                begin = M$indptr[i] + 1
                end = M$indptr[i+1]
                M$col[begin:end] = i ## index already starts at 1
            }
            M = Matrix::sparseMatrix(i=c(M$indices+1), j=c(M$col), x=c(M$data))
        }
        dimnames(M) = list(as.character(h5read(h5, "/var/_index")),
                           as.character(h5read(h5, "/obs/_index")))
        return(SingleCellExperiment(list(counts=M), colData=obs))
    }
}

make_pseudobulk <- function(x, u=NULL) {
    if (is.factor(x)) {
        u = levels(x)
    } else if (is.null(u)) {
        u = unique(x)
    }
    names(u) = u
    return(sapply(u, function(y) {
        1 * (x == y)
    }))
}
asform <- function(x) {
    as.formula(paste0(x, collapse=''))
}

deg <- function(sce,
                individual="biosample",
                model="NBGMM",
                NRUV=5,
                pathology="pathology",
                nebula.formula="~pathology + age + race + n_genes_by_counts",
                offset="total_counts") {
    counts = assays(sce)$counts
    mdata = as.data.frame(sce@colData)
    mdata[[pathology]] = as.factor(mdata[[pathology]])
    mdata = mdata[order(mdata[[individual]]),]
    mdata[[individual]] = as.factor(mdata[[individual]])
    data_ind = counts %*% make_pseudobulk(mdata[[individual]])
    uniq.obs = unique(mdata[,c(individual, pathology)])
    if (NRUV > 0) {
        design = model.matrix(asform(c("~", pathology)), data=uniq.obs)
        dgel = DGEList(data_ind, genes=rownames(data_ind))
        ### throw away low CPM genes
        to_keep = rowSums(edgeR::cpm(dgel) > 1) >= 3
        dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
        dgel = calcNormFactors(dgel, method="TMM")
        dgel = estimateGLMCommonDisp(dgel, design)
        dgel = estimateGLMTagwiseDisp(dgel, design)
        fit1 = glmFit(dgel, design)
        res1 = residuals(fit1, type="deviance")
        ruv_cov = RUVr(round(dgel$counts), as.character(rownames(dgel$counts)), k=NRUV, res1)

        ### Add RUVseq terms to formula
        uniq.obs = cbind(uniq.obs, ruv_cov$W)
        mdata$barcode_ = rownames(mdata)
        mdata = merge(mdata, uniq.obs, all.x=TRUE)
        rownames(mdata) = mdata$barcode_
        mdata = mdata[colnames(counts),]
        ruvw = paste0("W_", 1:NRUV)
        nebula.formula = c(nebula.formula, " + ", paste0(ruvw, collapse=" + "))
    }
    nebula.formula = asform(nebula.formula)
    print(nebula.formula)

    design = model.matrix(nebula.formula, data=mdata)
    if (is.na(offset) | !(offset %in% colnames(mdata))) {
        offset = rep(1, nrow(design))
    } else {
        print(paste0("Using offset ", offset))
        offset = mdata[rownames(design),][[offset]]
    }
    return(nebula(counts[,rownames(design)],
                  mdata[rownames(design),][[individual]],
                  design,
                  offset=offset,
                  model=model))
}

if (sys.nframe() == 0) {
    args = commandArgs(trailing=TRUE)
    if (length(args) == 3) {
        h5ad = args[1]
        ctcol = args[2]
        ct  = args[3]
        adata = read.h5ad.raw(h5ad)
        if (ctcol %in% colnames(colData(adata))) {
            print(paste0("Using cell type ", ctcol, " = ", ct))
            adata = adata[,colData(adata)[[ctcol]] == ct]
        }
        neb = deg(adata)
        saveRDS(neb, paste0("deg_", ctcol, "_", ct, ".Rds"))
        neb_df = neb$summary
        ovr = neb$overdispersion
        colnames(ovr) = paste0("overdispersion_", colnames(ovr))
        neb_df = cbind(neb_df, ovr)
        neb_df$convergence = neb$convergence
        neb_df$algorithm = neb$algorithm
        data.table::fwrite(neb_df, paste0("deg_", ctcol, "_", ct, ".tsv.gz"), sep="\t")
    }
}
