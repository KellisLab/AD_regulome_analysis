
loadSignac <- function(fragment.dir, peakMatrix, sep="#", suffix="_fragments.bed.gz", sample.col="Sample") {
    require(Signac)
    require(Seurat)
    us = table(peakMatrix[[sample.col]])
    us = names(us)[us > 1]
    s.obj.list = lapply(setNames(us, us), function(sample) {
        print(paste0(sample, ", ", match(sample, us), "/", length(us)))
        fragment = paste0(fragment.dir, "/", sample, suffix)
        mat = peakMatrix[,sample==peakMatrix[[sample.col]]]
        colnames(mat) = gsub(paste0("^", sample, sep), "", colnames(mat))
        print(paste0("createchromatinassay ", nrow(mat), " ", ncol(mat)))
        chrom_assay <- CreateChromatinAssay(
            counts=assays(mat)[[1]],
            sep=c(":", "-"),
            genome="GRCh38",
            fragments=fragment)
        cd = as.data.frame(colData(mat))
        s.obj = CreateSeuratObject(
            counts=chrom_assay,
            assay="peaks",
            project=sample,
            meta.data=cd)
        return(s.obj)
    })
    return(merge(s.obj.list[[1]], s.obj.list[-1]))
}

annotations <- GetGRangesFromEnsDb(ensdb =EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) = "UCSC"
genome(annotations) = "GRCh38"
# add the gene information to the object
Annotation(s.obj) = annotations



s.obj <- RunTFIDF(s.obj)
s.obj <- FindTopFeatures(s.obj, min.cutoff = 'q0')
s.obj <- RunSVD(s.obj)
s.obj = RunUMAP(s.obj, reduction="lsi", dims=2:30)
gene.activities <- GeneActivity(s.obj)
s.obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
s.obj = NormalizeData(s.obj, assay="RNA")
DefaultAssay(s.obj) = "RNA"


r.obj = NormalizeData(r.obj)
r.obj = ScaleData(r.obj)
