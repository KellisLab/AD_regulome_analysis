#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
load.promoters <- function(gtf="/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz") {
    gtf = readGFF(gtf)
    gtf = gtf[gtf$type == "gene",]
    gtf$tss = gtf$start
    gtf$tss[gtf$strand == "-"] = gtf$end[gtf$strand == "-"]
    gtf$to_left = 2000
    gtf$to_right = 100
    gtf$to_left[gtf$strand == "-"] = 100
    gtf$to_right[gtf$strand == "-"] = 2000
    d = duplicated(gtf$gene_name)
    gtf$gene_name[d] = paste0(gtf$gene_name[d], "-1")
    gr = with(gtf, GRanges(seqnames=seqid,
                           ranges=IRanges(tss - to_left, tss + to_right),
                           gene=gene_name,
                           tss=tss,
                           strand=strand))
    names(gr) = gr$gene
    return(gr)
}

load.plac <- function(interactome_file, enhancers_file, peaks_ranges, promoters_ranges) {
    idf = read.table(interactome_file, sep="\t", header=TRUE)
    igr1 = with(idf, GRanges(seqnames=chr1, ranges=IRanges(start1, end1)))
    igr2 = with(idf, GRanges(seqnames=chr2, ranges=IRanges(start2, end2)))

    egr = with(read.table(enhancers_file, sep="\t"), GRanges(seqnames=V1, ranges=IRanges(V2, V3)))
    eovp = findOverlaps(egr, peaks_ranges)
    egr = peaks_ranges[unique(subjectHits(eovp))] ### subset peaks to overlapping enhancers
    eovp_12 = with(findOverlaps(egr, igr1),
                   data.frame(enh_idx=queryHits,
                              idf_idx=subjectHits))
    eovp_21 = with(findOverlaps(egr, igr2),
                   data.frame(enh_idx=queryHits,
                              idf_idx=subjectHits))
    govp_21 = with(findOverlaps(promoters_ranges, igr1),
                   data.frame(prom_idx=queryHits,
                              idf_idx=subjectHits))
    govp_12 = with(findOverlaps(promoters_ranges, igr2),
                   data.frame(prom_idx=queryHits,
                              idf_idx=subjectHits))
    eg = rbind(merge(eovp_12, govp_12)[c("prom_idx", "idf_idx", "enh_idx")],
               merge(eovp_21, govp_21)[c("prom_idx", "idf_idx", "enh_idx")])
    d2n = distanceToNearest(egr[eg$enh_idx], promoters_ranges[eg$prom_idx])
    eg_distance = distance(egr[eg$enh_idx], promoters_ranges[eg$prom_idx])
    closest_gene_distance = as.data.frame(d2n)$distance
    plac = with(eg, data.frame(gene=names(promoters_ranges)[prom_idx],
                               fdr=idf$fdr[idf_idx],
                               ClusterSummit=idf$ClusterSummit[idf_idx],
                               ClusterNegLog10P=idf$ClusterNegLog10P[idf_idx],
                               ClusterLabel=idf$ClusterLabel[idf_idx],
                               ClusterSize=idf$ClusterSize[idf_idx],
                               distance=eg_distance,
                               is_closest_gene=closest_gene_distance >= eg_distance,
                               peak=names(egr)[enh_idx]))
    plac = plac[order(plac$fdr),]
    plac = plac[!duplicated(plac[c("peak", "gene")]),]
    rownames(plac) = paste0(plac$peak, "|", plac$gene)
    return(plac)
}

load.allpairs <- function(genes, peaks, distance=1000000) {
    extended.gr = with(genes, GRanges(seqnames=seqnames, ranges=IRanges(tss-distance, tss+distance)))
    ovp = findOverlaps(extended.gr, peaks)
    odf = with(ovp, data.frame(gene=names(genes)[queryHits],
                               peak=names(peaks)[subjectHits]))
    odf = odf[!duplicated(odf),]
    rownames(odf) = paste0(odf$peak, "|", odf$gene)
    return(odf)
}

get.training <- function() {
    in.frame = data.frame(interactome_file=paste0("/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/Plac_seq/",
                                       c("Oligo.interactome.plaqSeq.hg38.txt",
                                         "Neuronal.interactome.plaqSeq.hg38.txt",
                                         "Neuronal.interactome.plaqSeq.hg38.txt",
                                         "Microglia.interactome.plaqSeq.hg38.txt")),
               enhancer_file=paste0("/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/Plac_seq/",
                                    c("Oligo.enhancer.hg38.txt",
                                      "Neuronal.enhancer.hg38.txt",
                                      "Neuronal.enhancer.hg38.txt",
                                      "Microglia.enhancer.hg38.txt")),
               peak_file=paste0("/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/ATAC.matrix/matrix.by.major_celltype.20220816/",
                                c("Oligo", "Ex", "In", "Microglia"),
                                "/peak.annotation.txt"))
    rownames(in.frame) = c("Oligo", "Ex", "In", "Microglia")
    promoter.gr = load.promoters("/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
    L = lapply(rownames(in.frame), function(ct) {
        peak.gr = with(read.table(in.frame[ct, "peak_file"], sep="\t", header=TRUE),
                       GRanges(seqnames=seqnames, ranges=IRanges(start, end), peak_name=paste0(seqnames, ":", start, "-", end)))
        names(peak.gr) = peak.gr$peak_name
        load.plac(interactome_file=in.frame[ct,"interactome_file"],
                  enhancers_file=in.frame[ct, "enhancer_file"],
                  peaks_ranges=peak.gr,
                  promoters_ranges=promoter.gr)
    })
    names(L) = rownames(in.frame)
    return(L)
}
get.training.microglia <- function() {
    ### peaks
    peaks = read.table("/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/ATAC.matrix/matrix.by.major_celltype.20220816/Microglia/peak.annotation.txt",sep="\t", header=T)
    rownames(peaks) = with(peaks, paste0(seqnames, ":", start, "-", end))
    peak.gr = with(peaks, GRanges(seqnames=seqnames, ranges=IRanges(start, end)))
    names(peak.gr) = rownames(peaks)
    ### promoters
    promoter.gr = load.promoters("/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
    tdf = get.training(interactome_file="/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/Plac_seq/Microglia.interactome.plaqSeq.hg38.txt",
                       enhancers_file="/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/Plac_seq/Microglia.enhancer.hg38.txt",
                       peaks_ranges=peak.gr,
                       promoters_ranges=promoter.gr)
    return(tdf)
}

load.bed <- function(bed="/net/bmc-lab5/data/kellis/group/Xushen_Sharing/ToBen/AD.GWAS/AD.GWAS.LD.tagged.merged.hg38.bed") {
    with(read.table(bed), GRanges(seqnames=V1, ranges=IRanges(V2, V3)))
}

abc.normalize <- function(df, cols=c()) {
    if (length(cols) == 0) {
        cols = grep("^score_", colnames(df))
    }
    X = as.matrix(df[,C])
    for (i in grep("^score_[0-9]+$", colnames(X))) {
        X[,i] = atanh(X[,i])
        X[,i][X[,i] > 20] = 20 ### arctanh
    }
    df$base = apply(X, 1, prod)
    sf = as.data.frame(df %>% group_by(gene) %>% summarize(base_sum=sum(base)))
    rownames(sf) = sf$gene
    df$base_sum = sf[df$gene,"base_sum"]
    df$normalized_score = df$base / df$base_sum
    return(df)
}
