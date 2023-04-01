

library(dplyr)
library(GenomicRanges)
library(rtracklayer)

FL = list.files("~/data-group/Benjamin/AD_ATAC/modules/", pattern="*.bed.gz", full.names=T)
names(FL) = gsub(".bed.gz$", "", basename(FL))

df = dplyr::bind_rows(lapply(FL, function(x) {
    df = read.table(x);
    rownames(df) = with(df, paste0(V1, ":", V2, "-", V3));
    colnames(df) = c("chrom", "begin", "end", "module");
    df$feature = rownames(df)
    return(df)
}), .id="CellType")

df = merge(df, read.csv("/home/benjames/data-group/Benjamin/AD_ATAC/modules/feature_types.csv.gz"))
colnames(df) = gsub("^feature$", "peak", colnames(df))
gff = rtracklayer::readGFF("/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")

gr = with(gff[(gff$type == "gene") & (gff$gene_type == "protein_coding"),], GRanges(seqnames=seqid, ranges=IRanges(start, end), strand=strand, gene=make.unique(gene_name, sep="-")))
names(gr) = gr$gene
dr = with(df, GRanges(seqnames=chrom, ranges=IRanges(begin, end)))

df$nearest_gene = gr[nearest(dr, gr)]$gene
sf = df %>% group_by(CellType, module) %>% summarize(Exon=sum(feature_types == "Exon"),
                                                     Intron=sum(feature_types == "Intron"),
                                                     Distal=sum(feature_types == "Distal"),
                                                     Promoter=sum(feature_types == "Promoter"),
                                                     count=n(),
                                                     n_genes=n_distinct(nearest_gene),
                                                     n_chrom=n_distinct(chrom))

FL = list.files("~/data-group/Benjamin/AD_ATAC/linking/", pattern="^predicted.*.tsv.gz$", full.names=TRUE)
names(FL) = gsub("^predicted_", "", gsub(".tsv.gz$", "", basename(FL)))
links = dplyr::bind_rows(lapply(FL, function(x) {
    df = read.table(x, sep="\t", header=TRUE, row.names=1)
    return(df[df$score >= 0.001,])
}), .id="CellType")
lf = merge(links, df)
slf = lf %>% group_by(CellType, gene) %>% summarize(n_mod=n_distinct(module))
