#!/usr/bin/env bash
CT="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/TSS6/celltypes.csv"
GAcc="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/TSS6/GeneMatrix_estimated_0.h5ad"
RNA="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/RNA.h5ad"
GTF="/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz"
outdir="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/linking/pseudobulk/"
iter=2
while read line; do
    atac=$(echo $line | cut -d , -f 1)
    rna=$(echo $line | cut -d , -f 2)
    Peaks="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/CellTypeSpecificPeaks/${atac}.h5ad"
    integ="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/TSS6/integration/${rna}_${iter}.tsv.gz"
    if test -f "${Peaks}"; then
	echo $atac
	epiclust_pseudobulk.py -i Peaks="${Peaks}" RNA="${RNA}" \
			       -t "${integ}" -o "${outdir}/${atac}.h5ad" \
			       -c leiden_2 projid --tf-idf Peaks --gtf-annotation "${GTF}"
    fi
done  < "${CT}"
