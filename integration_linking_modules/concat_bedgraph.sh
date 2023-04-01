#!/usr/bin/env bash

chromsizes="/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star/chrNameLength.txt"
type -P bedtools >/dev/null || module load bedtools
type -P bedGraphToBigWig >/dev/null || module load ucsc-tools
while read celltype; do
    while read pathology; do
	echo "Pathology: ${pathology} celltype: ${celltype} summing bedgraphs"
	midput="${celltype}_${pathology}.bedgraph"
	output="${celltype}_${pathology}.bw"
	awk -v path="${pathology}" -v celltype="${celltype}" '$1 == path { print $2"_"celltype".bedgraph.gz" }' < pathology_samples.tsv \
	| xargs -n1 ls -d 2>/dev/null \
	| xargs bedtools unionbedg -i \
	| awk -v OFS="\t" '{ sum=0; for (i=4;i<=NF;i++) sum+=$i; print $1,$2,$3,sum }' \
              > "${midput}"
	echo "Pathology: ${pathology} celltype: ${celltype} writing to bigwig"
	bedGraphToBigWig "${midput}" "${chromsizes}" "${output}"
    done < <(cut -f 1 pathology_samples.tsv | sort | uniq)
done < <(cut -f 2 samples_celltypes.tsv | sort | uniq)
